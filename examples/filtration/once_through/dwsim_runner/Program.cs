using System.Diagnostics;
using System.Reflection;
using System.Runtime.Loader;
using System.Text.Json;
using DWSIM.Automation.FluentAPI;
using DWSIM.Automation.FluentAPI.Builders;

if (args.Length != 3) throw new ArgumentException("Usage: FiltrationRunner <distribution> <input.json|--inspect> <output-directory>");
var distribution = Path.GetFullPath(args[0]);
var inputPath = args[1] == "--inspect" ? args[1] : Path.GetFullPath(args[1]);
var outputDirectory = Path.GetFullPath(args[2]);
AssemblyLoadContext.Default.Resolving += (_, name) => {
    var path = Path.Combine(distribution, name.Name + ".dll");
    return File.Exists(path) ? AssemblyLoadContext.Default.LoadFromAssemblyPath(path) : null;
};
Directory.SetCurrentDirectory(distribution);
Filtration.Run(inputPath, outputDirectory);

record Inputs(double Temperature_K, double External_pressure_Pa, double TMP_Pa,
              double Channel_drop_Pa, double Feed_mass_flow_kg_s, double Recovery_fraction,
              double Pump_efficiency_fraction, string Scope);
record StreamState(double pressure_Pa, double temperature_K, double mass_flow_kg_s,
                   double actual_volume_flow_m3_s, double density_kg_m3,
                   double viscosity_Pa_s, double enthalpy_J_kg, double vapor_fraction,
                   double water_mass_fraction);

static class Filtration {
    static readonly JsonSerializerOptions Json = new() { WriteIndented = true };
    static StreamState State(MaterialStreamBuilder s) => new(
        s.PressurePa, s.TemperatureK, s.MassFlowKgPerSecond, s.VolumetricFlowM3PerSecond,
        s.Object.Phases[0].Properties.density ?? throw new Exception("Missing density"),
        s.Object.Phases[1].Properties.viscosity ?? throw new Exception("Missing liquid viscosity"),
        1000*(s.Object.Phases[0].Properties.enthalpy ?? throw new Exception("Missing enthalpy")),
        s.VaporFraction, s.OverallMassFraction("Water"));
    static double EnthalpyFlow(StreamState s) => s.mass_flow_kg_s*s.enthalpy_J_kg;

    public static void Run(string inputPath, string directory) {
        Directory.CreateDirectory(directory);
        if (inputPath == "--inspect") {
            var a=typeof(Flowsheet).Assembly;
            var names=new[]{"Flowsheet","SplitterBuilder","ValveBuilder","MixerBuilder","CrossflowUFBuilder"};
            var inventory=new { version=a.GetName().Version?.ToString(),
                types=a.GetExportedTypes().Where(t=>names.Contains(t.Name)).ToDictionary(
                    t=>t.FullName!, t=>t.GetMethods(BindingFlags.Public|BindingFlags.Instance|BindingFlags.DeclaredOnly)
                        .Where(m=>!m.IsSpecialName).Select(m=>m.ToString()).ToArray()) };
            File.WriteAllText(Path.Combine(directory,"installed_api.json"),JsonSerializer.Serialize(inventory,Json)+"\n");
            return;
        }
        var input=JsonSerializer.Deserialize<Inputs>(File.ReadAllText(inputPath)) ?? throw new ArgumentException("Missing inputs");
        double[] v=[input.Temperature_K,input.External_pressure_Pa,input.TMP_Pa,input.Channel_drop_Pa,
                    input.Feed_mass_flow_kg_s,input.Recovery_fraction,input.Pump_efficiency_fraction];
        if (v.Any(x=>!double.IsFinite(x)) || input.Temperature_K<=0 || input.External_pressure_Pa<=0 ||
            input.TMP_Pa<=0 || input.Channel_drop_Pa<0 || input.Feed_mass_flow_kg_s<=0 ||
            input.Recovery_fraction<=0 || input.Recovery_fraction>=1 ||
            input.Pump_efficiency_fraction<=0 || input.Pump_efficiency_fraction>1 ||
            !double.IsFinite(input.External_pressure_Pa+input.TMP_Pa+input.Channel_drop_Pa))
            throw new ArgumentException("Require finite positive T, external P, TMP and feed mass flow; channel drop >=0, recovery in (0,1), efficiency in (0,1].");
        var timer=Stopwatch.StartNew();
        var fs=Flowsheet.Create("Pore2Process: once-through hydraulic filtration")
            .WithCompound("Water").WithPropertyPackage(PropertyPackages.SteamTables);
        var feed=fs.AddMaterialStream("Feed").At(Q.Kelvin(input.Temperature_K),Q.Pascal(input.External_pressure_Pa))
            .WithMassFlow(Q.KgPerSecond(input.Feed_mass_flow_kg_s)).PositionAt(40,170);
        var discharge=fs.AddMaterialStream("Pump-discharge").PositionAt(200,170);
        var membraneFeed=fs.AddMaterialStream("Membrane-feed").PositionAt(380,170);
        var highP=fs.AddMaterialStream("Permeate-high").PositionAt(580,90);
        var highR=fs.AddMaterialStream("Retentate-high").PositionAt(580,280);
        var permeate=fs.AddMaterialStream("Permeate").PositionAt(800,90);
        var retentate=fs.AddMaterialStream("Retentate").PositionAt(800,280);
        var energy=fs.AddEnergyStream("Pump-power").PositionAt(115,290);
        var pump=fs.AddPump("P-1").WithOutletPressure(Q.Pascal(input.External_pressure_Pa+input.TMP_Pa+input.Channel_drop_Pa))
            .WithEfficiencyPercent(100*input.Pump_efficiency_fraction)
            .ConnectFeed(feed).ConnectProduct(discharge).ConnectEnergyFeed(energy,1).PositionAt(130,170);
        fs.AddValve("V-channel").WithPressureDrop(Q.Pascal(input.Channel_drop_Pa))
            .ConnectFeed(discharge).ConnectProduct(membraneFeed).PositionAt(295,170);
        fs.AddSplitter("S-1").WithSplitRatios(input.Recovery_fraction,1-input.Recovery_fraction)
            .ConnectFeed(membraneFeed).ConnectProduct(highP,0).ConnectProduct(highR,1).PositionAt(490,170);
        // DeltaP mode preserves the Julia-supplied hydraulic load; no Kv law is substituted.
        fs.AddValve("V-membrane").WithPressureDrop(Q.Pascal(input.TMP_Pa))
            .ConnectFeed(highP).ConnectProduct(permeate).PositionAt(690,90);
        fs.AddValve("V-retentate").WithPressureDrop(Q.Pascal(input.TMP_Pa))
            .ConnectFeed(highR).ConnectProduct(retentate).PositionAt(690,280);
        var initialization=timer.Elapsed.TotalSeconds;timer.Restart();
        fs.Solve();var solve=timer.Elapsed.TotalSeconds;
        var streams=new Dictionary<string,StreamState> {
            ["Feed"]=State(feed),["Pump-discharge"]=State(discharge),["Membrane-feed"]=State(membraneFeed),
            ["Permeate-high"]=State(highP),["Retentate-high"]=State(highR),
            ["Permeate"]=State(permeate),["Retentate"]=State(retentate)};
        double mass(params string[] names)=>names.Sum(n=>streams[n].mass_flow_kg_s);
        double h(params string[] names)=>names.Sum(n=>EnthalpyFlow(streams[n]));
        var power=1000*pump.PowerKW;
        var residuals=new Dictionary<string,double> {
            ["pump"]=h("Pump-discharge")-h("Feed")-power,
            ["channel"]=h("Membrane-feed")-h("Pump-discharge"),
            ["splitter"]=h("Permeate-high","Retentate-high")-h("Membrane-feed"),
            ["membrane"]=h("Permeate")-h("Permeate-high"),
            ["retentate_valve"]=h("Retentate")-h("Retentate-high"),
            ["overall"]=h("Permeate","Retentate")-h("Feed")-power};
        var massResiduals=new Dictionary<string,double> {
            ["pump"]=mass("Pump-discharge")-mass("Feed"),
            ["channel"]=mass("Membrane-feed")-mass("Pump-discharge"),
            ["splitter"]=mass("Permeate-high","Retentate-high")-mass("Membrane-feed"),
            ["membrane"]=mass("Permeate")-mass("Permeate-high"),
            ["retentate_valve"]=mass("Retentate")-mass("Retentate-high"),
            ["overall"]=mass("Permeate","Retentate")-mass("Feed")};
        var objects=fs.Inner.SimulationObjects.Values.ToDictionary(o=>o.GraphicObject.Tag,o=>o.GetType().FullName);
        fs.AddText("Hydraulic-scope",$"Once-through clean-water hydraulic surrogate | recovery {100*input.Recovery_fraction:G4}%\n"+
            $"Julia TMP {input.TMP_Pa:G7} Pa | assumed channel loss {input.Channel_drop_Pa:G7} Pa | DWSIM shaft power {power:G7} W",
            40,20);
        var path=Path.Combine(directory,"filtration.dwxmz");fs.Save(path,true);
        var reloaded=Flowsheet.Load(path);reloaded.Solve();
        var after=streams.Keys.ToDictionary(n=>n,n=>State(reloaded.MaterialStream(n)));
        foreach(var (name,before) in streams) {
            var s=after[name];
            if (Math.Abs(s.pressure_Pa-before.pressure_Pa)>1e-6 ||
                Math.Abs(s.mass_flow_kg_s-before.mass_flow_kg_s)>1e-10*input.Feed_mass_flow_kg_s ||
                Math.Abs(s.temperature_K-before.temperature_K)>1e-7 ||
                Math.Abs(s.enthalpy_J_kg-before.enthalpy_J_kg)>1e-3)
                throw new Exception("Reloaded stream changed: "+name);
        }
        var reloadPower=1000*reloaded.EnergyStream("Pump-power").EnergyFlowKW;
        if (Math.Abs(reloadPower-power)>1e-9*power) throw new Exception("Reloaded pump power changed");
        // The screenshot is a genuine render of the solved DWSIM flowsheet.
        fs.SaveScreenshot(Path.Combine(directory,"flowsheet.png"));
        var result=new {
            tool="DWSIM",version=typeof(Flowsheet).Assembly.GetName().Version?.ToString(),timestamp_utc=DateTimeOffset.UtcNow,
            inputs=input,property_package=PropertyPackages.SteamTables,compounds=new[]{"Water"},streams,objects,
            shaft_power_W=power,energy_stream_power_W=1000*energy.EnergyFlowKW,
            mass_residuals_kg_s=massResiduals,energy_residuals_W=residuals,
            reloaded_streams=after,reloaded_energy_stream_power_W=reloadPower,save_reload_solve_checked=true,
            initialization_s=initialization,solve_s=solve,
            units="JSON suffixes are SI; DWSIM kW and kJ/kg converted once to W and J/kg. Actual volume flow is flashed; recovery uses mass/common reference density."
        };
        File.WriteAllText(Path.Combine(directory,"summary.json"),JsonSerializer.Serialize(result,Json)+"\n");
        Console.WriteLine($"Solved and reloaded filtration: {power:G12} W, recovery {input.Recovery_fraction:G3}");
    }
}
