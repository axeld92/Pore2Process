using System.Diagnostics;
using System.Reflection;
using System.Runtime.Loader;
using System.Text.Json;
using System.Text.Json.Serialization;
using DWSIM.Automation.FluentAPI;
using DWSIM.Automation.FluentAPI.Builders;
using DWSIM.Interfaces.Enums;
using DWSIM.Interfaces.Enums.GraphicObjects;
using DWSIM.UnitOperations.SpecialOps;
using DWSIM.UnitOperations.UnitOperations;
using Flowsheet = DWSIM.Automation.FluentAPI.Flowsheet;

if (args.Length!=3) throw new ArgumentException("Usage: CrossflowRunner <distribution> <input.json|--inspect> <output-directory>");
var distribution=Path.GetFullPath(args[0]);
var inputPath=args[1]=="--inspect"?args[1]:Path.GetFullPath(args[1]);
var directory=Path.GetFullPath(args[2]);
AssemblyLoadContext.Default.Resolving+=(_,name)=>{
    var path=Path.Combine(distribution,name.Name+".dll");
    return File.Exists(path)?AssemblyLoadContext.Default.LoadFromAssemblyPath(path):null;
};
Directory.SetCurrentDirectory(distribution);
try { Crossflow.Run(inputPath,directory); }
catch(Exception error) { Console.Error.WriteLine(error);Environment.ExitCode=1; }

record Inputs(double Temperature_K,double External_pressure_Pa,double TMP_Pa,double Channel_drop_Pa,
    double Fresh_mass_flow_kg_s,double Module_to_permeate_ratio,double Recovery_fraction,
    double Feed_efficiency_fraction,double Circulation_efficiency_fraction,
    double Seed_module_mass_flow_kg_s,double Seed_temperature_K,
    double Mass_tolerance_kg_s,double Temperature_tolerance_K,double Pressure_tolerance_Pa,
    int Maximum_iterations,string Scope,
    [property:JsonIgnore(Condition=JsonIgnoreCondition.WhenWritingDefault)] bool Force_liquid=false);
record State(double pressure_Pa,double temperature_K,double mass_flow_kg_s,
    double actual_volume_flow_m3_s,double density_kg_m3,double viscosity_Pa_s,
    double enthalpy_J_kg,double vapor_fraction,double water_mass_fraction);

static class Crossflow {
    static readonly JsonSerializerOptions Json=new(){WriteIndented=true};
    static State Read(MaterialStreamBuilder s)=>new(s.PressurePa,s.TemperatureK,s.MassFlowKgPerSecond,
        s.VolumetricFlowM3PerSecond,s.Object.Phases[0].Properties.density??throw new Exception("Missing density"),
        s.Object.Phases[1].Properties.viscosity??throw new Exception("Missing viscosity"),
        1000*(s.Object.Phases[0].Properties.enthalpy??throw new Exception("Missing enthalpy")),
        s.VaporFraction,s.OverallMassFraction("Water"));
    static void Solve(Flowsheet fs,Recycle rec){
        try { fs.Solve(); }
        catch(Exception e) {
            Console.WriteLine(JsonSerializer.Serialize(Convergence(rec),Json));
            Console.WriteLine(JsonSerializer.Serialize(new {from=Read(fs.MaterialStream("Circulation-in")),to=Read(fs.MaterialStream("Tear-in"))},Json));
            throw new InvalidOperationException("RECYCLE_NOT_CONVERGED: DWSIM solve failed",e);
        }
        if (!rec.Converged) throw new InvalidOperationException("RECYCLE_NOT_CONVERGED: stock block did not converge");
    }
    static object Convergence(Recycle r)=>new {
        converged=r.Converged,iterations_taken=r.IterationsTaken,iteration_count=r.IterationCount,
        acceleration=r.AccelerationMethod.ToString(),legacy_mode=r.LegacyMode,maximum_iterations=r.MaximumIterations,
        temperature_tolerance_K=r.ConvergenceParameters.Temperatura,
        pressure_tolerance_Pa=r.ConvergenceParameters.Pressao,mass_tolerance_kg_s=r.ConvergenceParameters.VazaoMassica,
        history_temperature_error_K=r.ConvergenceHistory.TemperaturaE,
        history_pressure_error_Pa=r.ConvergenceHistory.PressaoE,
        history_mass_error_kg_s=r.ConvergenceHistory.VazaoMassicaE,
        errors=r.Errors.ToDictionary(v=>v.Key,v=>v.Value),values=r.Values.ToDictionary(v=>v.Key,v=>v.Value)
    };
    public static void Run(string inputPath,string directory){
        Directory.CreateDirectory(directory);
        if(inputPath=="--inspect"){
            Type[] ts=[typeof(Flowsheet),typeof(MixerBuilder),typeof(CoolerBuilder),typeof(GenericUnitOpBuilder),
                typeof(Recycle),typeof(DWSIM.UnitOperations.SpecialOps.Helpers.Recycle.ConvergenceParameters)];
            var inventory=new {version=typeof(Flowsheet).Assembly.GetName().Version?.ToString(),types=ts.ToDictionary(
                t=>t.FullName!,t=>new {methods=t.GetMethods(BindingFlags.Public|BindingFlags.Instance|BindingFlags.DeclaredOnly)
                    .Where(m=>!m.IsSpecialName).Select(m=>m.ToString()).ToArray(),
                    properties=t.GetProperties(BindingFlags.Public|BindingFlags.Instance|BindingFlags.DeclaredOnly).Select(p=>p.ToString()).ToArray(),
                    fields=t.GetFields(BindingFlags.Public|BindingFlags.Instance|BindingFlags.DeclaredOnly).Select(f=>f.ToString()).ToArray()})};
            File.WriteAllText(Path.Combine(directory,"installed_api.json"),JsonSerializer.Serialize(inventory,Json)+"\n");return;
        }
        var a=JsonSerializer.Deserialize<Inputs>(File.ReadAllText(inputPath))??throw new ArgumentException("Missing inputs");
        double[] positives=[a.Temperature_K,a.External_pressure_Pa,a.TMP_Pa,a.Channel_drop_Pa,a.Fresh_mass_flow_kg_s,
            a.Module_to_permeate_ratio,a.Recovery_fraction,a.Feed_efficiency_fraction,a.Circulation_efficiency_fraction,
            a.Seed_module_mass_flow_kg_s,a.Seed_temperature_K,a.Mass_tolerance_kg_s,a.Temperature_tolerance_K,a.Pressure_tolerance_Pa];
        if(positives.Any(v=>!double.IsFinite(v)||v<=0)||a.Recovery_fraction>=1||a.Module_to_permeate_ratio<=1/a.Recovery_fraction||
            a.Feed_efficiency_fraction>1||a.Circulation_efficiency_fraction>1||a.Maximum_iterations<1||a.Maximum_iterations>20000||
            !double.IsFinite(a.External_pressure_Pa+a.TMP_Pa+a.Channel_drop_Pa))
            throw new ArgumentException("Require finite positive inputs, recovery in (0,1), module/permeate > 1/recovery, efficiencies <=1 and 1..20000 iterations.");
        var watch=Stopwatch.StartNew();
        var fs=Flowsheet.Create("Pore2Process: clean-water crossflow recycle")
            .WithCompound("Water").WithPropertyPackage(PropertyPackages.SteamTables);
        var high=a.External_pressure_Pa+a.TMP_Pa;
        var fresh=fs.AddMaterialStream("Fresh").At(Q.Kelvin(a.Temperature_K),Q.Pascal(a.External_pressure_Pa))
            .WithMassFlow(Q.KgPerSecond(a.Fresh_mass_flow_kg_s)).PositionAt(20,220);
        var freshHigh=fs.AddMaterialStream("Fresh-high").PositionAt(200,220);
        var mixed=fs.AddMaterialStream("Mixed").PositionAt(380,220);
        var circIn=fs.AddMaterialStream("Circulation-in").At(Q.Kelvin(a.Seed_temperature_K),Q.Pascal(high))
            .WithMassFlow(Q.KgPerSecond(a.Seed_module_mass_flow_kg_s)).PositionAt(700,220);
        var circOut=fs.AddMaterialStream("Circulation-out").PositionAt(860,220);
        var membraneFeed=fs.AddMaterialStream("Membrane-feed").PositionAt(1020,220);
        var permeateHigh=fs.AddMaterialStream("Permeate-high").PositionAt(1190,100);
        var permeate=fs.AddMaterialStream("Permeate").PositionAt(1370,100);
        var retentate=fs.AddMaterialStream("Module-retentate").PositionAt(1190,340);
        var purgeHigh=fs.AddMaterialStream("Purge-high").PositionAt(1370,340);
        var purge=fs.AddMaterialStream("Purge").PositionAt(1530,340);
        var toTear=fs.AddMaterialStream("Tear-in").PositionAt(540,220);
        // Respect the specified TP state at the temperature-controlled tear.
        // DWSIM otherwise re-flashes connected pure-compound streams as PH.
        toTear.Object.OverrideSingleCompoundFlashBehavior=true;
        circIn.Object.OverrideSingleCompoundFlashBehavior=true;
        var retReturn=fs.AddMaterialStream("Return").PositionAt(1050,510).FlipHorizontal(true);
        var feedPower=fs.AddEnergyStream("Feed-power").PositionAt(90,350);
        var circPower=fs.AddEnergyStream("Circulation-power").PositionAt(760,110);
        var heat=fs.AddEnergyStream("Cooling-duty").PositionAt(450,350);
        var feedPump=fs.AddPump("P-feed").WithOutletPressure(Q.Pascal(high))
            .WithEfficiencyPercent(100*a.Feed_efficiency_fraction).ConnectFeed(fresh).ConnectProduct(freshHigh)
            .ConnectEnergyFeed(feedPower,1).PositionAt(115,220);
        fs.AddMixer("M-1").WithPressureBehavior(Mixer.PressureBehavior.Minimum)
            .ConnectFeed(freshHigh,0).ConnectFeed(retReturn,1).ConnectProduct(mixed).PositionAt(300,220);
        var cooler=fs.AddCooler("C-1").WithOutletTemperature(Q.Kelvin(a.Temperature_K)).WithPressureDrop(Q.Pascal(0))
            .WithEfficiencyPercent(100).ConnectFeed(mixed).ConnectProduct(toTear).ConnectEnergyProduct(heat,0).PositionAt(460,220);
        var circulationPump=fs.AddPump("P-circulation").WithOutletPressure(Q.Pascal(high+a.Channel_drop_Pa))
            .WithEfficiencyPercent(100*a.Circulation_efficiency_fraction).ConnectFeed(circIn).ConnectProduct(circOut)
            .ConnectEnergyFeed(circPower,1).PositionAt(780,220);
        fs.AddValve("V-channel").WithPressureDrop(Q.Pascal(a.Channel_drop_Pa))
            .ConnectFeed(circOut).ConnectProduct(membraneFeed).PositionAt(940,220);
        var alpha=1/a.Module_to_permeate_ratio;
        fs.AddSplitter("S-membrane").WithSplitRatios(alpha,1-alpha).ConnectFeed(membraneFeed)
            .ConnectProduct(permeateHigh,0).ConnectProduct(retentate,1).PositionAt(1110,220);
        fs.AddValve("V-membrane").WithPressureDrop(Q.Pascal(a.TMP_Pa))
            .ConnectFeed(permeateHigh).ConnectProduct(permeate).PositionAt(1290,100);
        var purgeFraction=(1/a.Recovery_fraction-1)/(a.Module_to_permeate_ratio-1);
        fs.AddSplitter("S-purge").WithSplitRatios(purgeFraction,1-purgeFraction).ConnectFeed(retentate)
            .ConnectProduct(purgeHigh,0).ConnectProduct(retReturn,1).PositionAt(1290,340);
        fs.AddValve("V-purge").WithPressureDrop(Q.Pascal(a.TMP_Pa))
            .ConnectFeed(purgeHigh).ConnectProduct(purge).PositionAt(1450,340);
        var block=fs.AddUnitOperation(ObjectType.OT_Recycle,"REC-1").ConnectFeed(toTear).ConnectProduct(circIn)
            .PositionAt(620,220);
        var recycle=(Recycle)block.Object;
        recycle.AccelerationMethod=AccelMethod.None;recycle.LegacyMode=true;recycle.MaximumIterations=a.Maximum_iterations;
        recycle.ConvergenceParameters.Temperatura=a.Temperature_tolerance_K;
        recycle.ConvergenceParameters.Pressao=a.Pressure_tolerance_Pa;
        recycle.ConvergenceParameters.VazaoMassica=a.Mass_tolerance_kg_s;
        recycle.ConvergenceParameters.VazaoMassicaRelativa=1e-12;
        MaterialStreamBuilder[] material=[fresh,freshHigh,mixed,circIn,circOut,membraneFeed,permeateHigh,permeate,retentate,purgeHigh,purge,toTear,retReturn];
        // Optional stock single-liquid flash route for the low-head showcase.
        // Keep older examples on their original default phase calculation.
        if(a.Force_liquid) foreach(var stream in material) stream.Object.ForcePhase=ForcedPhase.Liquid;
        var initialization=watch.Elapsed.TotalSeconds;watch.Restart();
        Solve(fs,recycle);var solve=watch.Elapsed.TotalSeconds;
        var convergence=Convergence(recycle);
        var streams=material.ToDictionary(s=>s.Object.GraphicObject.Tag,s=>Read(s));
        var saturationPressures=material.ToDictionary(s=>s.Object.GraphicObject.Tag,s=>
            ((DWSIM.Thermodynamics.PropertyPackages.PropertyPackage)s.Object.PropertyPackage).AUX_PVAPi(0,s.TemperatureK));
        if(a.Force_liquid && material.Any(s=>s.TemperatureK<=273.15 || s.PressurePa<=2*saturationPressures[s.Object.GraphicObject.Tag]))
            throw new InvalidOperationException("Forced liquid route requires water above freezing and pressure exceeding twice its stock saturation pressure.");
        var powers=new {feed_W=1000*feedPump.PowerKW,circulation_W=1000*circulationPump.PowerKW,
            feed_energy_stream_W=1000*feedPower.EnergyFlowKW,circulation_energy_stream_W=1000*circPower.EnergyFlowKW,
            cooling_removed_W=1000*cooler.HeatRemovedKW,cooling_energy_stream_W=1000*heat.EnergyFlowKW};
        var units=fs.Inner.SimulationObjects.Values.ToDictionary(o=>o.GraphicObject.Tag,o=>o.GetType().FullName);
        fs.AddText("Scope",$"Clean-water crossflow hydraulic surrogate | recovery {100*a.Recovery_fraction:G3}% | module/permeate {a.Module_to_permeate_ratio:G3}\n"+
            $"Julia TMP {a.TMP_Pa:G7} Pa | channel loss {a.Channel_drop_Pa:G7} Pa | stock recycle {recycle.IterationsTaken} iterations\n"+
            $"Fresh-feed pump {powers.feed_W:G7} W | circulation pump {powers.circulation_W:G7} W | cooling {powers.cooling_removed_W:G7} W",
            20,10);
        var path=Path.Combine(directory,"crossflow.dwxmz");fs.Save(path,true);
        var loaded=Flowsheet.Load(path);
        var rec2=(Recycle)loaded.Inner.SimulationObjects.Values.Single(o=>o.GraphicObject.Tag=="REC-1");
        Solve(loaded,rec2);
        var after=streams.Keys.ToDictionary(n=>n,n=>Read(loaded.MaterialStream(n)));
        var energies=new[]{"Feed-power","Circulation-power","Cooling-duty"}.ToDictionary(n=>n,n=>1000*loaded.EnergyStream(n).EnergyFlowKW);
        fs.SaveScreenshot(Path.Combine(directory,"flowsheet.png"));
        var result=new {tool="DWSIM",version=typeof(Flowsheet).Assembly.GetName().Version?.ToString(),timestamp_utc=DateTimeOffset.UtcNow,
            inputs=a,property_package=PropertyPackages.SteamTables,compounds=new[]{"Water"},streams,powers,units,
            forced_phases=material.ToDictionary(s=>s.Object.GraphicObject.Tag,s=>s.Object.ForcePhase.ToString()),
            reloaded_forced_phases=streams.Keys.ToDictionary(n=>n,n=>loaded.MaterialStream(n).Object.ForcePhase.ToString()),
            saturation_pressure_Pa=saturationPressures,
            recycle=convergence,reloaded_recycle=Convergence(rec2),reloaded_streams=after,reloaded_energy_streams_W=energies,
            tear_single_compound_flash_override=toTear.Object.OverrideSingleCompoundFlashBehavior&&circIn.Object.OverrideSingleCompoundFlashBehavior,
            save_reload_solve_checked=true,initialization_s=initialization,solve_s=solve,
            contract="Suffixes are SI. DWSIM kW/kJ per kg are converted once to W/J per kg. Recovery is mass based; actual volumes differ from fresh-density reference volumes. Cooling duty is heat removed, not electricity."};
        File.WriteAllText(Path.Combine(directory,"summary.json"),JsonSerializer.Serialize(result,Json)+"\n");
        Console.WriteLine($"Converged/reloaded: feed {powers.feed_W:G12} W, circulation {powers.circulation_W:G12} W, cooling {powers.cooling_removed_W:G12} W; iterations {recycle.IterationsTaken}");
    }
}
