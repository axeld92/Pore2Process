using System.Diagnostics;
using System.Runtime.Loader;
using System.Text.Json;
using DWSIM.Automation.FluentAPI;

// Install the resolver before JIT-compiling the method that references DWSIM.
if (args.Length != 3) throw new ArgumentException("Usage: DwsimProbe <distribution> <input.json> <output-directory>");
var distribution = Path.GetFullPath(args[0]);
var inputPath = Path.GetFullPath(args[1]);
var outputDirectory = Path.GetFullPath(args[2]);
AssemblyLoadContext.Default.Resolving += (_, name) => {
    var path = Path.Combine(distribution, name.Name + ".dll");
    return File.Exists(path) ? AssemblyLoadContext.Default.LoadFromAssemblyPath(path) : null;
};
Directory.SetCurrentDirectory(distribution);
Probe.Run(inputPath, outputDirectory);

record Inputs(double Temperature_K, double Inlet_pressure_Pa, double Outlet_pressure_Pa,
              double Mass_flow_kg_s, double Pump_efficiency_fraction,
              string Scope = "Independent water-pump compatibility probe; no OpenLB-derived input yet.");

static class Probe {
    public static void Run(string inputPath, string outputDirectory) {
        var input = JsonSerializer.Deserialize<Inputs>(File.ReadAllText(inputPath))
                    ?? throw new ArgumentException("Missing input snapshot.");
        var values = new [] { input.Temperature_K, input.Inlet_pressure_Pa, input.Outlet_pressure_Pa,
                             input.Mass_flow_kg_s, input.Pump_efficiency_fraction };
        if (values.Any(v => !double.IsFinite(v)) || input.Temperature_K <= 0 || input.Inlet_pressure_Pa <= 0 ||
            input.Outlet_pressure_Pa <= input.Inlet_pressure_Pa || input.Mass_flow_kg_s <= 0 ||
            input.Pump_efficiency_fraction <= 0 || input.Pump_efficiency_fraction > 1)
            throw new ArgumentException("Require finite, positive temperature/pressures/mass flow, outlet > inlet pressure, and efficiency in (0,1].");
        Directory.CreateDirectory(outputDirectory);
        var timer = Stopwatch.StartNew();
        var fs = Flowsheet.Create("Pore2Process compatibility: water pump")
            .WithCompound("Water").WithPropertyPackage(PropertyPackages.SteamTables);
        var feed = fs.AddMaterialStream("Water-in")
            .At(Q.Kelvin(input.Temperature_K), Q.Pascal(input.Inlet_pressure_Pa))
            .WithMassFlow(Q.KgPerSecond(input.Mass_flow_kg_s));
        var product = fs.AddMaterialStream("Water-out");
        var energy = fs.AddEnergyStream("Pump-power");
        var pump = fs.AddPump("P-1")
            .WithOutletPressure(Q.Pascal(input.Outlet_pressure_Pa))
            .WithEfficiencyPercent(100 * input.Pump_efficiency_fraction)
            .ConnectFeed(feed).ConnectProduct(product).ConnectEnergyFeed(energy, 1);
        fs.AutoLayout();
        var initializationSeconds = timer.Elapsed.TotalSeconds;
        timer.Restart();
        fs.Solve(); // Throws if the real DWSIM flowsheet solver returns errors.
        var solveSeconds = timer.Elapsed.TotalSeconds;
        var density = feed.Object.Phases[0].Properties.density
                      ?? throw new InvalidOperationException("No feed density returned.");
        var viscosity = feed.Object.Phases[1].Properties.viscosity
                      ?? throw new InvalidOperationException("No liquid dynamic viscosity returned.");
        var outletViscosity = product.Object.Phases[1].Properties.viscosity
                      ?? throw new InvalidOperationException("No outlet liquid viscosity returned.");
        var inletEnthalpy = feed.Object.Phases[0].Properties.enthalpy
                      ?? throw new InvalidOperationException("No feed enthalpy returned.");
        var outletEnthalpy = product.Object.Phases[0].Properties.enthalpy
                      ?? throw new InvalidOperationException("No product enthalpy returned.");
        // DWSIM stores power in kW and mass-specific enthalpy in kJ/kg.
        var powerW = 1000 * pump.PowerKW;
        var expectedPowerW = input.Mass_flow_kg_s * (input.Outlet_pressure_Pa - input.Inlet_pressure_Pa)
                             / (density * input.Pump_efficiency_fraction);
        var result = new {
            tool = "DWSIM", version = typeof(Flowsheet).Assembly.GetName().Version?.ToString(),
            timestamp_utc = DateTimeOffset.UtcNow, inputs = input,
            property_package = PropertyPackages.SteamTables, compounds = new [] {"Water"},
            inlet_pressure_Pa = feed.PressurePa, outlet_pressure_Pa = product.PressurePa,
            inlet_temperature_K = feed.TemperatureK, outlet_temperature_K = product.TemperatureK,
            inlet_mass_flow_kg_s = feed.MassFlowKgPerSecond,
            outlet_mass_flow_kg_s = product.MassFlowKgPerSecond,
            inlet_density_kg_m3 = density, inlet_enthalpy_J_kg = 1000 * inletEnthalpy,
            inlet_viscosity_Pa_s = viscosity, outlet_viscosity_Pa_s = outletViscosity,
            outlet_density_kg_m3 = product.Object.Phases[0].Properties.density,
            outlet_enthalpy_J_kg = 1000 * outletEnthalpy,
            inlet_vapor_fraction = feed.VaporFraction, outlet_vapor_fraction = product.VaporFraction,
            shaft_power_W = powerW, energy_stream_power_W = 1000 * energy.EnergyFlowKW,
            analytical_shaft_power_W = expectedPowerW,
            mass_balance_error_kg_s = product.MassFlowKgPerSecond - feed.MassFlowKgPerSecond,
            energy_balance_error_W = product.MassFlowKgPerSecond * 1000 * outletEnthalpy
                                  - feed.MassFlowKgPerSecond * 1000 * inletEnthalpy - powerW,
            pump_power_relative_error = Math.Abs(powerW - expectedPowerW) / expectedPowerW,
            save_reload_solve_checked = true,
            initialization_s = initializationSeconds, solve_s = solveSeconds,
            scope = input.Scope
        };
        var flowsheetPath = Path.Combine(outputDirectory, "water_pump.dwxmz");
        fs.Save(flowsheetPath, true);
        var reloaded = Flowsheet.Load(flowsheetPath);
        reloaded.Solve();
        var reloadedProduct = reloaded.MaterialStream("Water-out");
        if (Math.Abs(reloadedProduct.PressurePa - product.PressurePa) > 1e-6 ||
            Math.Abs(reloadedProduct.MassFlowKgPerSecond - product.MassFlowKgPerSecond) > 1e-9)
            throw new InvalidOperationException("Saved/reloaded flowsheet pressure or mass flow changed.");
        File.WriteAllText(Path.Combine(outputDirectory, "summary.json"),
            JsonSerializer.Serialize(result, new JsonSerializerOptions { WriteIndented = true }) + "\n");
        Console.WriteLine($"DWSIM solved: {powerW:F6} W; outlet {product.PressurePa:F1} Pa; {product.TemperatureK:F6} K");
    }
}
