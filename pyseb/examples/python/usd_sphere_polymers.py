"""Export a sphere with 100 short, freely oriented surface-linked polymers."""
import pyseb

world = pyseb.World("sphere_polymers")
graph = world.Add("SolidSphere", "aSphere")
for index in range(1, 101):
    graph = world.Link(
        "GaussianPolymer",
        f"polymer{index}.end1",
        f"aSphere.surface#polymer{index}",
    )
world.Add(graph, "assembly")

parameters = {"R_aSphere": 3.0}
parameters.update({f"Rg_polymer{index}": 0.2 for index in range(1, 101)})
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.surface_samples = 48
options.curve_samples = 96
world.export_usd("assembly", "sphere_polymers.usda", parameters, options)
