"""Export a sphere with three deterministic surface-attached rods."""
import pyseb

world = pyseb.World("sphere_rods")
graph = world.Add("SolidSphere", "aSphere")
for index, label in enumerate(("front0", "front18", "front29"), 1):
    graph = world.Link(
        "ThinRod",
        f"rod{index}.end1",
        f"aSphere.surface#{label}",
    )
world.Add(graph, "assembly")

parameters = {
    "R_aSphere": 3.0,
    "L_rod1": 3.0,
    "L_rod2": 3.0,
    "L_rod3": 3.0,
}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.surface_samples = 48
options.curve_samples = 2
world.export_usd("assembly", "sphere_rods.usda", parameters, options)
