"""Export 100 spheres linked from each center to the previous surface."""
import pyseb


world = pyseb.World("sphere_chain")
graph = world.Add("SolidSphere", "sphere1", "sphere")

for index in range(2, 101):
    previous = f"sphere{index - 1}"
    current = f"sphere{index}"
    graph = world.Link(
        "SolidSphere",
        f"{current}.center",
        f"{previous}.surface#link{index - 1}",
        "sphere",
    )

world.Add(graph, "chain")

parameters = {"R_sphere": 0.5}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 8
options.relaxation_sweeps = 1
options.minimum_clearance = 0.1
options.surface_samples = 16
world.export_usd("chain", "sphere_chain.usda", parameters, options)
