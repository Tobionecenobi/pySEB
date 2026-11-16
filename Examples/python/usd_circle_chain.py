"""Export 100 circles linked from each center to the previous contour."""
import pyseb


world = pyseb.World("circle_chain")
graph = world.Add("ThinCircle", "circle1", "circle")

for index in range(2, 101):
    previous = f"circle{index - 1}"
    current = f"circle{index}"
    graph = world.Link(
        "ThinCircle",
        f"{current}.center",
        f"{previous}.contour#link{index - 1}",
        "circle",
    )

world.Add(graph, "chain")

parameters = {"R_circle": 0.5}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 8
options.relaxation_sweeps = 1
options.minimum_clearance = 0.1
options.curve_samples = 32
world.export_usd("chain", "circle_chain.usda", parameters, options)
