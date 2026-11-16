"""Export 20 cylinders linked end-to-end through distributed end caps."""
import pyseb


world = pyseb.World("cylinder_chain")
graph = world.Add("SolidCylinder", "cylinder1", "cylinder")

for index in range(2, 21):
    previous = f"cylinder{index - 1}"
    current = f"cylinder{index}"
    graph = world.Link(
        "SolidCylinder",
        f"{current}.ends#bottom#joint{index - 1}",
        f"{previous}.ends#top#joint{index - 1}",
        "cylinder",
    )

world.Add(graph, "chain")

parameters = {"R_cylinder": 0.25, "L_cylinder": 3.0}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 16
options.relaxation_sweeps = 4
options.minimum_clearance = 0.15
options.surface_samples = 16
world.export_usd("chain", "cylinder_chain.usda", parameters, options)
