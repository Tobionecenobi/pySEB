"""Export 20 spheroids linked from north pole to south pole."""
import pyseb


world = pyseb.World("spheroid_chain")
graph = world.Add("Spheroid", "spheroid1", "spheroid")
for index in range(2, 21):
    graph = world.Link(
        "Spheroid",
        f"spheroid{index}.south",
        f"spheroid{index - 1}.north",
        "spheroid",
    )
world.Add(graph, "chain")

parameters = {"a_spheroid": 0.4, "c_spheroid": 1.2}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 32
options.relaxation_sweeps = 6
options.minimum_clearance = 0.15
options.surface_samples = 24
colors = (
    (0.12, 0.38, 0.72),
    (0.10, 0.62, 0.66),
    (0.32, 0.68, 0.38),
    (0.91, 0.63, 0.13),
    (0.82, 0.25, 0.18),
)
options.color_overrides = {
    f"spheroid{index}": colors[min((index - 1) // 4, len(colors) - 1)]
    for index in range(1, 21)
}
world.export_usd("chain", "spheroid_chain.usda", parameters, options)
