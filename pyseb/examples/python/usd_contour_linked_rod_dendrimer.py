"""Export a rod dendrimer whose links join two distributed contour points."""
import pyseb


world = pyseb.World("contour_linked_rod_dendrimer")
graph = world.Add("ThinRod", "rod1", "rod")

parents = ["rod1"]
rod_generation = {"rod1": 0}
next_index = 2
generation_count = 4
branches_per_parent = 2

for generation in range(1, generation_count + 1):
    children = []
    for parent in parents:
        for branch in range(1, branches_per_parent + 1):
            child = f"rod{next_index}"
            graph = world.Link(
                "ThinRod",
                f"{child}.contour#attachment",
                f"{parent}.contour#g{generation}_branch{branch}",
                "rod",
            )
            children.append(child)
            rod_generation[child] = generation
            next_index += 1
    parents = children

world.Add(graph, "dendrimer")

parameters = {"L_rod": 1.0}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 64
options.minimum_clearance = 0.08
options.curve_samples = 2
generation_colors = (
    (0.19, 0.20, 0.23),
    (0.15, 0.36, 0.61),
    (0.23, 0.56, 0.48),
    (0.88, 0.62, 0.24),
    (0.77, 0.27, 0.21),
)
options.color_overrides = {
    rod: generation_colors[generation]
    for rod, generation in rod_generation.items()
}
world.export_usd(
    "dendrimer",
    "contour_linked_rod_dendrimer.usda",
    parameters,
    options,
)
