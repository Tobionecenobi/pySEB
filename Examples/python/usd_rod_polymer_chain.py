"""Export an alternating chain of rods and Gaussian polymers joined end-to-end."""
import pyseb


world = pyseb.World("rod_polymer_chain")
members = []
for index in range(1, 11):
    members.append(("ThinRod", f"rod{index}", "rod"))
    members.append(("GaussianPolymer", f"polymer{index}", "polymer"))

first_model, first_name, first_tag = members[0]
graph = world.Add(first_model, first_name, first_tag)
for (model, name, tag), (_, previous, _) in zip(members[1:], members):
    graph = world.Link(model, f"{name}.end1", f"{previous}.end2", tag)

world.Add(graph, "chain")

parameters = {"L_rod": 2.0, "Rg_polymer": 0.55}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 32
options.relaxation_sweeps = 6
options.minimum_clearance = 0.1
options.curve_samples = 96
options.color_overrides = {
    name: (0.91, 0.42, 0.12) if model == "ThinRod" else (0.10, 0.47, 0.78)
    for model, name, _ in members
}
world.export_usd("chain", "rod_polymer_chain.usda", parameters, options)
