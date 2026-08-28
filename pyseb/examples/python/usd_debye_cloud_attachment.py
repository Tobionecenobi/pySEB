"""Export an intrinsic Debye sphere cloud attached to an analytic spheroid."""
import math

import pyseb


scatterers = []
golden_angle = math.pi * (3.0 - math.sqrt(5.0))
for index in range(28):
    fraction = index / 27.0
    angle = index * golden_angle
    radial = 0.35 + 0.3 * math.sin(math.pi * fraction)
    x = -2.8 + 2.7 * fraction
    y = radial * math.cos(angle)
    z = radial * math.sin(angle)
    radius = 0.0 if index == 7 else 0.13 + 0.025 * (index % 4)
    beta = -0.4 if index % 6 == 0 else 1.0
    scatterers.append(pyseb.SphereScatterer(x, y, z, radius, beta))
scatterers.append(pyseb.SphereScatterer(0.0, 0.0, 0.0, 0.2, 1.0))

cloud = pyseb.DebyeSphereCloud(scatterers)
cloud.addReferencePoint("join", 0.2, 0.0, 0.0)

world = pyseb.World("debye_cloud_attachment")
graph = world.Add(cloud, "cloud")
world.Link("Spheroid", "spheroid.south", "cloud.join", "spheroid")
world.Add(graph, "assembly")

parameters = {"a_spheroid": 0.65, "c_spheroid": 1.25}
options = pyseb.USDExportOptions(pyseb.LengthUnit.Angstrom)
options.seed = 42
options.layout_mode = pyseb.LayoutMode.Readable
options.orientation_trials = 64
options.relaxation_sweeps = 4
options.minimum_clearance = 0.15
options.surface_samples = 32
options.zero_radius_marker_size = 0.08
options.color_overrides = {
    "cloud": (0.10, 0.55, 0.78),
    "spheroid": (0.92, 0.43, 0.12),
}
world.export_usd("assembly", "debye_cloud_attachment.usda", parameters, options)
