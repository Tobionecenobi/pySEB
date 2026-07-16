"""Numerical Debye scattering from a small cloud of finite spheres."""

import pyseb


scatterers = [
    pyseb.SphereScatterer(0.0, 0.0, 0.0, 1.5, 1.0),
    pyseb.SphereScatterer(4.0, 0.0, 0.0, 1.0, 0.7),
]

cloud = pyseb.DebyeSphereCloud(scatterers)
cloud.addReferencePoint("left", 0.0, 0.0, 0.0)
cloud.addReferencePoint("right", 4.0, 0.0, 0.0)

world = pyseb.World()
world.Add(cloud, "cloud")

q = [0.01, 0.03, 0.1, 0.3]
form_factor = world.EvaluateFormFactor("cloud", {}, q)
left_amplitude = world.EvaluateFormFactorAmplitude("cloud.left", {}, q)
end_to_end_phase = world.EvaluatePhaseFactor(
    "cloud.left", "cloud.right", {}, q
)

print("q:", q)
print("F(q):", form_factor)
print("A_left(q):", left_amplitude)
print("Psi_left,right(q):", end_to_end_phase)
