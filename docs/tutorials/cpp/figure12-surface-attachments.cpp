#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("Figure12SurfaceAttachments");
    GraphID graph = world.Add("SolidSphericalShell", "core");
    const int pairCount = 6;

    for (int index = 0; index < pairCount; ++index) {
        const std::string number = std::to_string(index);
        world.Link(
            "GaussianPolymer",
            "polymer" + number + ".end1",
            "core.surfaceo#polymer" + number,
            "polymer"
        );
        world.Link(
            "ThinRod",
            "rod" + number + ".end1",
            "core.surfacei#rod" + number,
            "rod"
        );
    }
    world.Add(graph, "structure");

    ParameterList parameters;
    world.setParameter(parameters, "beta_core", 0.0);
    world.setParameter(parameters, "Ri_core", 8.0);
    world.setParameter(parameters, "Ro_core", 12.0);
    world.setParameter(parameters, "beta_polymer", 1.0);
    world.setParameter(parameters, "Rg_polymer", 1.0);
    world.setParameter(parameters, "beta_rod", 1.0);
    world.setParameter(parameters, "L_rod", std::sqrt(12.0));

    const DoubleVector q{0.05, 0.1, 0.25, 0.5, 1.0};
    const DoubleVector formFactor = world.EvaluateFormFactor(
        "structure", parameters, q
    );

    std::cout << "core=SolidSphericalShell,pairs=" << pairCount << '\n';
    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << formFactor[index] << '\n';
    }
}
