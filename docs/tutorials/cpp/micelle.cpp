#include <cstddef>
#include <iostream>
#include <string>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("MicelleTutorial");
    GraphID graph = world.Add("SolidSphere", "core");

    for (int index = 0; index < 6; ++index) {
        const std::string number = std::to_string(index);
        world.Link(
            "GaussianPolymer",
            "chain" + number + ".end1",
            "core.surface#anchor" + number,
            "chain"
        );
    }
    world.Add(graph, "micelle");

    ParameterList parameters;
    world.setParameter(parameters, "beta_core", 1.0);
    world.setParameter(parameters, "R_core", 30.0);
    world.setParameter(parameters, "beta_chain", 0.25);
    world.setParameter(parameters, "Rg_chain", 12.0);

    DoubleVector q{0.005, 0.02, 0.05, 0.1, 0.2};
    DoubleVector values = world.EvaluateFormFactor(
        "micelle", parameters, q
    );

    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << values[index] << '\n';
    }
}
