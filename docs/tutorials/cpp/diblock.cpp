#include <cstddef>
#include <iostream>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("DiblockTutorial");
    GraphID graph = world.Add("GaussianPolymer", "blockA");
    world.Link("GaussianPolymer", "blockB.end1", "blockA.end2");
    world.Add(graph, "diblock");

    ParameterList parameters;
    world.setParameter(parameters, "beta_blockA", 1.0);
    world.setParameter(parameters, "beta_blockB", 0.8);
    world.setParameter(parameters, "Rg_blockA", 5.0);
    world.setParameter(parameters, "Rg_blockB", 8.0);

    DoubleVector q{0.001, 0.01, 0.05, 0.1, 0.3};
    DoubleVector values = world.EvaluateFormFactor(
        "diblock", parameters, q
    );

    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << values[index] << '\n';
    }
}
