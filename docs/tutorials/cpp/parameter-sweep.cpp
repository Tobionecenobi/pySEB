#include <cstddef>
#include <iostream>
#include <vector>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("SweepTutorial");
    GraphID graph = world.Add("GaussianPolymer", "blockA");
    world.Link("GaussianPolymer", "blockB.end1", "blockA.end2");
    world.Add(graph, "diblock");

    const DoubleVector q{0.001, 0.01, 0.05, 0.1, 0.3};
    const std::vector<double> rgValues{4.0, 8.0, 16.0};

    std::cout << "Rg_blockB,q,F(q)\n";
    for (double rgB : rgValues) {
        ParameterList parameters;
        world.setParameter(parameters, "beta_blockA", 1.0);
        world.setParameter(parameters, "beta_blockB", 0.8);
        world.setParameter(parameters, "Rg_blockA", 5.0);
        world.setParameter(parameters, "Rg_blockB", rgB);

        DoubleVector values = world.EvaluateFormFactor(
            "diblock", parameters, q
        );
        for (std::size_t index = 0; index < q.size(); ++index) {
            std::cout
                << rgB << ','
                << q[index] << ','
                << values[index] << '\n';
        }
    }
}
