#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <string>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("Figure11Dendrimer");
    GraphID graph = world.Add("Point", "center");
    const int functionality = 3;
    int nextIndex = 0;

    std::function<void(int, const std::string&)> attach;
    attach = [&](int generation, const std::string& parentReference) {
        const int childCount = parentReference == "center.point"
            ? functionality
            : functionality - 1;

        for (int child = 0; child < childCount; ++child) {
            const int index = nextIndex++;
            const std::string name = "branch" + std::to_string(index);
            world.Link(
                "ThinRod",
                name + ".end1",
                parentReference,
                "branch"
            );
            if (generation > 1) {
                attach(generation - 1, name + ".end2");
            }
        }
    };

    attach(3, "center.point");
    world.Add(graph, "dendrimer");

    ParameterList parameters;
    world.setParameter(parameters, "beta_branch", 1.0);
    world.setParameter(parameters, "L_branch", std::sqrt(12.0));

    const DoubleVector q{0.01, 0.05, 0.1, 0.5, 1.0};
    const DoubleVector formFactor = world.EvaluateFormFactor(
        "dendrimer", parameters, q
    );

    std::cout << "generations=3,functionality=3,branches=" << nextIndex << '\n';
    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << formFactor[index] << '\n';
    }
}
