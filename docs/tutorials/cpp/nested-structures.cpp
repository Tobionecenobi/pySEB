#include <cstddef>
#include <iostream>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("NestedTutorial");

    GraphID diblock = world.Add("GaussianPolymer", "polyA");
    world.Link("GaussianPolymer", "polyB.end1", "polyA.end2");

    GraphID star = world.Add(diblock, "diblock1");
    world.Link(
        diblock,
        "diblock2:polyA.end1",
        "diblock1:polyA.end1"
    );

    GraphID chain = world.Add(star, "star1");
    world.Link(
        star,
        "star2:diblock1:polyB.end2",
        "star1:diblock2:polyB.end2"
    );
    world.Add(chain, "starChain");

    ParameterList parameters;
    world.setParameter(parameters, "beta_polyA", 1.0);
    world.setParameter(parameters, "beta_polyB", 0.75);
    world.setParameter(parameters, "Rg_polyA", 5.0);
    world.setParameter(parameters, "Rg_polyB", 8.0);

    DoubleVector q{0.01, 0.03, 0.1, 0.3};
    DoubleVector values = world.EvaluateFormFactor(
        "starChain", parameters, q
    );

    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << values[index] << '\n';
    }
}
