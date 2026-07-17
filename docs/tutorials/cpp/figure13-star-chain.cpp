#include <cstddef>
#include <iostream>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("Figure13Tutorial");

    GraphID diblock = world.Add("GaussianPolymer", "A");
    world.Link("GaussianPolymer", "B.end1", "A.end2");

    GraphID star = world.Add(diblock, "diblock1");
    world.Link(diblock, "diblock2:A.end1", "diblock1:A.end1");
    world.Link(diblock, "diblock3:A.end1", "diblock1:A.end1");
    world.Link(diblock, "diblock4:A.end1", "diblock1:A.end1");

    GraphID chain = world.Add(star, "star1");
    world.Link(
        star,
        "star2:diblock1:B.end2",
        "star1:diblock3:B.end2"
    );
    world.Link(
        star,
        "star3:diblock1:B.end2",
        "star2:diblock3:B.end2"
    );
    world.Link(
        star,
        "star4:diblock1:B.end2",
        "star3:diblock3:B.end2"
    );
    world.Link(
        star,
        "star5:diblock1:B.end2",
        "star4:diblock3:B.end2"
    );
    world.Add(chain, "chain");

    ParameterList parameters;
    world.setParameter(parameters, "beta_A", 1.0);
    world.setParameter(parameters, "beta_B", 1.0);
    world.setParameter(parameters, "Rg_A", 1.0);
    world.setParameter(parameters, "Rg_B", 1.0);

    const DoubleVector q{0.01, 0.05, 0.1, 0.5, 1.0};
    const DoubleVector formFactor = world.EvaluateFormFactor(
        "chain", parameters, q
    );

    std::cout << "stars=5,diblocks_per_star=4,subunits=40\n";
    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << formFactor[index] << '\n';
    }
}
