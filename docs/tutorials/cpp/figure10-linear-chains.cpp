#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

#include "SEB.hpp"
#include "Symbolic.hpp"

DoubleVector evaluateRodChain(bool contourLinks, const DoubleVector& q)
{
    World world(contourLinks ? "Figure10Contour" : "Figure10Endpoints");
    const int segmentCount = 12;
    GraphID graph = world.Add("ThinRod", "segment1", "segment");

    for (int index = 2; index <= segmentCount; ++index) {
        const std::string currentName = "segment" + std::to_string(index);
        const std::string previousName = "segment" + std::to_string(index - 1);
        std::string currentReference = currentName + ".end1";
        std::string previousReference = previousName + ".end2";

        if (contourLinks) {
            currentReference = currentName + ".contour#left" +
                std::to_string(index);
            previousReference = previousName + ".contour#right" +
                std::to_string(index - 1);
        }

        world.Link(
            "ThinRod",
            currentReference,
            previousReference,
            "segment"
        );
    }
    world.Add(graph, "chain");

    ParameterList parameters;
    world.setParameter(parameters, "beta_segment", 1.0);
    world.setParameter(parameters, "L_segment", std::sqrt(12.0));
    return world.EvaluateFormFactor("chain", parameters, q);
}

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    const DoubleVector q{0.01, 0.05, 0.1, 0.5, 1.0};
    const DoubleVector endpoints = evaluateRodChain(false, q);
    const DoubleVector contours = evaluateRodChain(true, q);

    std::cout << "segments=12,model=ThinRod\n";
    std::cout << "q,end_to_end,contour_to_contour\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << endpoints[index] << ','
                  << contours[index] << '\n';
    }
}
