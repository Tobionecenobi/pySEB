#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    if (!sebsym::set_backend("portable")) {
        std::cerr << "portable backend unavailable\n";
        return 1;
    }

    World world("RandomWalkRodTutorial");
    constexpr int segmentCount = 20;

    GraphID chain = world.Add("ThinRod", "rod1", "rod");
    for (int index = 2; index <= segmentCount; ++index) {
        const std::string current = "rod" + std::to_string(index);
        const std::string previous = "rod" + std::to_string(index - 1);
        world.Link(
            "ThinRod",
            current + ".end1",
            previous + ".end2",
            "rod"
        );
    }
    world.Add(chain, "rodChain");

    ParameterList rodParameters;
    world.setParameter(rodParameters, "beta_rod", 1.0);
    world.setParameter(rodParameters, "L_rod", std::sqrt(12.0));

    const double chainRg2 =
        world.EvaluateRadiusOfGyration2("rodChain", rodParameters);

    world.Add("GaussianPolymer", "gaussian");
    ParameterList gaussianParameters;
    world.setParameter(
        gaussianParameters,
        "beta_gaussian",
        static_cast<double>(segmentCount)
    );
    world.setParameter(gaussianParameters, "Rg_gaussian", std::sqrt(chainRg2));

    DoubleVector q{0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0};
    const DoubleVector rodValues =
        world.EvaluateFormFactor("rodChain", rodParameters, q);
    const DoubleVector gaussianValues =
        world.EvaluateFormFactor("gaussian", gaussianParameters, q);

    std::cout << "segments=" << segmentCount << '\n';
    std::cout << "rod_chain_Rg2=" << std::setprecision(10) << chainRg2 << '\n';
    std::cout << "q,rod_chain,gaussian_reference\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << rodValues[index] << ','
                  << gaussianValues[index] << '\n';
    }
}
