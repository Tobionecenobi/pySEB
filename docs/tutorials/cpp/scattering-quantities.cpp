#include <cstddef>
#include <iostream>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    DebyeSphereCloud* cloud = new DebyeSphereCloud({
        SphereScatterer(0.0, 0.0, 0.0, 0.0, 1.0),
        SphereScatterer(2.0, 0.0, 0.0, 0.0, 1.0)
    });
    cloud->addReferencePoint("left", 0.0, 0.0, 0.0);
    cloud->addReferencePoint("right", 2.0, 0.0, 0.0);

    World world("ScatteringQuantitiesTutorial");
    world.Add(cloud, "cloud");

    const ParameterList parameters;
    const DoubleVector q{0.0, 0.1, 0.3, 0.6, 1.0};
    const DoubleVector formFactor = world.EvaluateFormFactor(
        "cloud", parameters, q
    );
    const DoubleVector unnormalized = world.EvaluateFormFactorUnnormalized(
        "cloud", parameters, q
    );
    const DoubleVector amplitude = world.EvaluateFormFactorAmplitude(
        "cloud.left", parameters, q
    );
    const DoubleVector phase = world.EvaluatePhaseFactor(
        "cloud.left", "cloud.right", parameters, q
    );

    std::cout << "q,F(q),I(q),A_left(q),Psi_left_right(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ','
                  << formFactor[index] << ','
                  << unnormalized[index] << ','
                  << amplitude[index] << ','
                  << phase[index] << '\n';
    }

    std::cout << "Rg2="
              << world.EvaluateRadiusOfGyration2("cloud", parameters)
              << '\n';
    std::cout << "left_to_right_msd="
              << world.EvaluateSMSDRef2Ref(
                     "cloud.left", "cloud.right", parameters
                 )
              << '\n';
}
