#include <cstddef>
#include <iostream>
#include <iomanip>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    DebyeSphereCloud* cloud = new DebyeSphereCloud({
        SphereScatterer(-2.0, -0.65, -0.4, 0.0, 1.0),
        SphereScatterer(-1.2,  1.00,  0.3, 0.0, 2.0),
        SphereScatterer(-0.4, -1.25,  0.6, 0.0, 1.0),
        SphereScatterer( 0.4,  1.25, -0.6, 0.0, 1.0),
        SphereScatterer( 1.2, -1.00, -0.3, 0.0, 2.0),
        SphereScatterer( 2.0,  0.65,  0.4, 0.0, 1.0)
    });
    cloud->addReferencePoint("origin", 0.0, 0.0, 0.0);
    cloud->addReferencePoint("anchor", -3.0, 0.0, 0.0);

    World world("PointCloudFactorsTutorial");
    world.Add(cloud, "cloud");
    const ParameterList parameters; // Coordinates and weights are in the cloud.
    const DoubleVector q{0.0, 0.1, 0.3, 0.6, 1.0}; // nm^-1; coordinates in nm.
    const DoubleVector formFactor = world.EvaluateFormFactor("cloud", parameters, q);
    const DoubleVector amplitude = world.EvaluateFormFactorAmplitude(
        "cloud.anchor", parameters, q
    );
    const DoubleVector phase = world.EvaluatePhaseFactor(
        "cloud.origin", "cloud.anchor", parameters, q
    );

    std::cout << std::setprecision(12);
    std::cout << "q,F(q),A_anchor(q),Psi_origin_anchor(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << formFactor[index] << ','
                  << amplitude[index] << ',' << phase[index] << '\n';
    }
}
