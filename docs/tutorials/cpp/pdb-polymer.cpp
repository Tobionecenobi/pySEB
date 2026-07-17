#include <cstddef>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    const std::string pdbText =
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  "
        "1.00 20.00           N\n"
        "ATOM      2  C   GLY A   1       4.000   0.000   0.000  "
        "1.00 20.00           C\n"
        "END\n";

    std::istringstream input(pdbText);
    PDBParser parser;
    std::vector<AtomRecord> atoms = parser.Parse(input);

    AtomParameterProfile profile;
    profile.SetElement("N", 1.55, 7.0);
    profile.SetElement("C", 1.70, 6.0);

    AtomCloudBuildOptions build;
    build.referenceAtomSerials["attachment"] = 2;
    std::unique_ptr<DebyeSphereCloud> cloud =
        AtomCloudBuilder::Build(atoms, profile, build);

    World world("PDBTutorial");
    GraphID graph = world.Add(cloud.release(), "fragment");
    world.Link(
        "GaussianPolymer",
        "polymer.end1",
        "fragment.attachment"
    );
    world.Add(graph, "conjugate");

    ParameterList parameters;
    world.setParameter(parameters, "beta_polymer", 4.0);
    world.setParameter(parameters, "Rg_polymer", 6.0);

    DoubleVector q{0.0, 0.02, 0.05, 0.1, 0.2};
    DoubleVector values = world.EvaluateFormFactor(
        "conjugate", parameters, q
    );

    std::cout << "q,F(q)\n";
    for (std::size_t index = 0; index < q.size(); ++index) {
        std::cout << q[index] << ',' << values[index] << '\n';
    }
}
