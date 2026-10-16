// Profile complete vector intensities, keeping construction and geometry warm-up
// outside Callgrind's measured region. Ordinary builds need no Valgrind headers.
#include "World.hpp"
#include "Subunits/AtomicSAXSCloud.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#ifdef PYSEB_CALLGRIND
#include <valgrind/callgrind.h>
#endif

int main(int argc, char** argv) {
    const std::string kind = argc > 1 ? argv[1] : "chain";
    const int count = argc > 2 ? std::atoi(argv[2]) : 16;
    const int repeats = argc > 3 ? std::atoi(argv[3]) : 3;
    if (count < 1 || repeats < 1 ||
        (kind != "chain" && kind != "cloud" && kind != "hybrid")) return 1;
    World world("profile_numerical");
    ParameterList parameters{{"Rg_poly", 8.}, {"beta_poly", 10.},
                             {"R_sphere", 4.}, {"beta_sphere", 20.}};
    GraphID graph;
    if (kind == "chain") {
        graph = world.Add("GaussianPolymer", "p0", "poly");
        for (int i = 1; i < count; ++i)
            world.Link("GaussianPolymer", "p"+std::to_string(i)+".end1",
                       "p"+std::to_string(i-1)+".end2", "poly");
    } else {
        // Fixed coordinates and nonzero solvent/hydration make comparisons
        // reproducible while retaining all cloud pair terms and hybrid terms.
        std::vector<AtomicSAXSScatterer> atoms;
        for (int i = 0; i < count; ++i) {
            AtomicSAXSScatterer atom;
            atom.center = CartesianPoint3D(i%8, (i/8)%8, i/64);
            atom.vacuum = CromerMannCoefficients({{2.,1.,1.,1.}}, {{2.,4.,8.,16.}}, 1.);
            atom.excludedVolume = 12.;
            atom.solventAccessibility = .5;
            atoms.push_back(atom);
        }
        AtomicSAXSOptions options;
        options.includeHydration = true;
        options.hydrationShellContrast = .02;
        auto* cloud = new AtomicSAXSCloud(atoms, options);
        cloud->addReferencePoint("join", CartesianPoint3D(0.,0.,0.));
        graph = world.Add(cloud, "cloud");
        if (kind == "hybrid")
            world.Link("SolidSphere", "sphere.surface#join", "cloud.join", "sphere");
    }
    world.Add(graph, "particle");
    DoubleVector q;
    for (int i = 0; i < 31; ++i) q.push_back(.001 + .499*i/30.);
    world.EvaluateFormFactorUnnormalized("particle", parameters, q);
    const auto start = std::chrono::steady_clock::now();
    double checksum = 0.;
#ifdef PYSEB_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif
    for (int repeat = 0; repeat < repeats; ++repeat)
        for (double value : world.EvaluateFormFactorUnnormalized("particle", parameters, q))
            checksum += value;
#ifdef PYSEB_CALLGRIND
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
    std::cout << std::setprecision(17) << "{\"case\":\"" << kind
              << "\",\"count\":" << count << ",\"repeats\":" << repeats
              << ",\"seconds\":"
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << ",\"checksum\":" << checksum << "}\n";
}
