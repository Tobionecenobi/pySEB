// Compare complete native vector curves with the unbatched scalar reference.
// Construction, SASA and warm-up are excluded. No network data are needed.
#include "World.hpp"
#include "StructureIO/AtomicSAXSCloudBuilder.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#ifdef PYSEB_CALLGRIND
#include <valgrind/callgrind.h>
#endif

int main(int argc, char** argv) {
    if (argc < 3) return 1;
    const std::string mode = argv[1];
    const std::string source = argv[2];
    const int repeats = argc > 3 ? std::atoi(argv[3]) : 3;
    const bool hybrid = argc > 4 && std::string(argv[4]) == "hybrid";
    if ((mode != "scalar" && mode != "vector") || repeats < 1) return 1;
    AtomicSAXSOptions options;
    options.includeHydration = true;
    options.hydrationShellContrast = .02;
    std::unique_ptr<AtomicSAXSCloud> cloud;
    if (source == "protein") {
        // Select only model 1: different NMR conformations never coexist.
        StructureParseOptions parse;
        parse.modelNumber = 1;
        cloud = AtomicSAXSStructureLoader::LoadPDB(
            "docs/tutorials/data/1TIU-heavy-ensemble.pdb",
            AtomicSAXSParameterProfile::ProteinDefaults(), parse, options);
    } else {
        const int count = std::atoi(source.c_str());
        if (count < 1) return 1;
        std::vector<AtomicSAXSScatterer> atoms;
        for (int i = 0; i < count; ++i) {
            AtomicSAXSScatterer atom;
            atom.center = CartesianPoint3D(i%8, (i/8)%8, i/64);
            atom.vacuum = CromerMannCoefficients({{2.,1.,1.,1.}}, {{2.,4.,8.,16.}}, 1.);
            atom.excludedVolume = 12.;
            atom.solventAccessibility = .5;
            atoms.push_back(atom);
        }
        cloud.reset(new AtomicSAXSCloud(atoms, options));
    }
    const auto count = cloud->scattererCount();
    cloud->addReferencePoint("join", CartesianPoint3D(0.,0.,0.));
    World world("cloud_batch_benchmark");
    const auto graph = world.Add(cloud.release(), "cloud");
    if (hybrid) world.Link("SolidSphere", "sphere.center", "cloud.join");
    world.Add(graph, "particle");
    ParameterList parameters{{"R_sphere", 4.}, {"beta_sphere", -20.}};
    DoubleVector q;
    for (int i = 0; i < 31; ++i) q.push_back(.001+.499*i/30.);
    auto evaluate = [&]() {
        if (mode == "vector") return world.EvaluateFormFactorUnnormalized("particle", parameters, q);
        DoubleVector values;
        for (double value : q)
            values.push_back(world.EvaluateFormFactorUnnormalized("particle", parameters, value));
        return values;
    };
    evaluate();
    const auto start = std::chrono::steady_clock::now();
    DoubleVector curve;
#ifdef PYSEB_CALLGRIND
    CALLGRIND_START_INSTRUMENTATION;
#endif
    for (int repeat = 0; repeat < repeats; ++repeat) curve = evaluate();
#ifdef PYSEB_CALLGRIND
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout << std::setprecision(17) << "{\"source\":\""
              << (source == "protein" ? "protein" : "synthetic") << "\",\"mode\":\""
              << mode << "\",\"hybrid\":" << (hybrid ? "true" : "false")
              << ",\"atoms\":" << count << ",\"seconds_per_curve\":"
              << seconds/repeats << ",\"curve\":[";
    for (std::size_t i = 0; i < curve.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << curve[i];
    }
    std::cout << "]}\n";
}
