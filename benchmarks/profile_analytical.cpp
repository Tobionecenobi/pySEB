// Native profiling driver: repeated curve evaluation of a connected polymer chain.
// Instrument World.cpp and Subunit.cpp to resolve traversal and expression costs.
#include "World.hpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
int main(int argc, char** argv) {
    const int count = argc > 1 ? std::atoi(argv[1]) : 32;
    if (count < 1) return 1;
    World world("profile_analytical");
    GraphID graph = world.Add("GaussianPolymer", "p0", "poly");
    for (int i = 1; i < count; ++i)
        world.Link("GaussianPolymer", "p"+std::to_string(i)+".end1",
                   "p"+std::to_string(i-1)+".end2", "poly");
    world.Add(graph, "particle");
    ParameterList parameters{{"Rg_poly", 8.}, {"beta_poly", 10.}};
    DoubleVector q;
    for (int i = 0; i < 31; ++i) q.push_back(.001 + .499*i/30.);
    world.EvaluateFormFactorUnnormalized("particle", parameters, q);
    const auto start = std::chrono::steady_clock::now();
    double checksum = 0.;
    for (int repeat = 0; repeat < 3; ++repeat)
        for (double intensity : world.EvaluateFormFactorUnnormalized("particle", parameters, q))
            checksum += intensity;
    std::cout << "Elapsed seconds: "
              << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
              << "; checksum: " << checksum << '\n';
}
