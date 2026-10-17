// Separate native equation generation from rendering, avoiding Python parsing.
#include "World.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#ifdef PYSEB_CALLGRIND
#include <valgrind/callgrind.h>
#endif
int main(int argc, char** argv) {
    const int count = argc > 1 ? std::atoi(argv[1]) : 32;
    if (count < 1) return 1;
    World world("profile_symbolic");
    GraphID graph = world.Add("GaussianPolymer", "p0", "poly");
    for (int i=1; i<count; ++i)
        world.Link("GaussianPolymer", "p"+std::to_string(i)+".end1",
                   "p"+std::to_string(i-1)+".end2", "poly");
    world.Add(graph, "particle");
    ParameterList p{{"Rg_poly",3.}, {"beta_poly",2.}};
    std::vector<double> generation, rendering;
    std::size_t bytes=0;
    double value=0.;
    for (int repeat=0; repeat<4; ++repeat) {
        auto start=std::chrono::steady_clock::now();
#ifdef PYSEB_CALLGRIND
        // With --instr-atstart=no, profile generation without construction,
        // string rendering or evaluation. Ordinary builds need no Valgrind.
        CALLGRIND_START_INSTRUMENTATION;
#endif
        auto expression=world.FormFactor_Unnormalized("particle");
#ifdef PYSEB_CALLGRIND
        CALLGRIND_STOP_INSTRUMENTATION;
#endif
        auto generated=std::chrono::steady_clock::now();
        std::string text=expression.to_python();
        auto rendered=std::chrono::steady_clock::now();
        if (repeat) {
            generation.push_back(std::chrono::duration<double>(generated-start).count());
            rendering.push_back(std::chrono::duration<double>(rendered-generated).count());
        }
        bytes=text.size();
        value=world.Evaluate(expression,p,.1);
    }
    std::sort(generation.begin(),generation.end());
    std::sort(rendering.begin(),rendering.end());
    std::cout<<std::setprecision(17)<<"{\"segments\":"<<count
             <<",\"generation_seconds\":"<<generation[1]
             <<",\"render_seconds\":"<<rendering[1]
             <<",\"python_bytes\":"<<bytes<<",\"value_at_q_0_1\":"<<value<<"}\n";
}
