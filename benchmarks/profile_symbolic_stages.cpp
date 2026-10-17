// Embed Python to time the production SymPy converter separately from native
// generation, without adding profiling hooks to pySEB's public API.
#include <pybind11/embed.h>
#include "bindingsSymbolic.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#ifdef PYSEB_CALLGRIND
#include <valgrind/callgrind.h>
#endif
int main(int argc, char** argv) {
    const std::string topology = argc > 1 ? argv[1] : "chain";
    const int count = argc > 2 ? std::atoi(argv[2]) : 16;
    const int repeats = argc > 4 ? std::atoi(argv[4]) : 7;
    if (repeats < 1 || count < 1 || (topology != "chain" && topology != "star" && topology != "nested")) return 1;
    py::scoped_interpreter interpreter{};
    // Supply the environment's site-packages as argv[3] when embedding outside
    // an activated environment. The repository root provides the Python package.
    if (argc > 3) py::module::import("sys").attr("path").attr("insert")(0, argv[3]);
    py::module::import("sys").attr("path").attr("insert")(0, ".");
    py::module::import("pyseb.symbolic");
    auto sp = py::module::import("sympy");
    const bool profileConversion=argc>5 && std::string(argv[5]) == "profile";
    py::dict globals;
    globals["sp"] = sp;
    // A Python caller frame is required by SymPy's lambdify implementation.
    py::object prepare = py::eval("lambda e: sp.lambdify(sorted(e.free_symbols, key=str), e, modules=['numpy','scipy'], cse=True, docstring_limit=0)", globals);
    World world("independent_symbolic_profile");
    GraphID graph;
    ParameterList p;
    std::string first, last;
    if (topology == "nested") {
        for (int i=0; i<count; ++i) {
            const auto a="a"+std::to_string(i), b="b"+std::to_string(i), block="block"+std::to_string(i);
            const auto inner=world.Add("GaussianPolymer", a);
            world.Link("GaussianPolymer", b+".end1", a+".end2");
            if (!i) graph=world.Add(inner,block);
            else world.Link(inner, block+":"+a+".end1", "block"+std::to_string(i-1)+":b"+std::to_string(i-1)+".end2");
            p["Rg_"+a]=3.+.1*i; p["Rg_"+b]=4.+.2*i;
            p["beta_"+a]=2.; p["beta_"+b]=-1.;
        }
        first="particle:block0:a0.end1";
        last="particle:block"+std::to_string(count-1)+":b"+std::to_string(count-1)+".end2";
    } else {
        if (topology == "star") {
            graph=world.Add("SolidSphere","core");
            p["R_core"]=20.; p["beta_core"]=100.;
        }
        for (int i=0;i<count;++i) {
            const auto name="p"+std::to_string(i);
            if (topology == "chain" && !i) graph=world.Add("GaussianPolymer",name);
            else world.Link("GaussianPolymer",name+".end1",topology == "star"
                ? "core.surface#arm"+std::to_string(i) : "p"+std::to_string(i-1)+".end2");
            p["Rg_"+name]=3.+.1*i; p["beta_"+name]=i%2 ? -1. : 2.;
        }
        first=topology == "star" ? "particle:core.center" : "particle:p0.end1";
        last="particle:p"+std::to_string(count-1)+".end2";
    }
    world.Add(graph,"particle");
    for (const std::string operation : {"intensity","amplitude","phase"}) {
        std::vector<double> generation, conversion;
        double check=0., firstGeneration=0., firstConversion=0.;
        py::object converted;
        py::object profiler=py::none();
        if (profileConversion) profiler=py::module::import("cProfile").attr("Profile")();
        for (int repeat=0;repeat<=repeats;++repeat) {
            const auto start=std::chrono::steady_clock::now();
#ifdef PYSEB_CALLGRIND
            CALLGRIND_START_INSTRUMENTATION;
#endif
            Expression expression=operation == "intensity" ? world.FormFactor_Unnormalized("particle")
                : operation == "amplitude" ? world.FormFactorAmplitude(first) : world.PhaseFactor(first,last);
#ifdef PYSEB_CALLGRIND
            CALLGRIND_STOP_INSTRUMENTATION;
#endif
            const auto generated=std::chrono::steady_clock::now();
            if (profileConversion && !repeat) profiler.attr("enable")();
            converted=pyseb_expression_to_sympy(expression).attr("expr");
            const auto finished=std::chrono::steady_clock::now();
            if (profileConversion && !repeat) {
                profiler.attr("disable")();
                profiler.attr("dump_stats")("/tmp/pyseb-conversion-"+topology+"-"+
                    std::to_string(count)+"-"+operation+".pstats");
            }
            if (!repeat) {
                firstGeneration=std::chrono::duration<double>(generated-start).count();
                firstConversion=std::chrono::duration<double>(finished-generated).count();
            }
            if (repeat) {
                generation.push_back(std::chrono::duration<double>(generated-start).count());
                conversion.push_back(std::chrono::duration<double>(finished-generated).count());
            }
        }
        // Validate the final equation outside both timed stages. Preparing one
        // all-parameter function avoids expensive repeated SymPy substitution.
        const auto prepareStart=std::chrono::steady_clock::now();
        py::object fn=prepare(converted);
        const auto prepareEnd=std::chrono::steady_clock::now();
        // Include ordinary downstream uses: a conversion improvement must not
        // merely postpone expensive work until lambdify or the first printer.
        const auto textStart=std::chrono::steady_clock::now();
        const auto text=py::cast<std::string>(py::str(converted));
        const auto textEnd=std::chrono::steady_clock::now();
        const auto latex=py::cast<std::string>(sp.attr("latex")(converted));
        const auto latexEnd=std::chrono::steady_clock::now();
        py::list symbols=py::module::import("builtins").attr("sorted")(
            converted.attr("free_symbols"),py::arg("key")=py::module::import("builtins").attr("str"));
        for (double q : {.03,.1,.3}) {
            py::tuple arguments(py::len(symbols));
            for (std::size_t i=0;i<py::len(symbols);++i) {
                const auto name=py::cast<std::string>(py::str(symbols[i]));
                arguments[i]=name == "q" ? q : p.at(name);
            }
            const double expected=operation == "intensity" ? world.EvaluateFormFactorUnnormalized("particle",p,q)
                : operation == "amplitude" ? world.EvaluateFormFactorAmplitude(first,p,q) : world.EvaluatePhaseFactor(first,last,p,q);
            check=py::cast<double>(fn(*arguments));
            if (std::abs(check-expected)>1e-8+2e-8*std::abs(expected)) throw std::runtime_error("Equation parity failed");
        }
        py::tuple curveArguments(py::len(symbols));
        const auto qValues=py::module::import("numpy").attr("linspace")(.03,.3,31);
        for (std::size_t i=0;i<py::len(symbols);++i) {
            const auto name=py::cast<std::string>(py::str(symbols[i]));
            curveArguments[i]=name == "q" ? qValues : py::cast(p.at(name));
        }
        fn(*curveArguments); // Warm numerical dispatch before measuring reuse.
        const auto evaluateStart=std::chrono::steady_clock::now();
        for (int i=0;i<100;++i) fn(*curveArguments);
        const auto evaluateEnd=std::chrono::steady_clock::now();
        double simplifySeconds=-1.;
        if (argc>5 && std::string(argv[5]) == "simplify") {
            // Opt in on small cases: general symbolic simplification has no
            // practical runtime bound for a large independently tagged graph.
            const auto start=std::chrono::steady_clock::now();
            sp.attr("simplify")(converted);
            simplifySeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        }
        std::sort(generation.begin(),generation.end()); std::sort(conversion.begin(),conversion.end());
        std::cout<<std::setprecision(17)<<"{\"topology\":\""<<topology<<"\",\"size\":"<<count
                 <<",\"sympy_version\":\""<<py::cast<std::string>(sp.attr("__version__"))<<"\""
                 <<",\"sympy_cache_size\":\""<<py::cast<std::string>(py::str(
                     py::module::import("sympy.core.cache").attr("SYMPY_CACHE_SIZE")))<<"\""
                 <<",\"first_conversion_profiled\":"<<(profileConversion ? "true" : "false")
                 <<",\"operation\":\""<<operation<<"\",\"generation_seconds\":"<<generation[generation.size()/2]
                 <<",\"conversion_seconds\":"<<conversion[conversion.size()/2]
                 <<",\"first_generation_seconds\":"<<firstGeneration
                 <<",\"first_conversion_seconds\":"<<firstConversion
                 <<",\"prepare_seconds\":"<<std::chrono::duration<double>(prepareEnd-prepareStart).count()
                 <<",\"text_seconds\":"<<std::chrono::duration<double>(textEnd-textStart).count()
                 <<",\"latex_seconds\":"<<std::chrono::duration<double>(latexEnd-textEnd).count()
                 <<",\"text_bytes\":"<<text.size()<<",\"latex_bytes\":"<<latex.size()
                 <<",\"evaluate_seconds\":"<<std::chrono::duration<double>(evaluateEnd-evaluateStart).count()/100.
                 <<",\"simplify_seconds\":"<<simplifySeconds
                 <<",\"value_at_q_0_3\":"<<check<<"}"<<std::endl;
    }
}
