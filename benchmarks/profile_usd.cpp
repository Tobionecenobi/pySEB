// Native USD workloads. Model construction is excluded; exported scenes stay
// on disk for byte comparisons and optional inspection with OpenUSD's loader.
#include "World.hpp"
#include "USDExport.hpp"
#include "Subunits/DebyeSphereCloud.hpp"
#include "Subunits/AtomicSAXSCloud.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#ifdef __linux__
#include <sys/resource.h>
#endif

namespace {
using Clock = std::chrono::steady_clock;
std::map<std::string,double> stages;
std::string current;
Clock::time_point previous;
double envelopeSeconds = 0.;
}
// Called only by the benchmark's instrumented copy of USDExport.cpp.
void usdProfileStage(const char* stage) {
    const auto now = Clock::now();
    if (!current.empty()) stages[current] += std::chrono::duration<double>(now-previous).count();
    current = stage;
    previous = now;
}
void usdProfileEnvelope(double seconds) { envelopeSeconds += seconds; }

int main(int argc, char** argv) {
    if (argc != 7) {
        std::cerr << "usage: profile-usd analytical|debye|atomic count random|readable envelope|points repeats output.usda\n";
        return 1;
    }
    const std::string kind=argv[1], layout=argv[3], envelope=argv[4], path=argv[6];
    const int count=std::atoi(argv[2]), repeats=std::atoi(argv[5]);
    if (count<1 || repeats<1 || (kind!="analytical" && kind!="debye" && kind!="atomic")
        || (layout!="random" && layout!="readable") || (envelope!="envelope" && envelope!="points")) return 1;
    World world("usd_benchmark");
    GraphID graph;
    ParameterList parameters;
    if (kind=="analytical") {
        // A large assembly with distributed, end-to-end cylinder joints.
        graph=world.Add("SolidCylinder","c0","cylinder");
        for (int i=1; i<count; ++i)
            world.Link("SolidCylinder","c"+std::to_string(i)+".ends#bottom#joint"+std::to_string(i),
                "c"+std::to_string(i-1)+".ends#top#joint"+std::to_string(i),"cylinder");
        parameters={{"R_cylinder",.25},{"L_cylinder",3.}};
    } else {
        // Deterministic, irregularly spaced clouds exercise signed weights,
        // mixed radii, point instancing and cloud-to-analytical attachments.
        if (kind=="debye") {
            std::vector<SphereScatterer> points;
            for (int i=0; i<count; ++i) points.emplace_back(
                .7*(i%13), .9*((i/13)%13), .8*(i/169), i%7 ? .3 : 0., i%2 ? -1. : 2.);
            auto* cloud=new DebyeSphereCloud(points);
            cloud->addReferencePoint("join",0.,0.,0.);
            graph=world.Add(cloud,"cloud");
        } else {
            std::vector<AtomicSAXSScatterer> points;
            for (int i=0; i<count; ++i) {
                AtomicSAXSScatterer point;
                point.center=CartesianPoint3D(.7*(i%13),.9*((i/13)%13),.8*(i/169));
                point.excludedVolume=12.;
                point.vacuum=CromerMannCoefficients({{2.,1.,1.,1.}},{{2.,4.,8.,16.}},1.);
                points.push_back(point);
            }
            auto* cloud=new AtomicSAXSCloud(points);
            cloud->addReferencePoint("join",CartesianPoint3D());
            graph=world.Add(cloud,"cloud");
        }
        world.Link("SolidSphere","sphere.surface#anchor","cloud.join");
        parameters={{"R_sphere",2.},{"beta_sphere",-10.}};
    }
    world.Add(graph,"particle");
    USDExportOptions options(LengthUnit::Angstrom);
    options.seed=42;
    options.layoutMode=layout=="readable" ? USDLayoutMode::Readable : USDLayoutMode::Random;
    options.debyeEnvelope=envelope=="envelope";
    // Keep the normal quality settings, including 64 trials and four sweeps.
    auto run=[&]() { world.ExportUSD("particle",path,parameters,options); };
    run(); // Exclude lazy initialization and filesystem warm-up.
    for (int repeat=0; repeat<repeats; ++repeat) {
        stages.clear(); current.clear(); envelopeSeconds=0.;
        const auto start=Clock::now();
        run();
        const double total=std::chrono::duration<double>(Clock::now()-start).count();
        std::ifstream file(path,std::ios::binary|std::ios::ate);
        if (!file) throw std::runtime_error("missing exported scene");
        const double cleanup=total-stages["prepare"]-stages["layout"]-stages["serialize"]-stages["write"];
        std::cout << std::setprecision(17) << "{\"kind\":\""<<kind<<"\",\"count\":"<<count
            <<",\"layout\":\""<<layout<<"\",\"envelope\":"<<(options.debyeEnvelope?"true":"false")
            <<",\"repeat\":"<<repeat<<",\"bytes\":"<<file.tellg()<<",\"total_seconds\":"<<total
            <<",\"prepare_seconds\":"<<stages["prepare"]<<",\"layout_seconds\":"<<stages["layout"]
            <<",\"geometry_serialization_seconds\":"<<stages["serialize"]-envelopeSeconds
            <<",\"envelope_seconds\":"<<envelopeSeconds<<",\"write_seconds\":"<<stages["write"]
            <<",\"cleanup_seconds\":"<<cleanup;
#ifdef __linux__
        // Linux peak process RSS includes construction, warm-up and all prior
        // repeats; compare separate fresh processes with identical workloads.
        rusage usage{};
        getrusage(RUSAGE_SELF,&usage);
        std::cout <<",\"peak_rss_kib\":"<<usage.ru_maxrss;
#endif
        std::cout <<"}\n";
    }
}
