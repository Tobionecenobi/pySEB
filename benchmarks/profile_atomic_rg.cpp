#include "Subunits/AtomicSAXSCloud.hpp"
#include <iostream>
#include <vector>
// A small native workload for gprof: keep construction out of the repeated
// Rg calculation and use a fixed cloud so before/after results are comparable.
int main() {
    std::vector<AtomicSAXSScatterer> atoms;
    for (int i=0;i<512;++i) {
        AtomicSAXSScatterer a;
        a.center=CartesianPoint3D(i%8,(i/8)%8,i/64);
        a.vacuum=CromerMannCoefficients({{2.,1.,1.,1.}},{{2.,4.,8.,16.}},1.);
        a.excludedVolume=12.; a.solventAccessibility=.5;
        atoms.push_back(a);
    }
    AtomicSAXSOptions opts; opts.includeHydration=true; opts.hydrationShellContrast=.02;
    AtomicSAXSCloud cloud(atoms,opts);
    ParameterList p;
    double sum=0.;
    for(int i=0;i<100;++i) sum+=cloud.NumericRadiusOfGyration2(p);
    std::cout<<sum<<'\n';
}
