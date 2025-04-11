/*

    Samples scattering contributions from a thin spherical disk of radius R


*/ 

#include <iostream>
#include <random>
#include "DebyeSampler.hpp"
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

default_random_engine re;

#include "GenerateSamples.hpp"

void Sample(double R) {
    string dir = "ThinDisk_R"+to_string(R)+"/";
    fs::create_directory(dir);
    
    uniform_real_distribution<double> dist(-R,R);
    uniform_real_distribution<double> unit(0,1);
    
    // Reduce sample size for testing
    #ifdef TESTING
    int N = 1000;
    #else
    int N = 100000000;
    #endif
    
    vector<double> qvec;
    double qmin=0.1;
    double qmax=100;
 
    double q=qmin;
    while (q<qmax) {
        qvec.push_back(q);
        q*=1.03;
    }

    Sampler FF(qvec,                dir+"FF.q");
    Sampler FFAcenter(qvec,         dir+"FFA_center.q");
    Sampler FFArim(qvec,            dir+"FFA_rim.q");
    Sampler PF_contour_contour(qvec, dir+"PF_contour_contour.q");
    Sampler PF_center_contour(qvec,  dir+"PF_center_contour.q");
    Sampler Prim2rim(qvec,          dir+"P_rim2rim.q");
     
    for (int i=0; i<N; i++) {
        if (i % (N/100)==0)
            cout << 100.0*i/N << "% done\n";

        // Two points on the surface of the spherical disk      
        double x1,y1,z1, x2,y2,z2;
        
        RandomDisk(R,x1,y1,z1, dist);
        RandomDisk(R,x2,y2,z2, dist);

        // Two points on the spherical rim
        double rx1,ry1,rz1, rx2,ry2,rz2;
        RandomCircle(R,rx1,ry1,rz1, dist);
        RandomCircle(R,rx2,ry2,rz2, dist);

        // Form factor is the pair distance
        FF.add(x1-x2,y1-y2,z1-z2);

        // Form factor amplitudes         
        // Center to scatterer
        FFAcenter.add(x1,y1,z1);
        FFAcenter.add(x2,y2,z2);

        // Rim to scatterer
        FFArim.add(x1-rx1,y1-ry1,z1-rz1);
        FFArim.add(x2-rx2,y2-ry2,z2-rz2);
        FFArim.add(x1-rx2,y1-ry2,z1-rz2);
        FFArim.add(x2-rx1,y2-ry1,z2-rz1);

        // Phase factors
        PF_contour_contour.add(x2-x1,y2-y1,z2-z1);
        PF_center_contour.add(x1,y1,z1);
        PF_center_contour.add(x2,y2,z2);

        // Rim to rim
        Prim2rim.add(rx2-rx1,ry2-ry1,rz2-rz1);
    }
}

// Keep main for standalone usage
int main() {
    Sample(1.0);
    return 0;
}
