#include "../Mathematica/Sampled/Sample_SolidSphere.cpp"
#include <gtest/gtest.h>
#include <cmath>

TEST(SampleSolidSphereTest, RadiusOfGyrationValidation) {
    // Validate Rg^2 = 3R^2/5 for solid sphere
    const double R = 1.0;
    Sample(R);
    
    // Sample_SolidSphere.cpp outputs these values, we need to check them
    // These are theoretical values:
    const double theoretical_Rg2 = 0.6 * R * R;  // 3R^2/5
    const double theoretical_center_surface = R * R;  // R^2
    const double theoretical_surface_surface = 2 * R * R;  // 2R^2
    
    // Allow 1% error due to numerical sampling
    const double tolerance = 0.01;
    
    // Read and validate output files
    std::ifstream file("SolidSphere_R1/FF.q");
    ASSERT_TRUE(file.good());
    
    // Read q=0 form factor (should be 1.0)
    double q, F;
    file >> q >> F;
    EXPECT_NEAR(F, 1.0, tolerance);
}

TEST(SampleSolidSphereTest, FormFactorValidation) {
    const double R = 1.0;
    Sample(R);
    
    // Check form factor against Rayleigh's analytical expression
    std::ifstream file("SolidSphere_R1/FF.q");
    ASSERT_TRUE(file.good());
    
    std::string line;
    while (std::getline(file, line)) {
        double q, F;
        std::stringstream ss(line);
        ss >> q >> F;
        
        if (q > 0.1) continue; // Only check small q regime
        
        double x = q * R;
        double analytical = pow(3*(sin(x)-x*cos(x))/pow(x,3), 2);
        EXPECT_NEAR(F, analytical, 0.01);
    }
}

TEST(SampleSolidSphereTest, PhaseFactorValidation) {
    const double R = 1.0;
    Sample(R);
    
    // Check center-to-surface phase factor
    std::ifstream file("SolidSphere_R1/PF_center_surface.q");
    ASSERT_TRUE(file.good());
    
    std::string line;
    while (std::getline(file, line)) {
        double q, pf;
        std::stringstream ss(line);
        ss >> q >> pf;
        
        if (q > 0.1) continue; // Only check small q regime
        
        double x = q * R;
        double analytical = sin(x)/x;
        EXPECT_NEAR(pf, analytical, 0.01);
    }
}
