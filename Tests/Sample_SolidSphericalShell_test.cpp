#include "../Mathematica/Sampled/Sample_SolidSphericalShell.cpp"
#include <gtest/gtest.h>
#include <cmath>

TEST(SampleSolidSphericalShellTest, FormFactorValidation) {
    const double Ri = 2.0;
    const double Ro = 3.0;
    Sample(Ri, Ro);
    
    // Check form factor amplitude for shell
    std::ifstream file("SolidSphericalShell_Ri2_Ro3/FF.q");
    ASSERT_TRUE(file.good());
    
    std::string line;
    while (std::getline(file, line)) {
        double q, F;
        std::stringstream ss(line);
        ss >> q >> F;
        
        if (q > 0.1) continue; // Only check small q regime
        
        // Form factor should approach 1 at q→0
        if (q < 0.01) {
            EXPECT_NEAR(F, 1.0, 0.01);
        }
    }
}

TEST(SampleSolidSphericalShellTest, InnerOuterCorrelationValidation) {
    const double Ri = 2.0;
    const double Ro = 3.0;
    Sample(Ri, Ro);
    
    // Validate inner-outer surface correlations
    std::ifstream file("SolidSphericalShell_Ri2_Ro3/Pinner_outer.q");
    ASSERT_TRUE(file.good());
    
    std::string line;
    while (std::getline(file, line)) {
        double q, pf;
        std::stringstream ss(line);
        ss >> q >> pf;
        
        if (q > 0.1) continue;
        
        // For q→0, expect correlation to approach 1
        if (q < 0.01) {
            EXPECT_NEAR(pf, 1.0, 0.01);
        }
    }
}

TEST(SampleSolidSphericalShellTest, RadiusOfGyrationValidation) {
    const double Ri = 2.0;
    const double Ro = 3.0;
    Sample(Ri, Ro);
    
    // Theoretical Rg^2 for a spherical shell
    double theoretical_Rg2 = (3.0/5.0) * 
        (Ro*pow(Ri,3) + pow(Ri,4) + pow(Ri,2)*pow(Ro,2) + 
         Ri*pow(Ro,3) + pow(Ro,4))/(Ri*Ro + pow(Ri,2) + pow(Ro,2));
    
    // Read form factor at small q and compare to Guinier approximation
    std::ifstream file("SolidSphericalShell_Ri2_Ro3/FF.q");
    ASSERT_TRUE(file.good());
    
    // Allow 1% error due to numerical sampling
    const double tolerance = 0.01;
    
    double q, F;
    while (file >> q >> F) {
        if (q < 0.01) {
            // In Guinier regime: F(q) ≈ 1 - q²Rg²/3
            double expected_F = 1.0 - q*q*theoretical_Rg2/3.0;
            EXPECT_NEAR(F, expected_F, tolerance);
        }
    }
}
