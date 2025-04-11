#include "../Mathematica/Sampled/Sample_SolidCylinder.cpp"
#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

TEST(SampleSolidCylinderTest, FormFactorValidation) {
    const double R = 1.0;
    const double L = 5.0;
    Sample(R, L);
    
    // Check form factor normalization
    std::ifstream file("SolidCylinder_R1_L5/FF.q");
    ASSERT_TRUE(file.good());
    
    // Form factor should be 1 at q=0
    double q, F;
    file >> q >> F;
    EXPECT_NEAR(F, 1.0, 0.01);
}

TEST(SampleSolidCylinderTest, RadiusOfGyrationValidation) {
    const double R = 1.0;
    const double L = 5.0;
    Sample(R, L);
    
    // Known Rg^2 = (L^2 + 6R^2)/12 for solid cylinder
    const double theoretical_Rg2 = (L*L + 6*R*R)/12.0;
    
    // Check Guinier behavior
    std::ifstream file("SolidCylinder_R1_L5/FF.q");
    ASSERT_TRUE(file.good());
    
    std::string line;
    while (std::getline(file, line)) {
        double q, F;
        std::stringstream ss(line);
        ss >> q >> F;
        
        if (q < 0.01) {
            // In Guinier regime: F(q) ≈ 1 - q²Rg²/3
            double expected_F = 1.0 - q*q*theoretical_Rg2/3.0;
            EXPECT_NEAR(F, expected_F, 0.01);
        }
    }
}

TEST(SampleSolidCylinderTest, SurfaceCorrelationsValidation) {
    const double R = 1.0;
    const double L = 5.0;
    Sample(R, L);
    
    // Check hull-to-hull and end-to-end correlations
    std::ifstream hull_file("SolidCylinder_R1_L5/PF_hull_hull.q");
    std::ifstream end_file("SolidCylinder_R1_L5/PF_end_end.q");
    ASSERT_TRUE(hull_file.good());
    ASSERT_TRUE(end_file.good());
    
    // Both correlations should approach 1 at q=0
    double q, pf;
    hull_file >> q >> pf;
    EXPECT_NEAR(pf, 1.0, 0.01);
    
    end_file >> q >> pf;
    EXPECT_NEAR(pf, 1.0, 0.01);
    
    // Check mean square distances
    const double hull_hull_msd = L*L/6.0 + 2*R*R;  // Known analytical result
    const double end_end_msd = L*L/2.0 + R*R;      // Known analytical result
    
    // For small q, phase factor ≈ 1 - q²<R²>/2
    while (hull_file >> q >> pf && q < 0.01) {
        double expected_pf = 1.0 - q*q*hull_hull_msd/2.0;
        EXPECT_NEAR(pf, expected_pf, 0.01);
    }
    
    while (end_file >> q >> pf && q < 0.01) {
        double expected_pf = 1.0 - q*q*end_end_msd/2.0;
        EXPECT_NEAR(pf, expected_pf, 0.01);
    }
}
