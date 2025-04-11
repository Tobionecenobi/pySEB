#include "../Mathematica/Sampled/Sample_ThinDisk.cpp"
#include <gtest/gtest.h>
#include <cmath>
#include <fstream>

TEST(SampleThinDiskTest, FormFactorValidation) {
    const double R = 1.0;  // Disk radius
    Sample(R);
    
    // Check form factor normalization at q=0
    std::ifstream file("ThinDisk_R1/FF.q");
    ASSERT_TRUE(file.good());
    
    double q, F;
    file >> q >> F;
    EXPECT_NEAR(F, 1.0, 0.01);
}

TEST(SampleThinDiskTest, ContourCorrelationValidation) {
    const double R = 1.0;
    Sample(R);
    
    // Check contour-to-contour correlations
    std::ifstream file("ThinDisk_R1/PF_contour_contour.q");
    ASSERT_TRUE(file.good());
    
    double q, pf;
    file >> q >> pf;
    EXPECT_NEAR(pf, 1.0, 0.01);  // Should be 1 at q=0
    
    // For disk contour, mean square distance = R²
    while (file >> q >> pf) {
        if (q < 0.01) {
            double expected_pf = 1.0 - q*q*R*R/2.0;
            EXPECT_NEAR(pf, expected_pf, 0.01);
        }
    }
}

TEST(SampleThinDiskTest, CenterToContourValidation) {
    const double R = 1.0;
    Sample(R);
    
    // Check center-to-contour correlations
    std::ifstream file("ThinDisk_R1/PF_center_contour.q");
    ASSERT_TRUE(file.good());
    
    double q, pf;
    file >> q >> pf;
    EXPECT_NEAR(pf, 1.0, 0.01);  // Should be 1 at q=0
    
    // For center to contour, mean square distance = R²/2
    while (file >> q >> pf) {
        if (q < 0.01) {
            double expected_pf = 1.0 - q*q*R*R/4.0;
            EXPECT_NEAR(pf, expected_pf, 0.01);
        }
    }
}

TEST(SampleThinDiskTest, RadiusOfGyrationValidation) {
    const double R = 1.0;
    Sample(R);
    
    // For thin disk, Rg² = R²/2
    const double theoretical_Rg2 = R*R/2.0;
    
    std::ifstream file("ThinDisk_R1/FF.q");
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
