#include "gtest/gtest.h"
#include "World.hpp"
#include "GiNaCSymbolic.hpp"
#include "SEB.hpp"
#include "Exceptions.hpp"

// Global initialization of the SymbolicFactory
struct SymbolicFactoryInitializer {
    SymbolicFactoryInitializer() {
        GiNaCFactory* factory = new GiNaCFactory();
        SymbolicFactory::setInstance(factory);
    }
} g_initializer;

// Path to validation data files (relative to working directory, which is project root)
static const std::string VALIDATION_DIR = "Examples/Validation/";

// ===========================================================================================
// Validation_GaussianPolymer  (mirrors Examples/Validation_GaussianPolymer.cpp)
// ===========================================================================================
class GaussianPolymerValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new GaussianPolymer();
        w.Add(s, "poly");
        dir = VALIDATION_DIR + "GaussianPolymer/";
        w.setParameter(pl, "Rg_poly", 1.0);
    }
};

TEST_F(GaussianPolymerValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(GaussianPolymerValidation, ValidateSymbolically) {
    ASSERT_NO_THROW(s->ValidateSymbolically(pl));
}

TEST_F(GaussianPolymerValidation, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "F.dat"));
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_end1) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("end1", pl, dir + "Aend.dat"));
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_end2) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("end2", pl, dir + "Aend.dat"));
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_middle) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("middle", pl, dir + "Amiddle.dat"));
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_contour) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("contour", pl, dir + "F.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_end1_end2) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("end1", "end2", pl, dir + "Psi_end2end.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_end1_middle) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("end1", "middle", pl, dir + "Psi_middle2end.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_end2_middle) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("end2", "middle", pl, dir + "Psi_middle2end.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_end1) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "end1", pl, dir + "Aend.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_end2) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "end2", pl, dir + "Aend.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_middle) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "middle", pl, dir + "Psi_middle2contour.dat"));
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_contour) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "contour", pl, dir + "F.dat"));
}

// ===========================================================================================
// Validation_GaussianLoop  (mirrors Examples/Validation_GaussianLoop.cpp)
// ===========================================================================================
class GaussianLoopValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new GaussianLoop();
        w.Add(s, "loop");
        dir = VALIDATION_DIR + "GaussianLoop/";
        w.setParameter(pl, "Rg_loop", 1.0);
    }
};

TEST_F(GaussianLoopValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(GaussianLoopValidation, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "F.dat"));
}

TEST_F(GaussianLoopValidation, FormFactorAmplitude_contour) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("contour", pl, dir + "F.dat"));
}

TEST_F(GaussianLoopValidation, PhaseFactor_contour_contour) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "contour", pl, dir + "F.dat"));
}

// ===========================================================================================
// Validation_ThinRod  (mirrors Examples/Validation_ThinRod.cpp)
// ===========================================================================================
class ThinRodValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new ThinRod();
        w.Add(s, "rod");
        dir = VALIDATION_DIR + "ThinRod/";
        w.setParameter(pl, "L_rod", 1.0);
    }
};

TEST_F(ThinRodValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(ThinRodValidation, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "F.dat"));
}

TEST_F(ThinRodValidation, FormFactorAmplitude_contour) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("contour", pl, dir + "F.dat"));
}

TEST_F(ThinRodValidation, FormFactorAmplitude_end1) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("end1", pl, dir + "A_end.dat"));
}

TEST_F(ThinRodValidation, FormFactorAmplitude_end2) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("end2", pl, dir + "A_end.dat"));
}

TEST_F(ThinRodValidation, FormFactorAmplitude_middle) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("middle", pl, dir + "A_middle.dat"));
}

TEST_F(ThinRodValidation, PhaseFactor_end1_middle) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("end1", "middle", pl, dir + "Psi_end2middle.dat"));
}

TEST_F(ThinRodValidation, PhaseFactor_end2_middle) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("end2", "middle", pl, dir + "Psi_end2middle.dat"));
}

TEST_F(ThinRodValidation, PhaseFactor_contour_end1) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "end1", pl, dir + "A_end.dat"));
}

TEST_F(ThinRodValidation, PhaseFactor_contour_end2) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "end2", pl, dir + "A_end.dat"));
}

TEST_F(ThinRodValidation, PhaseFactor_contour_middle) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "middle", pl, dir + "Psi_contour2middle.dat"));
}

TEST_F(ThinRodValidation, PhaseFactor_contour_contour) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("contour", "contour", pl, dir + "F.dat"));
}

// ===========================================================================================
// Validation_ThinDisk  (mirrors Examples/Validation_ThinDisk.cpp)
// ===========================================================================================
class ThinDiskValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new ThinDisk();
        w.Add(s, "disk");
        dir = VALIDATION_DIR + "ThinDisk_R1/";
        w.setParameter(pl, "R_disk", 1.0);
    }
};

TEST_F(ThinDiskValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(ThinDiskValidation, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "FF.dat"));
}

TEST_F(ThinDiskValidation, FormFactorAmplitude_surface) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FF.dat"));
}

TEST_F(ThinDiskValidation, PhaseFactor_surface_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "FF.dat"));
}

TEST_F(ThinDiskValidation, FormFactorAmplitude_center) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"));
}

TEST_F(ThinDiskValidation, FormFactorAmplitude_rim) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("rim", pl, dir + "FFA_rim.dat"));
}

TEST_F(ThinDiskValidation, PhaseFactor_rim_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("rim", "surface", pl, dir + "FFA_rim.dat"));
}

TEST_F(ThinDiskValidation, PhaseFactor_center_rim) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "rim", pl, dir + "Psi_center2rim.dat"));
}

TEST_F(ThinDiskValidation, PhaseFactor_rim_rim) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("rim", "rim", pl, dir + "Psi_rim2rim.dat"));
}

// ===========================================================================================
// Validation_SolidSphere  (mirrors Examples/Validation_SolidSphere.cpp)
// ===========================================================================================
class SolidSphereValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;

    void SetUp() override {
        s = new SolidSphere();
        w.Add(s, "sph");
        pl = w.getParams();
    }
};

TEST_F(SolidSphereValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(SolidSphereValidation, ValidateSymbolically) {
    ASSERT_NO_THROW(s->ValidateSymbolically(pl));
}

// ===========================================================================================
// Validation_SolidCylinder  (mirrors Examples/Validation_SolidCylinder.cpp)
// ===========================================================================================

// --- Parameter set 1: R=1, L=1.5 ---
class SolidCylinderValidation_R1_L1p5 : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new SolidCylinder();
        w.Add(s, "c");
        dir = VALIDATION_DIR + "SolidCylinder_R1_L1.5/";
        w.setParameter(pl, "R_c", 1.0);
        w.setParameter(pl, "L_c", 1.5);
    }
};

TEST_F(SolidCylinderValidation_R1_L1p5, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "F.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_center) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_ends) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("ends", pl, dir + "FFA_ends.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_hull) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("hull", pl, dir + "FFA_hull.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_surface) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FFA_surface.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_center_ends) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "ends", pl, dir + "PF_center2ends.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_center_hull) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "hull", pl, dir + "PF_center2hull.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_center_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "surface", pl, dir + "PF_center2surface.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_ends_ends) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("ends", "ends", pl, dir + "PF_end2end.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_ends_hull) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("ends", "hull", pl, dir + "PF_end2hull.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_ends_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("ends", "surface", pl, dir + "PF_end2surface.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_hull_hull) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("hull", "hull", pl, dir + "PF_hull2hull.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_hull_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("hull", "surface", pl, dir + "PF_hull2surface.dat"));
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_surface_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "PF_surface2surface.dat"));
}

// --- Parameter set 2: R=2, L=0.5 ---
class SolidCylinderValidation_R2_L0p5 : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new SolidCylinder();
        w.Add(s, "c");
        dir = VALIDATION_DIR + "SolidCylinder_R2_L0.5/";
        w.setParameter(pl, "R_c", 2.0);
        w.setParameter(pl, "L_c", 0.5);
    }
};

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "F.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_center) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_ends) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("ends", pl, dir + "FFA_ends.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_hull) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("hull", pl, dir + "FFA_hull.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_surface) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FFA_surface.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_center_ends) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "ends", pl, dir + "PF_center2ends.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_center_hull) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "hull", pl, dir + "PF_center2hull.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_center_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "surface", pl, dir + "PF_center2surface.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_ends_ends) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("ends", "ends", pl, dir + "PF_end2end.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_ends_hull) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("ends", "hull", pl, dir + "PF_end2hull.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_ends_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("ends", "surface", pl, dir + "PF_end2surface.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_hull_hull) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("hull", "hull", pl, dir + "PF_hull2hull.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_hull_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("hull", "surface", pl, dir + "PF_hull2surface.dat"));
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_surface_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "PF_surface2surface.dat"));
}

// ===========================================================================================
// Validation_SolidSphericalShell  (mirrors Examples/Validation_SolidSphericalShell.cpp)
// ===========================================================================================
class SolidSphericalShellValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;
    std::string dir;

    void SetUp() override {
        s = new SolidSphericalShell();
        w.Add(s, "shell");
        dir = VALIDATION_DIR + "SolidSphericalShell_Ri2.33_Ro3.44/";
        w.setParameter(pl, "Ri_shell", 2.33);
        w.setParameter(pl, "Ro_shell", 3.44);
    }
};

TEST_F(SolidSphericalShellValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(SolidSphericalShellValidation, FormFactor) {
    EXPECT_TRUE(s->ValidateFormFactorFile(pl, dir + "FF.dat"));
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_center) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"));
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_surfacei) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("surfacei", pl, dir + "FFA_inner.dat"));
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_surfaceo) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("surfaceo", pl, dir + "FFA_outer.dat"));
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_surface) {
    EXPECT_TRUE(s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FFA_surface.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_center_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("center", "surface", pl, dir + "PF_center_surface.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surfacei_surfacei) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surfacei", "surfacei", pl, dir + "PF_inner_inner.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surfacei_surfaceo) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surfacei", "surfaceo", pl, dir + "PF_inner_outer.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surfaceo_surfaceo) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surfaceo", "surfaceo", pl, dir + "PF_outer_outer.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surface_surfacei) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surface", "surfacei", pl, dir + "PF_inner_surface.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surface_surfaceo) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surface", "surfaceo", pl, dir + "PF_outer_surface.dat"));
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surface_surface) {
    EXPECT_TRUE(s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "PF_surface_surface.dat"));
}

// ===========================================================================================
// Validation_ThinSphericalShell  (mirrors Examples/Validation_ThinSphericalShell.cpp)
// ===========================================================================================
class ThinSphericalShellValidation : public ::testing::Test {
protected:
    World w;
    SubUnit* s;
    ParameterList pl;

    void SetUp() override {
        s = new ThinSphericalShell();
        w.Add(s, "shell");
        pl = w.getParams();
    }
};

TEST_F(ThinSphericalShellValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(ThinSphericalShellValidation, ValidateSymbolically) {
    ASSERT_NO_THROW(s->ValidateSymbolically(pl));
}