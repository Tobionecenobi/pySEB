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

#ifndef VALIDATION_DIR
#define VALIDATION_DIR "Examples/Validation/"
#endif

static const std::string VALIDATION_DATA_DIR = VALIDATION_DIR;

template <typename Fn>
void RunValidationDiagnostic(Fn fn) {
    try {
        (void)fn();
    } catch (const SEBException& e) {
        std::cout << "WARNING validation diagnostic threw:\n" << e << std::endl;
    }
}

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
        dir = VALIDATION_DATA_DIR + "GaussianPolymer/";
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
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "F.dat"); });
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_end1) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("end1", pl, dir + "Aend.dat"); });
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_end2) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("end2", pl, dir + "Aend.dat"); });
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_middle) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("middle", pl, dir + "Amiddle.dat"); });
}

TEST_F(GaussianPolymerValidation, FormFactorAmplitude_contour) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("contour", pl, dir + "F.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_end1_end2) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("end1", "end2", pl, dir + "Psi_end2end.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_end1_middle) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("end1", "middle", pl, dir + "Psi_middle2end.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_end2_middle) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("end2", "middle", pl, dir + "Psi_middle2end.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_end1) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "end1", pl, dir + "Aend.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_end2) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "end2", pl, dir + "Aend.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_middle) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "middle", pl, dir + "Psi_middle2contour.dat"); });
}

TEST_F(GaussianPolymerValidation, PhaseFactor_contour_contour) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "contour", pl, dir + "F.dat"); });
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
        dir = VALIDATION_DATA_DIR + "GaussianLoop/";
        w.setParameter(pl, "Rg_loop", 1.0);
    }
};

TEST_F(GaussianLoopValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(GaussianLoopValidation, FormFactor) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "F.dat"); });
}

TEST_F(GaussianLoopValidation, FormFactorAmplitude_contour) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("contour", pl, dir + "F.dat"); });
}

TEST_F(GaussianLoopValidation, PhaseFactor_contour_contour) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "contour", pl, dir + "F.dat"); });
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
        dir = VALIDATION_DATA_DIR + "ThinRod/";
        w.setParameter(pl, "L_rod", 1.0);
    }
};

TEST_F(ThinRodValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(ThinRodValidation, FormFactor) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "F.dat"); });
}

TEST_F(ThinRodValidation, FormFactorAmplitude_contour) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("contour", pl, dir + "F.dat"); });
}

TEST_F(ThinRodValidation, FormFactorAmplitude_end1) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("end1", pl, dir + "A_end.dat"); });
}

TEST_F(ThinRodValidation, FormFactorAmplitude_end2) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("end2", pl, dir + "A_end.dat"); });
}

TEST_F(ThinRodValidation, FormFactorAmplitude_middle) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("middle", pl, dir + "A_middle.dat"); });
}

TEST_F(ThinRodValidation, PhaseFactor_end1_middle) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("end1", "middle", pl, dir + "Psi_end2middle.dat"); });
}

TEST_F(ThinRodValidation, PhaseFactor_end2_middle) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("end2", "middle", pl, dir + "Psi_end2middle.dat"); });
}

TEST_F(ThinRodValidation, PhaseFactor_contour_end1) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "end1", pl, dir + "A_end.dat"); });
}

TEST_F(ThinRodValidation, PhaseFactor_contour_end2) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "end2", pl, dir + "A_end.dat"); });
}

TEST_F(ThinRodValidation, PhaseFactor_contour_middle) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "middle", pl, dir + "Psi_contour2middle.dat"); });
}

TEST_F(ThinRodValidation, PhaseFactor_contour_contour) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("contour", "contour", pl, dir + "F.dat"); });
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
        dir = VALIDATION_DATA_DIR + "ThinDisk_R1/";
        w.setParameter(pl, "R_disk", 1.0);
    }
};

TEST_F(ThinDiskValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(ThinDiskValidation, FormFactor) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "FF.dat"); });
}

TEST_F(ThinDiskValidation, FormFactorAmplitude_surface) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FF.dat"); });
}

TEST_F(ThinDiskValidation, PhaseFactor_surface_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "FF.dat"); });
}

TEST_F(ThinDiskValidation, FormFactorAmplitude_center) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"); });
}

TEST_F(ThinDiskValidation, FormFactorAmplitude_rim) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("rim", pl, dir + "FFA_rim.dat"); });
}

TEST_F(ThinDiskValidation, PhaseFactor_rim_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("rim", "surface", pl, dir + "FFA_rim.dat"); });
}

TEST_F(ThinDiskValidation, PhaseFactor_center_rim) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "rim", pl, dir + "Psi_center2rim.dat"); });
}

TEST_F(ThinDiskValidation, PhaseFactor_rim_rim) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("rim", "rim", pl, dir + "Psi_rim2rim.dat"); });
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
        dir = VALIDATION_DATA_DIR + "SolidCylinder_R1_L1.5/";
        w.setParameter(pl, "R_c", 1.0);
        w.setParameter(pl, "L_c", 1.5);
    }
};

TEST_F(SolidCylinderValidation_R1_L1p5, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactor) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "F.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_center) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_ends) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("ends", pl, dir + "FFA_ends.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_hull) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("hull", pl, dir + "FFA_hull.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, FormFactorAmplitude_surface) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FFA_surface.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_center_ends) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "ends", pl, dir + "PF_center2ends.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_center_hull) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "hull", pl, dir + "PF_center2hull.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_center_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "surface", pl, dir + "PF_center2surface.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_ends_ends) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("ends", "ends", pl, dir + "PF_end2end.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_ends_hull) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("ends", "hull", pl, dir + "PF_end2hull.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_ends_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("ends", "surface", pl, dir + "PF_end2surface.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_hull_hull) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("hull", "hull", pl, dir + "PF_hull2hull.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_hull_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("hull", "surface", pl, dir + "PF_hull2surface.dat"); });
}

TEST_F(SolidCylinderValidation_R1_L1p5, PhaseFactor_surface_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "PF_surface2surface.dat"); });
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
        dir = VALIDATION_DATA_DIR + "SolidCylinder_R2_L0.5/";
        w.setParameter(pl, "R_c", 2.0);
        w.setParameter(pl, "L_c", 0.5);
    }
};

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactor) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "F.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_center) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_ends) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("ends", pl, dir + "FFA_ends.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_hull) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("hull", pl, dir + "FFA_hull.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, FormFactorAmplitude_surface) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FFA_surface.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_center_ends) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "ends", pl, dir + "PF_center2ends.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_center_hull) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "hull", pl, dir + "PF_center2hull.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_center_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "surface", pl, dir + "PF_center2surface.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_ends_ends) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("ends", "ends", pl, dir + "PF_end2end.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_ends_hull) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("ends", "hull", pl, dir + "PF_end2hull.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_ends_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("ends", "surface", pl, dir + "PF_end2surface.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_hull_hull) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("hull", "hull", pl, dir + "PF_hull2hull.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_hull_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("hull", "surface", pl, dir + "PF_hull2surface.dat"); });
}

TEST_F(SolidCylinderValidation_R2_L0p5, PhaseFactor_surface_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "PF_surface2surface.dat"); });
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
        dir = VALIDATION_DATA_DIR + "SolidSphericalShell_Ri2.33_Ro3.44/";
        w.setParameter(pl, "Ri_shell", 2.33);
        w.setParameter(pl, "Ro_shell", 3.44);
    }
};

TEST_F(SolidSphericalShellValidation, ValidateDefinedTerms) {
    ASSERT_NO_THROW(s->ValidateDefinedTerms());
}

TEST_F(SolidSphericalShellValidation, FormFactor) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorFile(pl, dir + "FF.dat"); });
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_center) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("center", pl, dir + "FFA_center.dat"); });
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_surfacei) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("surfacei", pl, dir + "FFA_inner.dat"); });
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_surfaceo) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("surfaceo", pl, dir + "FFA_outer.dat"); });
}

TEST_F(SolidSphericalShellValidation, FormFactorAmplitude_surface) {
    RunValidationDiagnostic([&]() { return s->ValidateFormFactorAmplitudeFile("surface", pl, dir + "FFA_surface.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_center_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("center", "surface", pl, dir + "PF_center_surface.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surfacei_surfacei) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surfacei", "surfacei", pl, dir + "PF_inner_inner.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surfacei_surfaceo) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surfacei", "surfaceo", pl, dir + "PF_inner_outer.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surfaceo_surfaceo) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surfaceo", "surfaceo", pl, dir + "PF_outer_outer.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surface_surfacei) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surface", "surfacei", pl, dir + "PF_inner_surface.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surface_surfaceo) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surface", "surfaceo", pl, dir + "PF_outer_surface.dat"); });
}

TEST_F(SolidSphericalShellValidation, PhaseFactor_surface_surface) {
    RunValidationDiagnostic([&]() { return s->ValidatePhaseFactorFile("surface", "surface", pl, dir + "PF_surface_surface.dat"); });
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
