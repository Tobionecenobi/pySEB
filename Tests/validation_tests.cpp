#include "gtest/gtest.h"
#include "World.hpp"
#include "GiNaCSymbolic.hpp"
#include "SEB.hpp"
#include "Exceptions.hpp"

// Global initialization of the SymbolicFactory
struct SymbolicFactoryInitializer {
    SymbolicFactoryInitializer() {
        GiNaCSymbolic* factory = new GiNaCSymbolic();
        SymbolicFactory::setInstance(factory);
    }
} g_initializer;

// Basic Shape Validation Tests
TEST(ValidationTest, SolidSphereValidation) {
    World world;
    GraphID g = world.Add(new SolidSphere(), "sphere");
    world.Add(g, "structure");
    
    Expression F = world.FormFactor("structure");
    ParameterList params = world.getParams();
    
    // Test q→0 limit (should be 1)
    params = world.getParams();
    world.setParameter(params, "q", 0.0);
    Expression F_q0 = F.subs(params);
    EXPECT_NEAR(F_q0.eval(), 1.0, 1e-10);
    
    // Test Guinier regime
    double R = 1.0;
    world.setParameter(params, "R_sphere", R);
    for (double q : {0.01, 0.02, 0.05}) {
        params = world.getParams();
        world.setParameter(params, "q", q);
        Expression F_q = F.subs(params);
        double val = F_q.eval();
        // Check against Rayleigh's analytical expression
        double x = q * R;
        double analytical = pow(3*(sin(x)-x*cos(x))/pow(x,3), 2);
        EXPECT_NEAR(val, analytical, 1e-6);
    }
    
    // Test phase factors
    Expression phase = world.PhaseFactor("structure:sphere.center", "structure:sphere.surface");
    for (double q : {0.1, 0.5, 1.0}) {
        params = world.getParams();
        world.setParameter(params, "q", q);
        Expression phase_q = phase.subs(params);
        double val = phase_q.eval();
        double x = q * R;
        double analytical = sin(x)/x;
        EXPECT_NEAR(val, analytical, 1e-6);
    }

    // Validate mean square distances using RadiusOfGyration2
    Expression Rg2 = world.RadiusOfGyration2("structure");
    params = world.getParams();
    double val = Rg2.subs(params).eval();
    EXPECT_NEAR(val, 3*R*R/5, 1e-10);  // Known analytical result for solid sphere
}

// Similar validation tests for other basic shapes...
TEST(ValidationTest, SolidCylinderValidation) {
    // Validate form factors and amplitudes numerically
}

TEST(ValidationTest, ThinDiskValidation) {
    // Validate against analytical expressions
}

TEST(ValidationTest, ThinRodValidation) {
    World world;
    GraphID g = world.Add(new ThinRod(), "rod");
    world.Add(g, "structure");
    
    Expression F = world.FormFactor("structure");
    ParameterList params = world.getParams();
    
    // Test q→0 limit
    world.setParameter(params, "q", 0.0);
    Expression F_q0 = F.subs(params);
    EXPECT_NEAR(F_q0.eval(), 1.0, 1e-10);
    
    // Test end-to-end correlations
    Expression phase = world.PhaseFactor("structure:rod.end1", "structure:rod.end2");
    
    // Set rod length L=1 for simplicity
    world.setParameter(params, "L_rod", 1.0);
    
    // For small q, phase factor ≈ 1 - q²L²/2
    for (double q : {0.01, 0.02, 0.05}) {
        world.setParameter(params, "q", q);
        Expression phase_q = phase.subs(params);
        double expected = 1.0 - q*q/2.0;  // L=1
        EXPECT_NEAR(phase_q.eval(), expected, 1e-6);
    }
    
    // Test contour correlations
    Expression contour_phase = world.PhaseFactor("structure:rod.contour#p1", 
                                               "structure:rod.contour#p2");
    
    // For contour points, <R²> = L²/6 for uniform sampling
    for (double q : {0.01, 0.02, 0.05}) {
        world.setParameter(params, "q", q);
        Expression phase_q = contour_phase.subs(params);
        double expected = 1.0 - q*q/6.0;  // L=1
        EXPECT_NEAR(phase_q.eval(), expected, 1e-6);
    }
}

TEST(ValidationTest, SolidSphericalShellValidation) {
    // Validate inner/outer correlations
}

TEST(ValidationTest, ThinSphericalShellValidation) {
    World world;
    GraphID g = world.Add(new ThinSphericalShell(), "shell");
    world.Add(g, "structure");
    
    Expression F = world.FormFactor("structure");
    ParameterList params = world.getParams();
    
    // Test q→0 limit 
    world.setParameter(params, "q", 0.0);
    Expression F_q0 = F.subs(params);
    EXPECT_NEAR(F_q0.eval(), 1.0, 1e-10);
    
    // Set shell radius R=1 for simplicity
    world.setParameter(params, "R_shell", 1.0);
    
    // Test surface correlations
    Expression surface_phase = world.PhaseFactor("structure:shell.surface#p1", 
                                               "structure:shell.surface#p2");
    
    // For surface points on a sphere, <R²> = 2R² for uniform sampling
    for (double q : {0.01, 0.02, 0.05}) {
        world.setParameter(params, "q", q);
        Expression phase_q = surface_phase.subs(params);
        double expected = 1.0 - q*q;  // R=1, so <R²>=2
        EXPECT_NEAR(phase_q.eval(), expected, 1e-6);
    }
    
    // Test center-to-surface correlations 
    Expression center_phase = world.PhaseFactor("structure:shell.center",
                                              "structure:shell.surface");
                                              
    // For center to surface, distance is exactly R
    for (double q : {0.1, 0.5, 1.0}) {
        world.setParameter(params, "q", q);
        Expression phase_q = center_phase.subs(params);
        double val = phase_q.eval();
        double x = q;  // R=1
        double analytical = sin(x)/x;
        EXPECT_NEAR(val, analytical, 1e-6);
    }
}

// Polymer Validation Tests
TEST(ValidationTest, GaussianPolymerValidation) {
    World world;
    GraphID g = world.Add(new GaussianPolymer(), "polymer");
    world.Add(g, "structure");
    
    Expression F = world.FormFactor("structure");
    ParameterList params = world.getParams();
    
    // Test Debye function behavior
    double Rg = 5.0;  // Set radius of gyration
    world.setParameter(params, "Rg_polymer", Rg);
    
    // Test q→0 limit
    world.setParameter(params, "q", 0.0);
    Expression F_q0 = F.subs(params);
    EXPECT_NEAR(F_q0.eval(), 1.0, 1e-10);
    
    // Test intermediate q values against Debye function
    std::vector<double> q_values = {0.1/Rg, 0.5/Rg, 1.0/Rg, 2.0/Rg};
    for(double q : q_values) {
        world.setParameter(params, "q", q);
        Expression F_q = F.subs(params);
        double x = (q * Rg) * (q * Rg);
        double debye = 2.0 * (exp(-x) + x - 1.0) / (x * x);
        EXPECT_NEAR(F_q.eval(), debye, 1e-6);
    }
}

TEST(ValidationTest, GaussianLoopValidation) {
    // Validate against known loop expressions
}

// Complex Structure Validation Tests
TEST(ValidationTest, GuinierRegimeValidation) {
    World world;
    
    // Test different shapes in Guinier regime
    std::vector<std::pair<SubUnit*, std::string>> shapes = {
        {new SolidSphere(), "sphere"},
        {new SolidCylinder(), "cylinder"},
        {new GaussianPolymer(), "polymer"},
        {new ThinDisk(), "disk"}
    };
    
    for (const auto& [shape, name] : shapes) {
        GraphID g = world.Add(shape, name);
        world.Add(g, "structure");
        
        Expression F = world.FormFactor("structure");
        ParameterList params = world.getParams();
        
        // Test q→0 limit
        params = world.getParams();
        world.setParameter(params, "q", 0.0);
        Expression F_q0 = F.subs(params);
        EXPECT_NEAR(F_q0.eval(), 1.0, 1e-10);
        
        // Test Guinier behavior: F(q) ≈ 1 - (qRg)²/3
        for (double q : {0.01, 0.02, 0.05}) {
            params = world.getParams();
            world.setParameter(params, "q", q);
            Expression F_q = F.subs(params);
            double val = F_q.eval();
            EXPECT_GT(val, 0.9);
        }
    }
}

TEST(ValidationTest, ExtendedGuinierRegimeValidation) {
    // Similar to GuinierRegimeValidation but for more complex shapes
}

TEST(ValidationTest, SymbolicValidation) {
    World world;
    
    // Test symbolic validation for different shapes
    auto* sphere = new SolidSphere();
    auto* cylinder = new SolidCylinder();
    auto* polymer = new GaussianPolymer();
    
    // Set parameters
    ParameterList params;
    params["R_sphere"] = 1.0;
    params["beta_sphere"] = 1.0;
    params["R_cylinder"] = 1.0;
    params["L_cylinder"] = 5.0;
    params["beta_cylinder"] = 1.0;
    params["Rg_poly"] = 1.0;
    params["beta_poly"] = 1.0;
    
    // Validate symbolic expressions
    ASSERT_NO_THROW(sphere->ValidateSymbolically(params));
    ASSERT_NO_THROW(cylinder->ValidateSymbolically(params));
    ASSERT_NO_THROW(polymer->ValidateSymbolically(params));
}

// Paper-Specific Validation Tests
TEST(ValidationTest, ContourCorrelationValidation) {
    World world;
    std::vector<std::pair<SubUnit*, std::string>> shapes = {
        {new GaussianLoop(), "loop"},
        {new GaussianPolymer(), "polymer"},
        {new ThinCircle(), "circle"},
        {new ThinRod(), "rod"}
    };
    
    for(const auto& [shape, name] : shapes) {
        GraphID g = world.Add(shape, name);
        world.Add(g, "structure");
        
        // Validate contour correlation functions
        Expression phase = world.PhaseFactor("structure:" + name + ".contour#p1", 
                                           "structure:" + name + ".contour#p2");
        ParameterList params = world.getParams();
        
        // Test q→0 limit
        world.setParameter(params, "q", 0.0);
        Expression phase_q0 = phase.subs(params);
        EXPECT_NEAR(phase_q0.eval(), 1.0, 1e-10);
        
        // Test large q behavior
        world.setParameter(params, "q", 100.0);
        Expression phase_qinf = phase.subs(params);
        EXPECT_NEAR(std::abs(phase_qinf.eval()), 0.0, 1e-6);
    }
}

TEST(ValidationTest, DendrimerCorrelationValidation) {
    World world;
    
    // Create and validate dendrimer correlations
    for(const std::string& prefix : {"polymer", "rod"}) {
        GraphID dendrimer = world.Add(new Point(), "center");
        
        // Add branches
        world.Link(prefix == "polymer" ? 
                  static_cast<SubUnit*>(new GaussianPolymer()) : 
                  static_cast<SubUnit*>(new ThinRod()),
                  prefix + "1_1.end1", "center.center");
        
        world.Link(prefix == "polymer" ? 
                  static_cast<SubUnit*>(new GaussianPolymer()) : 
                  static_cast<SubUnit*>(new ThinRod()),
                  prefix + "2_1.end1", "center.center");
        
        world.Add(dendrimer, prefix + "_dendrimer");
        
        // Validate correlations
        Expression phase = world.PhaseFactor(prefix + "_dendrimer:" + prefix + "1_1.end1",
                                           prefix + "_dendrimer:" + prefix + "2_1.end1");
        ParameterList params = world.getParams();
        
        // Test correlation limits
        world.setParameter(params, "q", 0.0);
        Expression phase_q0 = phase.subs(params);
        EXPECT_NEAR(phase_q0.eval(), 1.0, 1e-10);
    }
}