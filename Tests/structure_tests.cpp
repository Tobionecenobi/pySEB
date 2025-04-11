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

// Test Complex Structure Creation
TEST(StructureTest, CrossTermStructure) {
    World world;
    
    // Test different combinations of subunits
    GraphID g1 = world.Add(new SolidSphere(), "sphere");
    world.Link(new GaussianPolymer(), "poly.end1", "sphere.surface#p1", "poly");
    world.Link(new ThinDisk(), "disk.center", "sphere.surface#p2");
    world.Add(g1, "structure");
    
    EXPECT_TRUE(world.hasName("structure"));
    ASSERT_NO_THROW(world.FormFactor("structure"));
}

TEST(StructureTest, ShellComplexStructure) {
    World world;
    
    // Create a structure with shells and other shapes
    GraphID g1 = world.Add(new SolidSphericalShell(), "solid_shell");
    world.Link(new ThinRod(), "rod.end1", "solid_shell.outer#p1", "rod");
    world.Link(new ThinSphericalShell(), "thin_shell.center", "solid_shell.outer#p2");
    world.Add(g1, "structure");
    
    EXPECT_TRUE(world.hasName("structure"));
    ASSERT_NO_THROW(world.FormFactor("structure"));
}

TEST(StructureTest, DistributedPointStructure) {
    World world;
    
    // Test different shapes with distributed points
    std::vector<std::pair<SubUnit*, std::string>> shapes = {
        {new GaussianPolymer(), "contour"},
        {new ThinDisk(), "contour"},
        {new SolidSphere(), "surface"}
    };
    
    for (const auto& [shape, name] : shapes) {
        GraphID g = world.Add(shape, name);
        world.Add(g, "structure");
        EXPECT_TRUE(world.hasName("structure"));
    }
}

TEST(StructureTest, ChainOfRodsStructure) {
    World world;
    GraphID chain = world.Add(new ThinRod(), "rod1");
    
    for(int i = 2; i <= 5; i++) {
        world.Link(new ThinRod(), "rod" + std::to_string(i) + ".end1", 
                  "rod" + std::to_string(i-1) + ".end2");
    }
    world.Add(chain, "rod_chain");
    
    EXPECT_TRUE(world.hasName("rod_chain"));
    ASSERT_NO_THROW(world.FormFactor("rod_chain"));
}

// Existing DendrimerCreation test updated with more structure checks
TEST(StructureTest, DendrimerCreation) {
    World world;
    
    // Create a simple dendrimer with functionality f=3 and generation g=2
    const int f = 3;
    const int g = 2;
    
    // Create core
    GraphID dendrimer = world.Add(new Point(), "core");
    
    // First generation
    for(int i = 1; i <= f; i++) {
        world.Link(new GaussianPolymer(), "poly1_" + std::to_string(i) + ".end1", "core.center");
        
        // Second generation
        for(int j = 1; j <= f-1; j++) {
            world.Link(new GaussianPolymer(), 
                      "poly2_" + std::to_string(i) + "_" + std::to_string(j) + ".end1",
                      "poly1_" + std::to_string(i) + ".end2");
        }
    }
    world.Add(dendrimer, "dendrimer");
    
    EXPECT_TRUE(world.hasName("dendrimer"));
    ASSERT_NO_THROW(world.FormFactor("dendrimer"));
}

TEST(StructureTest, DiBlockStarChainStructure) {
    World world;
    
    // Create a diblock star chain
    GraphID chain = world.Add(new Point(), "center1");
    
    // Create first star
    for(int i = 1; i <= 3; i++) {
        world.Link(new GaussianPolymer(), "poly1_" + std::to_string(i) + ".end1", "center1.center");
    }
    
    // Add second star
    world.Link(new Point(), "center2.center", "poly1_1.end2");
    for(int i = 2; i <= 3; i++) {
        world.Link(new GaussianPolymer(), "poly2_" + std::to_string(i) + ".end1", "center2.center");
    }
    
    world.Add(chain, "star_chain");
    EXPECT_TRUE(world.hasName("star_chain"));
    ASSERT_NO_THROW(world.FormFactor("star_chain"));
}

TEST(StructureTest, MicelleStructure) {
    World world;
    
    // Create micelle with N polymers attached to spherical core
    const int N = 8;
    GraphID micelle = world.Add(new SolidSphere(), "core");
    
    for(int i = 1; i <= N; i++) {
        world.Link(new GaussianPolymer(), "poly" + std::to_string(i) + ".end1", 
                  "core.surface#p" + std::to_string(i));
    }
    world.Add(micelle, "micelle");
    
    EXPECT_TRUE(world.hasName("micelle"));
    ASSERT_NO_THROW(world.FormFactor("micelle"));
}

TEST(StructureTest, StarPolymerStructure) {
    World world;
    
    // Create star with f arms
    const int f = 6;
    GraphID star = world.Add(new Point(), "center");
    
    for(int i = 1; i <= f; i++) {
        world.Link(new GaussianPolymer(), "arm" + std::to_string(i) + ".end1", "center.center");
    }
    world.Add(star, "star");
    
    EXPECT_TRUE(world.hasName("star"));
    ASSERT_NO_THROW(world.FormFactor("star"));
}

TEST(StructureTest, TriBlockCopolymerStructure) {
    World world;
    
    // Create ABC triblock copolymer
    GraphID chain = world.Add(new GaussianPolymer(), "blockA");
    world.Link(new GaussianPolymer(), "blockB.end1", "blockA.end2");
    world.Link(new GaussianPolymer(), "blockC.end1", "blockB.end2");
    world.Add(chain, "triblock");
    
    EXPECT_TRUE(world.hasName("triblock"));
    ASSERT_NO_THROW(world.FormFactor("triblock"));
}