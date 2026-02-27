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
    GraphID g1 = world.Add(new SolidSphericalShell(), "solidShell");
    world.Link(new ThinRod(), "rod.end1", "solidShell.surfaceo#p1", "rod");
    world.Link(new ThinSphericalShell(), "thinShell.center", "solidShell.surfaceo#p2");
    world.Add(g1, "structure");
    
    EXPECT_TRUE(world.hasName("structure"));
    ASSERT_NO_THROW(world.FormFactor("structure"));
}

TEST(StructureTest, DistributedPointStructure) {
    // Each shape tested independently in its own world (names like "contour" would clash if shared)
    {
        World w;
        GraphID g = w.Add(new GaussianPolymer(), "polymer");
        w.Add(g, "polymerStruct");
        EXPECT_TRUE(w.hasName("polymerStruct"));
        ASSERT_NO_THROW(w.FormFactor("polymerStruct"));
    }
    {
        World w;
        GraphID g = w.Add(new ThinDisk(), "disk");
        w.Add(g, "diskStruct");
        EXPECT_TRUE(w.hasName("diskStruct"));
        ASSERT_NO_THROW(w.FormFactor("diskStruct"));
    }
    {
        World w;
        GraphID g = w.Add(new SolidSphere(), "sphere");
        w.Add(g, "sphereStruct");
        EXPECT_TRUE(w.hasName("sphereStruct"));
        ASSERT_NO_THROW(w.FormFactor("sphereStruct"));
    }
}

TEST(StructureTest, ChainOfRodsStructure) {
    World world;
    GraphID chain = world.Add(new ThinRod(), "rod1");
    
    for(int i = 2; i <= 5; i++) {
        world.Link(new ThinRod(), "rod" + std::to_string(i) + ".end1", 
                  "rod" + std::to_string(i-1) + ".end2");
    }
    world.Add(chain, "rodChain");
    
    EXPECT_TRUE(world.hasName("rodChain"));
    ASSERT_NO_THROW(world.FormFactor("rodChain"));
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
        world.Link(new GaussianPolymer(), "poly1" + std::to_string(i) + ".end1", "core.point");
        
        // Second generation
        for(int j = 1; j <= f-1; j++) {
            world.Link(new GaussianPolymer(), 
                      "poly2" + std::to_string(i) + std::to_string(j) + ".end1",
                      "poly1" + std::to_string(i) + ".end2");
        }
    }
    world.Add(dendrimer, "dendrimer");
    
    EXPECT_TRUE(world.hasName("dendrimer"));
    ASSERT_NO_THROW(world.FormFactor("dendrimer"));
}

TEST(StructureTest, DiBlockStarChainStructure) {
    World world;

    // Define diblock copolymer
    GraphID diblock = world.Add(new GaussianPolymer(), "polyA");
    world.Link(new GaussianPolymer(), "polyB.end1", "polyA.end2");

    // Define star
    GraphID star = world.Add(diblock, "diblock1");
    world.Link(diblock, "diblock2:polyA.end1", "diblock1:polyA.end1");
    world.Link(diblock, "diblock3:polyA.end1", "diblock1:polyA.end1");
    world.Link(diblock, "diblock4:polyA.end1", "diblock1:polyA.end1");

    // Define chain of stars
    GraphID chain = world.Add(star, "star1");
    world.Link(star, "star2:diblock1:polyB.end2", "star1:diblock3:polyB.end2");
    world.Link(star, "star3:diblock1:polyB.end2", "star2:diblock3:polyB.end2");
    world.Link(star, "star4:diblock1:polyB.end2", "star3:diblock3:polyB.end2");
    world.Link(star, "star5:diblock1:polyB.end2", "star4:diblock3:polyB.end2");

    world.Add(chain, "starchain");
    EXPECT_TRUE(world.hasName("starchain"));
    ASSERT_NO_THROW(world.FormFactor("starchain"));
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
        world.Link(new GaussianPolymer(), "arm" + std::to_string(i) + ".end1", "center.point");
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

// Error message diagnostics: bad reference point names should list valid ones
TEST(StructureTest, BadReferencePointListsValidOnes) {
    // ThinRod has: end1, end2, middle
    // Using an unknown reference point should throw and name all valid ones.
    {
        World w;
        w.Add(new ThinRod(), "rod1");
        try {
            w.Link(new ThinRod(), "rod2.badref", "rod1.end2");
            FAIL() << "Expected SEBException for unknown reference point on new subunit";
        } catch (const SEBException& e) {
            std::string msg = e.what();
            EXPECT_NE(msg.find("badref"), std::string::npos) << "Message should mention the bad name";
            EXPECT_NE(msg.find("end1"),   std::string::npos) << "Message should list 'end1'";
            EXPECT_NE(msg.find("end2"),   std::string::npos) << "Message should list 'end2'";
            EXPECT_NE(msg.find("middle"), std::string::npos) << "Message should list 'middle'";
        }
    }
    {
        // Same check for the *existing* subunit (second throw site)
        World w;
        w.Add(new ThinRod(), "rod1");
        try {
            w.Link(new ThinRod(), "rod2.end1", "rod1.badref");
            FAIL() << "Expected SEBException for unknown reference point on existing subunit";
        } catch (const SEBException& e) {
            std::string msg = e.what();
            EXPECT_NE(msg.find("badref"), std::string::npos) << "Message should mention the bad name";
            EXPECT_NE(msg.find("end1"),   std::string::npos) << "Message should list 'end1'";
            EXPECT_NE(msg.find("end2"),   std::string::npos) << "Message should list 'end2'";
            EXPECT_NE(msg.find("middle"), std::string::npos) << "Message should list 'middle'";
        }
    }
    {
        // Third throw site: testPathSyntax called via PhaseFactor on a structure path
        World w;
        GraphID g = w.Add(new ThinRod(), "rod1");
        w.Add(g, "chain");
        try {
            w.PhaseFactor("chain:rod1.badref", "chain:rod1.end2");
            FAIL() << "Expected SEBException for unknown reference point in path syntax check";
        } catch (const SEBException& e) {
            std::string msg = e.what();
            EXPECT_NE(msg.find("badref"), std::string::npos) << "Message should mention the bad name";
            EXPECT_NE(msg.find("end1"),   std::string::npos) << "Message should list 'end1'";
            EXPECT_NE(msg.find("end2"),   std::string::npos) << "Message should list 'end2'";
            EXPECT_NE(msg.find("middle"), std::string::npos) << "Message should list 'middle'";
        }
    }
}