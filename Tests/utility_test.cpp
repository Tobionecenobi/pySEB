#include "gtest/gtest.h"
#include "World.hpp"
#include "utility.cpp"
#include "GiNaCSymbolic.hpp"

// Global initialization of the SymbolicFactory
struct SymbolicFactoryInitializer {
    SymbolicFactoryInitializer() {
        GiNaCSymbolic* factory = new GiNaCSymbolic();
        SymbolicFactory::setInstance(factory);
    }
} g_initializer;

// Test for creating a random walk polymer
TEST(UtilityTest, CreateRandomWalkPolymer) {
    World world;
    createRandomWalkPolymer(world, 10, "TestRod");

    // Verify that the structure was created
    EXPECT_TRUE(world.hasName("TestRod1"));
    EXPECT_TRUE(world.hasName("TestRodStructure"));
}

// Test for evaluating form factors
TEST(UtilityTest, EvaluateFormFactor) {
    World world;
    createRandomWalkPolymer(world, 10, "TestRod");

    // Evaluate form factor and save to a file
    ASSERT_NO_THROW(evaluateFormFactor(world, "TestRodStructure", "test_output.q"));

    // Verify that the output file is created (mock or check file existence in real scenarios)
    // For now, we assume the function runs without exceptions
}