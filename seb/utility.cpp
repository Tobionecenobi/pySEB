#include "World.hpp"
#include <string>
#include <vector>

// Utility function to create a random walk polymer
void createRandomWalkPolymer(World& world, int numRods, const std::string& baseName) {
    // Define the first rod in the chain
    GraphID firstRod = world.Add(new ThinRod(), baseName + "1", "rod");

    // Add and link the remaining rods
    for (int i = 2; i <= numRods; ++i) {
        std::string currentRod = baseName + std::to_string(i);
        std::string previousRod = baseName + std::to_string(i - 1);
        world.Link(new ThinRod(), currentRod + ".end1", previousRod + ".end2", "rod");
    }

    // Name the structure formed by the rods
    world.Add(firstRod, baseName + "Structure");
}

// Utility function to evaluate form factors and save results
void evaluateFormFactor(World& world, const std::string& structureName, const std::string& outputFileName) {
    // Obtain the form factor expression
    ex formFactor = world.FormFactor(structureName, WORLDMAXDEPTH, QVAR);

    // Set parameters
    ParameterList params;
    world.setParameter(params, "beta_rod", 1);
    world.setParameter(params, "L_rod", sqrt(12));

    // Generate q-values
    DoubleVector qValues = world.logspace(0.01, 10.0, 200);

    // Evaluate and save the form factor
    world.Evaluate(formFactor, params, qValues, outputFileName, "Form factor evaluation");
}