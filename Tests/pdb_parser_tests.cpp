#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "SEB.hpp"

namespace {

const char* pdbModel =
    "MODEL        1\n"
    "ATOM      1  N   ALA A   1      11.104  13.207   9.234  0.50 20.00           N  \n"
    "ATOM      2  CA AALA A   1      12.560  13.100   9.500  1.00 21.00           C  \n"
    "ATOM      3  CA BALA A   1      12.700  13.200   9.600  1.00 22.00           C  \n"
    "HETATM    4  O   HOH A 101      15.000  16.000  17.000  1.00 30.00           O  \n"
    "ENDMDL\n"
    "MODEL        2\n"
    "ATOM      5  N   GLY B   2      21.000  22.000  23.000  1.00 40.00           N  \n"
    "ENDMDL\n";

}

// Verify fixed-column parsing and the conservative defaults: model 1,
// ATOM records only, primary alternate locations, and no HETATM water.
TEST(PDBParserTest, ParsesPrimaryAtomsFromFirstModel)
{
    std::istringstream input(pdbModel);
    const PDBParser parser;
    const std::vector<AtomRecord> atoms =
        parser.Parse(input, StructureParseOptions(), "example.pdb");

    ASSERT_EQ(atoms.size(), 2u);
    EXPECT_EQ(atoms[0].serial, 1);
    EXPECT_EQ(atoms[0].atomName, "N");
    EXPECT_EQ(atoms[0].residueName, "ALA");
    EXPECT_EQ(atoms[0].chainId, "A");
    EXPECT_EQ(atoms[0].residueNumber, 1);
    EXPECT_EQ(atoms[0].element, "N");
    EXPECT_DOUBLE_EQ(atoms[0].occupancy, 0.5);
    EXPECT_NEAR(atoms[0].coordinate.x, 11.104, 1e-12);

    EXPECT_EQ(atoms[1].serial, 2);
    EXPECT_EQ(atoms[1].alternateLocation, "A");
    EXPECT_EQ(atoms[1].element, "C");
}

// Verify that parser options select another MODEL and may include HETATM
// records, water, and every alternate conformation.
TEST(PDBParserTest, HonorsModelHetatmWaterAndAlternateOptions)
{
    const PDBParser parser;

    StructureParseOptions allFirstModel;
    allFirstModel.includeHetatm = true;
    allFirstModel.includeWater = true;
    allFirstModel.alternateLocations = AlternateLocationPolicy::All;
    std::istringstream firstInput(pdbModel);
    const std::vector<AtomRecord> firstModel =
        parser.Parse(firstInput, allFirstModel, "example.pdb");
    ASSERT_EQ(firstModel.size(), 4u);
    EXPECT_EQ(firstModel.back().recordType, "HETATM");
    EXPECT_EQ(firstModel.back().residueName, "HOH");

    StructureParseOptions secondModel;
    secondModel.modelNumber = 2;
    std::istringstream secondInput(pdbModel);
    const std::vector<AtomRecord> selected =
        parser.Parse(secondInput, secondModel, "example.pdb");
    ASSERT_EQ(selected.size(), 1u);
    EXPECT_EQ(selected[0].serial, 5);
    EXPECT_EQ(selected[0].modelNumber, 2);
    EXPECT_EQ(selected[0].chainId, "B");
}

// Verify profile precedence and cloud construction. A residue/atom-specific
// rule must override the element rule, occupancy optionally scales beta, and
// selected atom serials become named cloud reference points.
TEST(AtomCloudBuilderTest, AppliesProfileOccupancyAndReferences)
{
    std::istringstream input(pdbModel);
    const PDBParser parser;
    const std::vector<AtomRecord> atoms =
        parser.Parse(input, StructureParseOptions(), "example.pdb");

    AtomParameterProfile profile;
    profile.SetElement("N", 1.0, 2.0);
    profile.SetElement("C", 1.5, -1.0);
    profile.SetAtom("ALA", "CA", 1.7, 3.0);

    AtomCloudBuildOptions options;
    options.referenceAtomSerials["n_term"] = 1;
    std::unique_ptr<DebyeSphereCloud> cloud =
        AtomCloudBuilder::Build(atoms, profile, options);

    ASSERT_EQ(cloud->scattererCount(), 2u);
    EXPECT_NEAR(cloud->NumericTotalBeta(ParameterList()), 4.0, 1e-12);

    const std::map<refPoint, CartesianPoint3D> references =
        cloud->getReferenceCoordinates();
    ASSERT_NE(references.find("n_term"), references.end());
    EXPECT_NEAR(references.at("n_term").x, 11.104, 1e-12);

    World world;
    world.Add(cloud.release(), "protein");
    EXPECT_NEAR(
        world.EvaluateFormFactor("protein", ParameterList(), 0.0),
        1.0,
        1e-12
    );
}

// Verify that malformed fixed-column numeric fields produce an actionable
// parser exception containing the source filename and line number.
TEST(PDBParserTest, ReportsSourceAndLineForMalformedCoordinates)
{
    const std::string malformed =
        "ATOM      1  N   ALA A   1      INVALID 13.207   9.234  1.00 20.00           N  \n";
    std::istringstream input(malformed);
    const PDBParser parser;

    try {
        parser.Parse(input, StructureParseOptions(), "broken.pdb");
        FAIL() << "Expected malformed PDB coordinate to throw";
    } catch (const SEBException& exception) {
        EXPECT_NE(
            std::string(exception.what()).find("broken.pdb:1"),
            std::string::npos
        );
    }
}
