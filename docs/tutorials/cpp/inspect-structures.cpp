#include <iostream>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    if (!sebsym::set_backend("portable")) {
        std::cerr << "portable backend unavailable\n";
        return 1;
    }

    World world("StructureInspectionTutorial");

    GraphID diblock = world.Add("GaussianPolymer", "A");
    world.Link("GaussianPolymer", "B.end1", "A.end2");

    GraphID pair = world.Add(diblock, "left");
    world.Link(diblock, "right:A.end1", "left:A.end1");
    world.Add(pair, "pair");

    std::cout << std::boolalpha;
    std::cout << "has pair: " << world.hasName("pair") << '\n';
    std::cout << "pair is structure: " << world.isStructure("pair") << '\n';
    std::cout << "A is subunit: " << world.isSubunit("A") << '\n';
    std::cout << "graph containing pair: " << world.getGraphID("pair") << "\n\n";
    std::cout << "nested instances linked: "
              << world.isLinked("left:A.end1", "right:A.end1") << "\n\n";

    std::cout << "Links\n";
    world.printLinks();

    std::cout << "\nGraph IDs\n";
    world.printGraphIDs();

    std::cout << "\nStructure tree\n";
    world.folderprint("pair");
}
