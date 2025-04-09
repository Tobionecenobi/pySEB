// Standard C++ headers
#include<iostream>

// Include SEB functionality
#include "SEB.hpp"
#include "GiNaCSymbolic.hpp"
#include <memory>

// Global initialization of the SymbolicFactory
struct SymbolicFactoryInitializer {
    SymbolicFactoryInitializer() {
        std::cout << "Setting SymbolicFactory instance..." << std::endl;
        GiNaCSymbolic* factory = new GiNaCSymbolic();
        std::cout << "Created factory: " << (factory != nullptr ? "OK" : "NULL") << std::endl;
        SymbolicFactory::setInstance(factory);
        std::cout << "SymbolicFactory instance set." << std::endl;
    }
} g_initializer;

/*

    In this example we generate a tri-block copolymer and derive various expressions

                A                    B                     C
        x----------------x = x-----------------x = x-----------------x
         end1        end2     end1         end2     end1         end2

*/

int main()
{
 try{
    World w("World");

    // Define diblock copolymer
    GraphID diblock = w.Add(new GaussianPolymer(), "A");
    w.Link(new GaussianPolymer(),  "B.end1", "A.end2");

    // Name diblock structure
    w.Add(diblock, "diblock");


    Expression F = w.FormFactor("diblock");

    cout << "\nForm Factor (default output) --------------------------------------\n";
    cout << F << endl;

    cout << "\nForm Factor (latex formatted) -------------------------------------------\n";
    cout << F.to_latex() << endl;

    cout << "\nForm Factor (c formatted) -----------------------------------------------\n";
    cout << F.to_cform() << endl;

    cout << "\nForm Factor (Python) ----------------------------------------------------\n";
    cout << F.to_python() << endl;

}
catch (const SEBException e)
{
    std::cout << e;                    // Print what the error was, and where it was triggered.
}



}
