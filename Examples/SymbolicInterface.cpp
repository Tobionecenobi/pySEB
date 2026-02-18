#include <iostream>
#include "SEB.hpp"
#include "GiNaCSymbolic.hpp"

int main() {
    try {
        // Initialize the SymbolicFactory with GiNaC
        SymbolicFactory::setInstance(new GiNaCSymbolic());
        std::cout << "Using GiNaC for symbolic computations" << std::endl;
        
        // Create a world
        World w("GiNaCWorld");
        
        // Add a diblock copolymer
        w.Add("diblock", "diblock");
        
        // Get the form factor expression
        Expression F = w.FormFactor("diblock");
        
        // Print in different formats
        std::cout << "\nForm Factor (default):" << std::endl;
        std::cout << F << std::endl;
        
        std::cout << "\nForm Factor (LaTeX):" << std::endl;
        std::cout << F.to_latex() << std::endl;
        
        std::cout << "\nForm Factor (C code):" << std::endl;
        std::cout << F.to_cform() << std::endl;
        
        std::cout << "\nForm Factor (Python):" << std::endl;
        std::cout << F.to_python() << std::endl;
        
        // Evaluate at q=0.1
        // This would require substituting q with 0.1 and evaluating
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
