#include <iostream>
#include "SEB.hpp"
#include "Symbolic.hpp"
#include "SymbolInterface.hpp"

int main() {
    try {
        sebsym::initialize();
        if (!sebsym::set_backend("ginac")) {
            throw std::runtime_error("This example requires the GiNaC symbolic backend");
        }

        SymbolInterface *symbols = SymbolInterface::instance();

        Expression alpha = symbols->getSymbol("alpha");
        Expression F = symbols->getSymbol("F", "chain");
        Expression A = symbols->getSymbol("A", "polymer", "end1");
        Expression P = symbols->getSymbol("Psi", "polymer", "end1", "end2");

        Expression test = alpha*F + 2*A*P*A;

        std::cout << "Backend: " << sebsym::active_backend() << std::endl;
        std::cout << "Symbolic expression: " << test << std::endl;
        std::cout << "Symbolic expression (LaTeX): " << test.to_latex() << std::endl;
        std::cout << "Symbolic expression (Python): " << test.to_python() << std::endl;
        std::cout << "Symbolic expression (C code): " << test.to_cform() << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
