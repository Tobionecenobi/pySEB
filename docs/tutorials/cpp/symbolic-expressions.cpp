#include <iostream>

#include "SEB.hpp"
#include "Symbolic.hpp"

int main()
{
    sebsym::initialize();
    sebsym::set_backend("portable");

    World world("SymbolicTutorial");
    world.Add("GaussianPolymer", "polymer");

    Expression formFactor = world.FormFactor("polymer");

    std::cout << "backend=" << sebsym::active_backend() << '\n';
    std::cout << "default=" << formFactor << '\n';
    std::cout << "python=" << formFactor.to_python() << '\n';
    std::cout << "latex=" << formFactor.to_latex() << '\n';
    std::cout << "c=" << formFactor.to_cform() << '\n';

    ParameterList parameters;
    world.setParameter(parameters, "beta_polymer", 1.0);
    world.setParameter(parameters, "Rg_polymer", 5.0);
    std::cout << "F(q=0.1)="
              << world.Evaluate(formFactor, parameters, 0.1)
              << '\n';
}
