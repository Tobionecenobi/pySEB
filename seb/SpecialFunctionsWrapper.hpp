#ifndef INCLUDE_SPECIALFUNCTIONS_WRAPPER_HPP
#define INCLUDE_SPECIALFUNCTIONS_WRAPPER_HPP

#include "Expression.hpp"
#include "SymbolicInterface.hpp"

#ifdef USE_GINAC
#include "GiNaCSymbolic.hpp"

// Only define if not already defined by GiNaC
#ifndef REGISTER_FUNCTION
#define REGISTER_FUNCTION(name, options) \
    namespace { \
        struct name##_register_helper { \
            name##_register_helper() { \
                std::cout << "Registering function " << #name << std::endl; \
            } \
        } name##_register_helper_instance; \
    }
#endif

#endif // USE_GINAC

// Function declarations that work with our Expression type
Expression BesselJ0_wrapper(const Expression& x);
Expression BesselJ1_wrapper(const Expression& x);
Expression DawsonF_wrapper(const Expression& x);
Expression Erf_wrapper(const Expression& x);
Expression Erfc_wrapper(const Expression& x);

#endif // INCLUDE_SPECIALFUNCTIONS_WRAPPER_HPP
