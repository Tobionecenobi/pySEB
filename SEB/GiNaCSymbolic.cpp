#include "GiNaCSymbolic.hpp"

// Implementation of GiNaCExpression methods that couldn't be inlined
// (Most methods are already implemented in the header file)

// Helper function to create a GiNaCExpression from a GiNaC expression
SymExprPtr make_ginac_expr(const GiNaC::ex& expr) {
    return std::make_shared<GiNaCExpression>(expr);
}

// Create a global instance of the GiNaC factory
static GiNaCFactory ginacFactory;

// Initialize the factory at program startup
struct GiNaCFactoryInitializer {
    GiNaCFactoryInitializer() {
        SymbolicFactory::setInstance(&ginacFactory);
    }
};

// Create a static instance of the initializer
static GiNaCFactoryInitializer initializer;
