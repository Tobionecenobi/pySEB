#include "GiNaCSymbolic.hpp"

#ifdef USE_GINAC_IMPL

// Implementation of GiNaCExpression methods that couldn't be inlined
// (Most methods are already implemented in the header file)

// Helper function to create a GiNaCExpression from a GiNaC expression
SymExprPtr make_ginac_expr(const GiNaC::ex& expr) {
    return std::make_shared<GiNaCExpression>(expr);
}

// Implementation of GiNaCFactory methods
SymExprPtr GiNaCFactory::createSymbol(const std::string& name) {
    return std::make_shared<GiNaCExpression>(GiNaC::symbol(name));
}

SymExprPtr GiNaCFactory::createConstant(double value) {
    return std::make_shared<GiNaCExpression>(GiNaC::numeric(value));
}

SymExprPtr GiNaCFactory::createPi() {
    return std::make_shared<GiNaCExpression>(GiNaC::Pi);
}

SymExprPtr GiNaCFactory::createE() {
    return std::make_shared<GiNaCExpression>(GiNaC::exp(1));
}

SymExprPtr GiNaCFactory::createI() {
    return std::make_shared<GiNaCExpression>(GiNaC::I);
}

// Create a global instance of the GiNaC factory
static GiNaCFactory ginacFactory;

// Initialize the factory at program startup
struct GiNaCSymbolicInitializer {
    GiNaCSymbolicInitializer() {
        SymbolicFactory::registerBackend("ginac", &ginacFactory);
        SymbolicFactory::setBackend("ginac");
    }
};

// Create a static instance of the initializer
static GiNaCSymbolicInitializer initializer;

#endif // USE_GINAC_IMPL
