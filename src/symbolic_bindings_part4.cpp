// String representation
std::string PySymbolicExpression::to_string() const override {
    PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_string, );
}

std::string PySymbolicExpression::to_latex() const override {
    PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_latex, );
}

std::string PySymbolicExpression::to_python() const override {
    PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_python, );
}

std::string PySymbolicExpression::to_cform() const override {
    PYBIND11_OVERRIDE_PURE(std::string, SymbolicExpression, to_cform, );
}

// Forward declaration of the Python symbolic factory class
class PySymbolicFactory : public SymbolicFactory {
public:
    using SymbolicFactory::SymbolicFactory;

    SymExprPtr createSymbol(const std::string& name) override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createSymbol, name);
    }

    SymExprPtr createConstant(double value) override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createConstant, value);
    }

    SymExprPtr createPi() override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createPi, );
    }

    SymExprPtr createE() override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createE, );
    }

    SymExprPtr createI() override {
        PYBIND11_OVERRIDE_PURE(SymExprPtr, SymbolicFactory, createI, );
    }
};
