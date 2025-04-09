void init_symbolic(py::module& m) {
    // Bind the SymbolicExpression class
    py::class_<SymbolicExpression, PySymbolicExpression, std::shared_ptr<SymbolicExpression>>(m, "SymbolicExpression")
        .def(py::init<>())
        .def("add", &SymbolicExpression::add)
        .def("subtract", &SymbolicExpression::subtract)
        .def("multiply", &SymbolicExpression::multiply)
        .def("divide", &SymbolicExpression::divide)
        .def("negate", &SymbolicExpression::negate)
        .def("pow", &SymbolicExpression::pow)
        .def("exp", &SymbolicExpression::exp)
        .def("log", &SymbolicExpression::log)
        .def("sin", &SymbolicExpression::sin)
        .def("cos", &SymbolicExpression::cos)
        .def("tan", &SymbolicExpression::tan)
        .def("asin", &SymbolicExpression::asin)
        .def("acos", &SymbolicExpression::acos)
        .def("atan", &SymbolicExpression::atan)
        .def("sinh", &SymbolicExpression::sinh)
        .def("cosh", &SymbolicExpression::cosh)
        .def("tanh", &SymbolicExpression::tanh)
        .def("sqrt", &SymbolicExpression::sqrt)
        .def("abs", &SymbolicExpression::abs)
        .def("bessel_j0", &SymbolicExpression::bessel_j0)
        .def("bessel_j1", &SymbolicExpression::bessel_j1)
        .def("dawson", &SymbolicExpression::dawson)
        .def("erf", &SymbolicExpression::erf)
        .def("erfc", &SymbolicExpression::erfc)
        .def("subs", py::overload_cast<const std::string&, const SymExprPtr&>(&SymbolicExpression::subs, py::const_))
        .def("subs", py::overload_cast<const ParameterMap&>(&SymbolicExpression::subs, py::const_))
        .def("eval", &SymbolicExpression::eval)
        .def("is_numeric", &SymbolicExpression::is_numeric)
        .def("to_string", &SymbolicExpression::to_string)
        .def("to_latex", &SymbolicExpression::to_latex)
        .def("to_python", &SymbolicExpression::to_python)
        .def("to_cform", &SymbolicExpression::to_cform)
        .def("__str__", &SymbolicExpression::to_string)
        .def("__repr__", &SymbolicExpression::to_string);

    // Bind the SymbolicFactory class
    py::class_<SymbolicFactory, PySymbolicFactory>(m, "SymbolicFactory")
        .def(py::init<>())
        .def("createSymbol", &SymbolicFactory::createSymbol)
        .def("createConstant", &SymbolicFactory::createConstant)
        .def("createPi", &SymbolicFactory::createPi)
        .def("createE", &SymbolicFactory::createE)
        .def("createI", &SymbolicFactory::createI)
        .def_static("instance", &SymbolicFactory::instance)
        .def_static("setInstance", &SymbolicFactory::setInstance);

#ifdef USE_GINAC
    // Bind the GiNaCExpression class
    py::class_<GiNaCExpression, SymbolicExpression, std::shared_ptr<GiNaCExpression>>(m, "GiNaCExpression")
        .def(py::init<const GiNaC::ex&>());

    // Bind the GiNaCFactory class
    py::class_<GiNaCFactory, SymbolicFactory>(m, "GiNaCFactory")
        .def(py::init<>());
#endif

    // Bind helper functions
    m.def("symbol", &symbol);
    m.def("constant", &constant);
    m.def("pi", &pi);
    m.def("e", &e);
    m.def("i", &i);
}
