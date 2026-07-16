#include "Symbolic.hpp"
#include "ExpressionFunctions.hpp"
#include "SpecialFunctions.hpp"
#include "SymbolInterface.hpp"
#include "Types.hpp"

#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <type_traits>

static_assert(std::is_same<Expression, sebsym::Expression>::value,
              "global Expression should remain a compatibility alias");
static_assert(std::is_same<ex, sebsym::Expression>::value,
              "SEB expression typedef should use sebsym::Expression");
static_assert(std::is_same<ExpressionMap, std::map<sebsym::Expression, sebsym::Expression>>::value,
              "SEB expression maps should use sebsym::Expression");

namespace {

struct SebSymbolicFixture : public ::testing::Test {
    void SetUp() override
    {
        sebsym::initialize();
        ASSERT_TRUE(sebsym::set_backend("portable"));
    }
};

sebsym::Expression parity_expression()
{
    auto x = sebsym::symbol("x");
    auto y = sebsym::symbol("y");

    return (x * sebsym::constant(2.0) + y).sin()
        + (x / sebsym::constant(3.0)).exp()
        + (y + sebsym::constant(4.0)).log()
        + (x + sebsym::e()).sqrt()
        + sebsym::pi().sin();
}

double evaluate_parity_expression_with_ginac()
{
    if (!sebsym::set_backend("ginac")) {
        throw std::runtime_error("backend is not available: ginac");
    }

    return parity_expression()
        .subs("x", 1.25)
        .subs("y", 2.5)
        .eval();
}

sebsym::Expression special_function_expression()
{
    auto x = sebsym::symbol("x");
    return (x / sebsym::constant(2.0)).erf()
        + (x / sebsym::constant(5.0)).erfc()
        + x.bessel_j0()
        + x.bessel_j1()
        + (x / sebsym::constant(3.0)).dawson();
}

double evaluate_special_function_expression_with_ginac()
{
    if (!sebsym::set_backend("ginac")) {
        throw std::runtime_error("backend is not available: ginac");
    }

    return special_function_expression()
        .subs("x", 1.25)
        .eval();
}

} // namespace

TEST_F(SebSymbolicFixture, PortableBackendBuildsAndExports)
{
    auto x = sebsym::symbol("x");
    auto expr = (x * sebsym::constant(2.0) + sebsym::constant(3.0)).sin();

    EXPECT_EQ(sebsym::active_backend(), "portable");
    EXPECT_NE(expr.to_python().find("sin"), std::string::npos);
    EXPECT_NEAR(expr.subs("x", 4.0).eval(), std::sin(11.0), 1e-12);
    EXPECT_THROW(expr.eval(), std::runtime_error);
}

TEST_F(SebSymbolicFixture, RegistryReportsPortableCapabilities)
{
    auto backends = sebsym::available_backends();
    EXPECT_NE(std::find(backends.begin(), backends.end(), "portable"), backends.end());

    auto caps = sebsym::active_capabilities();
    EXPECT_TRUE(caps.numeric_evaluation);
    EXPECT_TRUE(caps.python_output);
    EXPECT_FALSE(caps.series_expansion);
}

TEST_F(SebSymbolicFixture, PortableBackendRejectsUnsupportedSeries)
{
    auto x = sebsym::symbol("x");
    auto expr = x.sin();

    try {
        expr.series(x, sebsym::constant(0.0), 3);
        FAIL() << "portable backend should reject series expansion";
    } catch (const std::runtime_error& error) {
        std::string message = error.what();
        EXPECT_NE(message.find("portable"), std::string::npos);
        EXPECT_NE(message.find("series expansion"), std::string::npos);
    }

    EXPECT_THROW(expr.coeff(x, 1), std::runtime_error);
}

TEST_F(SebSymbolicFixture, GiNaCBackendReportsSeriesCapabilityWhenAvailable)
{
    auto backends = sebsym::available_backends();
    if (std::find(backends.begin(), backends.end(), "ginac") == backends.end()) {
        GTEST_SKIP() << "GiNaC backend is not available";
    }

    ASSERT_TRUE(sebsym::set_backend("ginac"));
    auto caps = sebsym::active_capabilities();

    EXPECT_TRUE(caps.numeric_evaluation);
    EXPECT_TRUE(caps.series_expansion);
    EXPECT_TRUE(caps.symbolic_simplification);
}

TEST_F(SebSymbolicFixture, GiNaCEvaluatesCommonNumericSubset)
{
    auto backends = sebsym::available_backends();
    if (std::find(backends.begin(), backends.end(), "ginac") == backends.end()) {
        GTEST_SKIP() << "GiNaC backend is not available";
    }

    double ginac_value = evaluate_parity_expression_with_ginac();

    EXPECT_NEAR(ginac_value, std::sin(5.0) + std::exp(1.25 / 3.0) +
                             std::log(6.5) + std::sqrt(1.25 + std::exp(1.0)),
                1e-12);
}

TEST_F(SebSymbolicFixture, PortableExportsCommonSubsetToSymPySyntax)
{
    ASSERT_TRUE(sebsym::set_backend("portable"));
    std::string python = parity_expression().to_python();

    EXPECT_NE(python.find("sin"), std::string::npos);
    EXPECT_NE(python.find("exp"), std::string::npos);
    EXPECT_NE(python.find("log"), std::string::npos);
    EXPECT_NE(python.find("sqrt"), std::string::npos);
}

TEST_F(SebSymbolicFixture, PortableExportsSpecialFunctionsToSymPySyntax)
{
    ASSERT_TRUE(sebsym::set_backend("portable"));
    std::string python = special_function_expression().to_python();

    EXPECT_NE(python.find("erf"), std::string::npos);
    EXPECT_NE(python.find("erfc"), std::string::npos);
    EXPECT_NE(python.find("besselj(0"), std::string::npos);
    EXPECT_NE(python.find("besselj(1"), std::string::npos);
    EXPECT_NE(python.find("dawson"), std::string::npos);
}

TEST_F(SebSymbolicFixture, PortableExportsSineIntegralInsteadOfSinOverX)
{
    ASSERT_TRUE(sebsym::set_backend("portable"));
    auto x = sebsym::symbol("x");

    std::string python = Six(x).to_python();

    EXPECT_NE(python.find("Si"), std::string::npos);
    EXPECT_EQ(python.find("sin"), std::string::npos);
}

TEST_F(SebSymbolicFixture, PortableExportsUnevaluatedIntegralToSymPySyntax)
{
    ASSERT_TRUE(sebsym::set_backend("portable"));
    auto t = sebsym::symbol("t");
    auto expr = integrate(t, 0, 1, t * t);

    EXPECT_NE(expr.to_python().find("Integral"), std::string::npos);
    EXPECT_THROW(expr.eval(), std::runtime_error);
}

TEST_F(SebSymbolicFixture, PortableExportsHypergeometricAndStruveToSymPySyntax)
{
    ASSERT_TRUE(sebsym::set_backend("portable"));
    auto x = sebsym::symbol("x");

    EXPECT_NE(Hypergeometric0F1Regularized(sebsym::constant(2.0), x).to_python().find("hyper"),
              std::string::npos);
    EXPECT_NE(StruveH0(x).to_python().find("struve(0"), std::string::npos);
    EXPECT_NE(StruveH1(x).to_python().find("struve(1"), std::string::npos);
}

TEST_F(SebSymbolicFixture, GiNaCEvaluatesSpecialFunctions)
{
    auto backends = sebsym::available_backends();
    if (std::find(backends.begin(), backends.end(), "ginac") == backends.end()) {
        GTEST_SKIP() << "GiNaC backend is not available";
    }

    double ginac_value = evaluate_special_function_expression_with_ginac();

    EXPECT_TRUE(std::isfinite(ginac_value));
}

TEST_F(SebSymbolicFixture, SymbolInterfaceBuildsExampleExpressionWithGiNaCBackend)
{
    auto backends = sebsym::available_backends();
    if (std::find(backends.begin(), backends.end(), "ginac") == backends.end()) {
        GTEST_SKIP() << "GiNaC backend is not available";
    }

    ASSERT_TRUE(sebsym::set_backend("ginac"));
    SymbolInterface* symbols = SymbolInterface::instance();

    sebsym::Expression alpha = symbols->getSymbol("alpha");
    sebsym::Expression alpha_again = symbols->getSymbol("alpha");
    sebsym::Expression form_factor = symbols->getSymbol("F", "chain");
    sebsym::Expression amplitude = symbols->getSymbol("A", "polymer", "end1");
    sebsym::Expression phase = symbols->getSymbol("Psi", "polymer", "end1", "end2");

    sebsym::Expression expression = alpha * form_factor + 2 * amplitude * phase * amplitude;

    EXPECT_EQ(alpha.to_string(), alpha_again.to_string());
    EXPECT_NE(expression.to_string().find("alpha"), std::string::npos);
    EXPECT_NE(expression.to_string().find("F_chain"), std::string::npos);
    EXPECT_NE(expression.to_string().find("A_polymer:end1"), std::string::npos);
    EXPECT_NE(expression.to_string().find("Psi_polymer:end1,end2"), std::string::npos);
    EXPECT_FALSE(expression.to_latex().empty());
    EXPECT_FALSE(expression.to_python().empty());
    EXPECT_FALSE(expression.to_cform().empty());
}
