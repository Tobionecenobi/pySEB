#include "Symbolic.hpp"
#include "Types.hpp"
#ifdef USE_GINAC_IMPL
#include "GiNaCSymbolic.hpp"
#endif

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

void register_optional_ginac_backend()
{
#ifdef USE_GINAC_IMPL
    static GiNaCFactory ginac_factory;
    SymbolicFactory::registerBackend("ginac", &ginac_factory);
#endif
}

struct SebSymbolicFixture : public ::testing::Test {
    void SetUp() override
    {
        sebsym::initialize();
        register_optional_ginac_backend();
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

double evaluate_parity_expression(const std::string& backend)
{
    if (!sebsym::set_backend(backend)) {
        throw std::runtime_error("backend is not available: " + backend);
    }

    return parity_expression()
        .subs("x", 1.25)
        .subs("y", 2.5)
        .eval();
}

} // namespace

TEST_F(SebSymbolicFixture, PortableBackendBuildsSubstitutesAndEvaluates)
{
    auto x = sebsym::symbol("x");
    auto expr = (x * sebsym::constant(2.0) + sebsym::constant(3.0)).sin();
    auto numeric = expr.subs("x", 4.0);

    EXPECT_EQ(sebsym::active_backend(), "portable");
    EXPECT_NE(numeric.to_python().find("sin"), std::string::npos);
    EXPECT_NEAR(numeric.eval(), std::sin(11.0), 1e-12);
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

TEST_F(SebSymbolicFixture, PortableAndGiNaCAgreeOnCommonNumericSubset)
{
    auto backends = sebsym::available_backends();
    if (std::find(backends.begin(), backends.end(), "ginac") == backends.end()) {
        GTEST_SKIP() << "GiNaC backend is not available";
    }

    double portable_value = evaluate_parity_expression("portable");
    double ginac_value = evaluate_parity_expression("ginac");

    EXPECT_NEAR(portable_value, ginac_value, 1e-12);
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
    auto x = sebsym::symbol("x");
    std::string python = (x.erf() + x.erfc() + x.bessel_j0() + x.bessel_j1()).to_python();

    EXPECT_NE(python.find("erf"), std::string::npos);
    EXPECT_NE(python.find("erfc"), std::string::npos);
    EXPECT_NE(python.find("besselj(0"), std::string::npos);
    EXPECT_NE(python.find("besselj(1"), std::string::npos);
}
