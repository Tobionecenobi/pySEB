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
