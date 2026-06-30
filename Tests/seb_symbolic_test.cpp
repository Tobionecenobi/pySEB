#include "Symbolic.hpp"
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

    EXPECT_THROW(expr.series(x, sebsym::constant(0.0), 3), std::runtime_error);
}
