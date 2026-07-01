#include <iostream>
#include <memory>
#include "Expression.hpp"
#include "GiNaCSymbolic.hpp"

// Test basic operations of the Expression class
void test_basic_operations() {
    std::cout << "Testing basic operations..." << std::endl;
    
    // Create a symbolic factory
    auto factory = std::make_shared<GiNaCSymbolic>();
    
    // Create symbols
    Expression x = factory->createSymbol("x");
    Expression y = factory->createSymbol("y");
    Expression z = factory->createSymbol("z");
    
    // Test basic arithmetic operations
    Expression sum = x + y;
    Expression diff = x - y;
    Expression prod = x * y;
    Expression quot = x / y;
    Expression neg = -x;
    
    std::cout << "x + y = " << sum.to_string() << std::endl;
    std::cout << "x - y = " << diff.to_string() << std::endl;
    std::cout << "x * y = " << prod.to_string() << std::endl;
    std::cout << "x / y = " << quot.to_string() << std::endl;
    std::cout << "-x = " << neg.to_string() << std::endl;
    
    // Test functions
    Expression exp_x = x.exp();
    Expression log_x = x.log();
    Expression sin_x = x.sin();
    Expression cos_x = x.cos();
    
    std::cout << "exp(x) = " << exp_x.to_string() << std::endl;
    std::cout << "log(x) = " << log_x.to_string() << std::endl;
    std::cout << "sin(x) = " << sin_x.to_string() << std::endl;
    std::cout << "cos(x) = " << cos_x.to_string() << std::endl;
    
    // Test substitution
    Expression expr = x * x + y;
    Expression subst = expr.subs("x", constant(2));
    
    std::cout << "x^2 + y = " << expr.to_string() << std::endl;
    std::cout << "After substituting x=2: " << subst.to_string() << std::endl;
    
    // Test evaluation
    Expression num_expr = constant(3) * constant(4);
    double result = num_expr.eval();
    
    std::cout << "3 * 4 = " << result << std::endl;
    
    std::cout << "Basic operations test completed." << std::endl;
}

// Test series expansion and coefficient extraction
void test_series_expansion() {
    std::cout << "\nTesting series expansion and coefficient extraction..." << std::endl;
    
    // Create a symbolic factory
    auto factory = std::make_shared<GiNaCSymbolic>();
    
    // Create symbols
    Expression x = factory->createSymbol("x");
    Expression y = factory->createSymbol("y");
    
    // Test series expansion
    Expression expr = x.sin();
    Expression series_expr = expr.series(x, constant(0), 5);
    
    std::cout << "sin(x) = " << expr.to_string() << std::endl;
    std::cout << "Series expansion of sin(x) around x=0 up to order 5: " << series_expr.to_string() << std::endl;
    
    // Test series to polynomial conversion
    Expression poly = series_expr.series_to_poly(series_expr);
    std::cout << "Series as polynomial: " << poly.to_string() << std::endl;
    
    // Test coefficient extraction
    Expression coeff1 = poly.coeff(x, 1);
    Expression coeff3 = poly.coeff(x, 3);
    
    std::cout << "Coefficient of x^1: " << coeff1.to_string() << std::endl;
    std::cout << "Coefficient of x^3: " << coeff3.to_string() << std::endl;
    
    std::cout << "Series expansion test completed." << std::endl;
}

// Test numeric evaluation and conversion
void test_numeric_evaluation() {
    std::cout << "\nTesting numeric evaluation and conversion..." << std::endl;
    
    // Create a symbolic factory
    auto factory = std::make_shared<GiNaCSymbolic>();
    
    // Create constants
    Expression pi_val = factory->createPi();
    Expression e_val = factory->createE();
    
    std::cout << "pi = " << pi_val.to_string() << std::endl;
    std::cout << "e = " << e_val.to_string() << std::endl;
    
    // Test numeric evaluation
    Expression expr = pi_val / constant(4);
    Expression eval_expr = expr.evalf();
    
    std::cout << "pi/4 = " << expr.to_string() << std::endl;
    std::cout << "Evaluated: " << eval_expr.to_string() << std::endl;
    
    // Test to_double conversion
    double pi_4 = expr.to_double();
    std::cout << "pi/4 as double: " << pi_4 << std::endl;
    
    // Test is_zero and is_numeric
    Expression zero = constant(0);
    Expression one = constant(1);
    Expression sym = factory->createSymbol("z");
    
    std::cout << "Is 0 zero? " << (zero.is_zero() ? "Yes" : "No") << std::endl;
    std::cout << "Is 1 zero? " << (one.is_zero() ? "Yes" : "No") << std::endl;
    std::cout << "Is 0 numeric? " << (zero.is_numeric() ? "Yes" : "No") << std::endl;
    std::cout << "Is z numeric? " << (sym.is_numeric() ? "Yes" : "No") << std::endl;
    
    std::cout << "Numeric evaluation test completed." << std::endl;
}

int main() {
    std::cout << "Expression Class Test Program" << std::endl;
    std::cout << "============================" << std::endl;
    
    try {
        test_basic_operations();
        test_series_expansion();
        test_numeric_evaluation();
        
        std::cout << "\nAll tests completed successfully!" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
