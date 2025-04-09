#include "SymbolicInterface.hpp"
#include <stdexcept>

// Initialize the static instance pointer
SymbolicFactory* SymbolicFactory::_instance = nullptr;

// Get the singleton instance
SymbolicFactory* SymbolicFactory::instance() {
    if (_instance == nullptr) {
        throw std::runtime_error("SymbolicFactory instance not set");
    }
    return _instance;
}

// Set the factory implementation
void SymbolicFactory::setInstance(SymbolicFactory* factory) {
    _instance = factory;
}
