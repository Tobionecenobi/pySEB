#include "SymbolicInterface.hpp"
#include "SymbolicPortable.hpp"
#include <stdexcept>

#ifdef USE_GINAC_IMPL
namespace sebsym {
void registerGiNaCBackend();
}
#endif

// Initialize the static instance pointer
SymbolicFactory* SymbolicFactory::_instance = nullptr;
std::string SymbolicFactory::_activeBackend;
std::map<std::string, SymbolicFactory*> SymbolicFactory::_backends;

namespace {
std::map<std::string, SymbolicFactory*>& backend_registry()
{
    static std::map<std::string, SymbolicFactory*> backends;
    return backends;
}

std::string& active_backend_name()
{
    static std::string name;
    return name;
}
}

// Get the singleton instance
SymbolicFactory* SymbolicFactory::instance() {
    if (_instance == nullptr) {
        sebsym::registerDefaultBackends();
        setBackend("portable");
    }
    if (_instance == nullptr) {
        throw std::runtime_error("SymbolicFactory instance not set");
    }
    return _instance;
}

// Set the factory implementation
void SymbolicFactory::setInstance(SymbolicFactory* factory) {
    _instance = factory;
    if (factory) {
        active_backend_name() = factory->backendName();
        backend_registry()[active_backend_name()] = factory;
    }
}

void SymbolicFactory::registerBackend(const std::string& name, SymbolicFactory* factory) {
    if (!factory) {
        throw std::runtime_error("Cannot register null symbolic backend: " + name);
    }
    backend_registry()[name] = factory;
}

bool SymbolicFactory::setBackend(const std::string& name) {
    auto& backends = backend_registry();
    auto it = backends.find(name);
    if (it == backends.end()) return false;
    _instance = it->second;
    active_backend_name() = name;
    return true;
}

std::vector<std::string> SymbolicFactory::availableBackends() {
    sebsym::registerDefaultBackends();
    std::vector<std::string> names;
    for (const auto& backend : backend_registry()) {
        names.push_back(backend.first);
    }
    return names;
}

std::string SymbolicFactory::activeBackendName() {
    if (_instance == nullptr) {
        instance();
    }
    return active_backend_name();
}

SymbolicCapabilities SymbolicFactory::activeCapabilities() {
    return instance()->capabilities();
}

namespace sebsym {

void registerDefaultBackends()
{
    registerPortableBackend();
#ifdef USE_GINAC_IMPL
    registerGiNaCBackend();
#endif
}

}
