#ifndef INCLUDE_GUARD_SEB_SYMBOLIC
#define INCLUDE_GUARD_SEB_SYMBOLIC

#include "Expression.hpp"
#include "SymbolicInterface.hpp"
#include "SymbolicPortable.hpp"

namespace sebsym {

using Backend = ::SymbolicFactory;
using BackendCapabilities = ::SymbolicCapabilities;

class BackendRegistry {
public:
    static void register_backend(const std::string& name, Backend* backend)
    {
        ::SymbolicFactory::registerBackend(name, backend);
    }

    static bool set_backend(const std::string& name)
    {
        return ::SymbolicFactory::setBackend(name);
    }

    static std::vector<std::string> available_backends()
    {
        return ::SymbolicFactory::availableBackends();
    }

    static std::string active_backend()
    {
        return ::SymbolicFactory::activeBackendName();
    }

    static BackendCapabilities active_capabilities()
    {
        return ::SymbolicFactory::activeCapabilities();
    }
};

inline void initialize()
{
    registerPortableBackend();
}

inline bool set_backend(const std::string& name)
{
    return BackendRegistry::set_backend(name);
}

inline std::vector<std::string> available_backends()
{
    return BackendRegistry::available_backends();
}

inline std::string active_backend()
{
    return BackendRegistry::active_backend();
}

inline BackendCapabilities active_capabilities()
{
    return BackendRegistry::active_capabilities();
}

} // namespace sebsym

#endif // INCLUDE_GUARD_SEB_SYMBOLIC
