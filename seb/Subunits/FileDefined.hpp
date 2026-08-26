#ifndef PYSEB_SUBUNITS_FILE_DEFINED_HPP
#define PYSEB_SUBUNITS_FILE_DEFINED_HPP

#include "Subunit.hpp"
#include "SubunitIO/SubunitDefinition.hpp"

class FileDefinedSubunit : public SubUnit {
public:
    explicit FileDefinedSubunit(const pyseb::SubunitDefinition& definition);
    ~FileDefinedSubunit() override = default;

    void Init(string name, string tag, SymbolInterface* GEX) override;

    Expression FormFactor(
        ParameterList& betas,
        ParameterList& params,
        int varForm = GENERIC) override;
    Expression FormFactorAmplitude(
        refPoint reference,
        ParameterList& betas,
        ParameterList& params,
        int varForm = GENERIC) override;
    double NumericTotalBeta(const ParameterList& values) override;

    const pyseb::SubunitDefinition& getDefinition() const { return definition_; }

protected:
    pyseb::SubunitDefinition definition_;
};

#endif
