#ifndef PYSEB_SUBUNITS_FILE_DEFINED_HPP
#define PYSEB_SUBUNITS_FILE_DEFINED_HPP

#include <map>
#include <memory>

#include "SubunitIO/SubunitDefinition.hpp"
#include "Subunits/Integrated.hpp"

class FileDefinedSubunit : public IntegratedSubunit {
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

private:
    std::map<std::string, std::unique_ptr<NumericalIntegrator>> integrators_;

    bool requiresParsedNumericalEvaluation(
        const pyseb::ParsedExpression& expression) const;
    double evaluateNumerically(
        const pyseb::ParsedExpression& expression,
        double q,
        const ParameterList& values);
};

#endif
