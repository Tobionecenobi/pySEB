// Parameter handling, numerical scattering, and Guinier coefficients.
#include "World.hpp"
#include "ExpressionFunctions.hpp"

#include <algorithm>
#include <cmath>

using sebsym::Expression;


// Methods for handling parameters, and evaluating expressions.



// Returns a copy of the parameter lists.
ParameterList World::getParams()
{
   ParameterList pl(betas);
   pl.insert(params.begin(), params.end());
   return pl;
}

// Returns a copy of the parameter lists.
ParameterList World::getParamsq()
{
   ParameterList pl(betas);
   pl.insert(params.begin(), params.end());
   Expression q=GLEX->getSymbol("q");
   pl[q.to_string()]=0;
   return pl;
}

ParameterList World::getq(double q) {
   ParameterList pl;
   Expression Q=GLEX->getSymbol("q");
   pl[Q.to_string()]=q;
   return pl;
}

ParameterList World::q(double q) { return getq(q); }

double World::setParameter(ParameterList& pl, string str, double x)
{
   // Throw if str not known??
   Expression e=GLEX->get(str);
   pl[e.to_string()]=x;
   return x;
}


DoubleVector World::linspace(double q1,double q2,int n)
{
   if (q2<q1) swap(q1,q2);
   if (q1 <=0) throw SEBException("Negative values not allowed","World::linspace("+to_string(q1)+","+to_string(q2)+","+to_string(n)+ ")");
   if (q2 <=0) throw SEBException("Negative values not allowed","World::linspace("+to_string(q1)+","+to_string(q2)+","+to_string(n)+ ")");
   if (n <=1)  throw SEBException("Number of points must be larger than 1","World::linspace("+to_string(q1)+","+to_string(q2)+","+to_string(n)+ ")");

   DoubleVector q;
   for (int i=0;i<n;i++)
         q.push_back( (q2-q1)*i/(n-1) + q1 );

   return q;
}

DoubleVector World::logspace(double q1,double q2,int n)
{
   if (q2<q1) swap(q1,q2);
   if (q1 <=0) throw SEBException("Negative values not allowed","World::logspace("+to_string(q1)+","+to_string(q2)+","+to_string(n)+ ")");
   if (q2 <=0) throw SEBException("Negative values not allowed","World::logspace("+to_string(q1)+","+to_string(q2)+","+to_string(n)+ ")");
   if (n <=1)  throw SEBException("Number of points must be larger than 1","World::logspace("+to_string(q1)+","+to_string(q2)+","+to_string(n)+ ")");

   DoubleVector q;
   for (int i=0;i<n;i++)
       {
         q.push_back( pow( 10.0,  log10(q2/q1)*i/(n-1) + log10(q1)  ) );
       }
   return q;
}

DoubleVector World::logspace(double q1,double q2, double base)
{
   if (q2<q1) swap(q1,q2);
   if (q1 <=0) throw SEBException("Negative values not allowed","World::logspace("+to_string(q1)+","+to_string(q2)+","+to_string(base)+ ")");
   if (q2 <=0) throw SEBException("Negative values not allowed","World::logspace("+to_string(q1)+","+to_string(q2)+","+to_string(base)+ ")");
   if (base<=1)  throw SEBException("Base must be larger than one","World::logspace("+to_string(q1)+","+to_string(q2)+","+to_string(base)+ ")");

   DoubleVector q;
   double x=q1;
   while (x<=q2)
     {
        q.push_back(x);
        x*=base;
     }

   return q;
}

double World::Evaluate(Expression e, ParameterList& pl)
{
   Expression result = e.subs(pl).evalf();

   // Check if the expression is numeric
   if (result.is_numeric())
        return result.eval();
    else
        throw SEBException("Expression did not evaluate to number, since not all parameters were specified", "World::Evaluate(Expression,ParameterList)");
}

double World::Evaluate(Expression e, ParameterList& pl, double q)
{
   Expression result = e.subs(pl).subs(getq(q)).evalf();

   // Check if the expression is numeric
   if (result.is_numeric())
        return result.eval();
    else
        throw SEBException("Expression did not evaluate to number, since not all parameters were specified", "World::Evaluate(Expression,ParameterList, double)");
}

DoubleVector World::Evaluate(Expression e, ParameterList& pl, DoubleVector& q)
{
   DoubleVector I(q.size());
   Expression eval = e.subs(pl);
   for (int i=0; i<q.size() ; i++)
      {
        Expression result = eval.subs(getq(q[i])).evalf();

        if (result.is_numeric()) I[i] = result.eval();
        else
        throw SEBException("Expression did not evaluate to number, since not all parameters were specified", "World::Evaluate(Expression,ParameterList, DoubleVector)");
      }

   return I;
}

DoubleVector World::Evaluate(Expression e, ParameterList& pl, DoubleVector& q, string fname, string comment, string prefix)
{
   DoubleVector I=Evaluate(e,pl,q);

   ofstream fo(fname);
   if (!fo.is_open()) throw SEBException("Could not open file "+fname, "World::Evaluate(Expression, Parameterlist&, DoubleVector&, string "+fname+" )");

   fo << prefix << "File geneated by SEB\n";
   fo << prefix << comment << "\n";
   fo << prefix << "Expression = " << e << "\n";
   fo << prefix << "Parameters:\n";
   for (auto const& p : pl)
       fo << prefix << "\t" << p.first << "=" << p.second << "\n";

   for (int i=0;i<q.size(); i++)
      fo << q[i] << " " << I[i] << "\n";
   fo.close();

   return I;
}

double World::EvaluateFormFactor(
    string name,
    const ParameterList& values,
    double q)
{
    const double beta = GenerateTotalBetaNumeric(name, values);
    if (std::abs(beta) <= 1e-14) {
        throw SEBException(
            "Cannot normalize a form factor when the total excess scattering length is zero",
            "World::EvaluateFormFactor()"
        );
    }
    return GenerateAllToAllNumeric(name, q, values) / (beta * beta);
}

DoubleVector World::EvaluateFormFactor(
    string name,
    const ParameterList& values,
    const DoubleVector& q)
{
    DoubleVector result;
    result.reserve(q.size());
    if (q.empty()) return result;
    const double beta = GenerateTotalBetaNumeric(name, values);
    if (std::abs(beta) <= 1e-14) {
        throw SEBException(
            "Cannot normalize a form factor when the total excess scattering length is zero",
            "World::EvaluateFormFactor()"
        );
    }
    NumericEvaluationContext context;
    for (const double value : q) {
        context.beginSample();
        result.push_back(GenerateAllToAllNumeric(name, value, values, &context) / (beta * beta));
    }
    return result;
}

double World::EvaluateFormFactorUnnormalized(
    string name,
    const ParameterList& values,
    double q)
{
    return GenerateAllToAllNumeric(name, q, values);
}

DoubleVector World::EvaluateFormFactorUnnormalized(
    string name,
    const ParameterList& values,
    const DoubleVector& q)
{
    DoubleVector result;
    result.reserve(q.size());
    NumericEvaluationContext context;
    for (const double value : q) {
        context.beginSample();
        result.push_back(GenerateAllToAllNumeric(name, value, values, &context));
    }
    return result;
}

std::pair<DoubleVector, DoubleVector> World::EvaluateFormFactorAndAmplitudeUnnormalized(
    string name, refPoint reference, const ParameterList& values, const DoubleVector& q)
{
    if (prefix(reference) != name) {
        throw SEBException("Amplitude reference must belong to the intensity target",
                           "World::EvaluateFormFactorAndAmplitudeUnnormalized()");
    }
    std::pair<DoubleVector, DoubleVector> result;
    result.first.reserve(q.size());
    result.second.reserve(q.size());
    NumericEvaluationContext context;
    for (const double value : q) {
        context.beginSample();
        // The complete intensity includes internal interference. Its leaf
        // amplitudes can also serve the independent reference-amplitude query.
        result.first.push_back(GenerateAllToAllNumeric(name, value, values, &context));
        result.second.push_back(GenerateRefToAllNumeric(reference, value, values, &context));
    }
    return result;
}

double World::EvaluateFormFactorAmplitude(
    refPoint ref,
    const ParameterList& values,
    double q)
{
    const double beta = GenerateTotalBetaNumeric(prefix(ref), values);
    if (std::abs(beta) <= 1e-14) {
        throw SEBException(
            "Cannot normalize a form-factor amplitude when the total excess scattering length is zero",
            "World::EvaluateFormFactorAmplitude()"
        );
    }
    return GenerateRefToAllNumeric(ref, q, values) / beta;
}

DoubleVector World::EvaluateFormFactorAmplitude(
    refPoint ref,
    const ParameterList& values,
    const DoubleVector& q)
{
    DoubleVector result;
    result.reserve(q.size());
    if (q.empty()) return result;
    const double beta = GenerateTotalBetaNumeric(prefix(ref), values);
    if (std::abs(beta) <= 1e-14) {
        throw SEBException(
            "Cannot normalize a form-factor amplitude when the total excess scattering length is zero",
            "World::EvaluateFormFactorAmplitude()"
        );
    }
    NumericEvaluationContext context;
    for (const double value : q) {
        context.beginSample();
        result.push_back(GenerateRefToAllNumeric(ref, value, values, &context) / beta);
    }
    return result;
}

double World::EvaluateFormFactorAmplitudeUnnormalized(
    refPoint ref,
    const ParameterList& values,
    double q)
{
    return GenerateRefToAllNumeric(ref, q, values);
}

DoubleVector World::EvaluateFormFactorAmplitudeUnnormalized(
    refPoint ref,
    const ParameterList& values,
    const DoubleVector& q)
{
    DoubleVector result;
    result.reserve(q.size());
    NumericEvaluationContext context;
    for (const double value : q) {
        context.beginSample();
        result.push_back(GenerateRefToAllNumeric(ref, value, values, &context));
    }
    return result;
}

double World::EvaluatePhaseFactor(
    refPoint r1,
    refPoint r2,
    const ParameterList& values,
    double q)
{
    return GenerateRefToRefNumeric(r1, r2, q, values);
}

DoubleVector World::EvaluatePhaseFactor(
    refPoint r1,
    refPoint r2,
    const ParameterList& values,
    const DoubleVector& q)
{
    DoubleVector result;
    result.reserve(q.size());
    NumericEvaluationContext context;
    for (const double value : q) {
        context.beginSample();
        result.push_back(GenerateRefToRefNumeric(r1, r2, value, values, &context));
    }
    return result;
}

double World::EvaluateRadiusOfGyration2(
    string name,
    const ParameterList& values)
{
    const GuinierTerm term = GenerateAllToAllGuinier(name, values);
    if (std::abs(term.constant) <= 1e-14) {
        throw SEBException(
            "Cannot calculate radius of gyration when the total excess scattering length is zero",
            "World::EvaluateRadiusOfGyration2()"
        );
    }
    return 3.0 * term.q2Coefficient / term.constant;
}

double World::EvaluateSMSDRef2Scat(
    refPoint ref,
    const ParameterList& values)
{
    const GuinierTerm term = GenerateRefToAllGuinier(ref, values);
    if (std::abs(term.constant) <= 1e-14) {
        throw SEBException(
            "Cannot calculate reference-to-scatterer sigma<R^2> when total beta is zero",
            "World::EvaluateSMSDRef2Scat()"
        );
    }
    return 6.0 * term.q2Coefficient / term.constant;
}

double World::EvaluateSMSDRef2Ref(
    refPoint r1,
    refPoint r2,
    const ParameterList& values)
{
    const GuinierTerm term = GenerateRefToRefGuinier(r1, r2, values);
    if (std::abs(term.constant) <= 1e-14) {
        throw SEBException(
            "Reference-to-reference Guinier normalization is zero",
            "World::EvaluateSMSDRef2Ref()"
        );
    }
    return 6.0 * term.q2Coefficient / term.constant;
}

double World::GenerateTotalBetaNumeric(
    string name,
    const ParameterList& values)
{
    if (isSubunit(name)) {
        return getSubunit(name)->NumericTotalBeta(values);
    }

    if (isStructure(name)) {
        double beta = 0.0;
        const GraphID gid = getStructure(name)->getGraphID();
        for (auto child = subgraph_cbegin(gid); child != subgraph_cend(gid); ++child) {
            beta += GenerateTotalBetaNumeric(*child, values);
        }
        return beta;
    }

    throw SEBException(
        "Unknown sub-unit or structure " + name,
        "World::GenerateTotalBetaNumeric()"
    );
}

bool World::LinkedForEvaluation(
    const refPoint& first, const refPoint& second, NumericEvaluationContext* context)
{
    if (links.empty()) return false;
    // Bound optional topology acceleration and retain the original path for
    // calls without a context or worlds too large for this local cache.
    if (!context || links.size() > 4096) return isLinked(first, second);
    if (!context->linksReady) {
        context->linkedReferences.assign(links.begin(), links.end());
        std::sort(context->linkedReferences.begin(), context->linkedReferences.end());
        context->linksReady = true;
    }
    const auto& low = second < first ? second : first;
    const auto& high = second < first ? first : second;
    // Compare references directly: nested names can exceed string's inline
    // capacity, so constructing an owning lookup key would allocate per edge.
    const auto key = std::tie(low, high);
    const auto found = std::lower_bound(context->linkedReferences.begin(),
        context->linkedReferences.end(), key, [](const auto& link, const auto& value) {
            return std::tie(link.first, link.second) < value;
        });
    return found != context->linkedReferences.end() && found->first == low && found->second == high;
}

ReferencePointList World::CachedPath(
    const string& first, const string& second, bool check,
    NumericEvaluationContext* context)
{
    if (!context) return findPath(first, second, check, nullptr);
    const auto key = std::make_tuple(first, second, check);
    const auto found = context->paths.find(key);
    if (found != context->paths.end()) return found->second;
    ReferencePointList path = findPath(first, second, check, context);
    // Long-chain pair paths can require cubic storage. Bound this optional
    // acceleration and fall back to uncached search when the shared path and
    // adjacency budget is full.
    if (context->paths.size() < 4096 &&
        path.size() <= 262144 - context->referenceEntries) {
        context->paths.emplace(key, path);
        context->referenceEntries += path.size();
    }
    return path;
}

double World::GenerateRefToRefNumeric(
    refPoint r1,
    refPoint r2,
    double q,
    const ParameterList& values,
    NumericEvaluationContext* context)
{
    if (!std::isfinite(q)) {
        throw SEBException("q must be finite", "World::GenerateRefToRefNumeric()");
    }

    const string myself = prefix(r1);
    if (myself != prefix(r2)) {
        throw SEBException(
            "Both reference points must start at the same structure level",
            "World::GenerateRefToRefNumeric()"
        );
    }

    if (isStructure(myself)) {
        if (LinkedForEvaluation(r1, r2, context) || r1 == r2) {
            return 1.0;
        }

        if (prefix(postfix(r1)) == prefix(postfix(r2))) {
            return GenerateRefToRefNumeric(postfix(r1), postfix(r2), q, values, context);
        }

        ReferencePointList path = CachedPath(r1, r2, false, context);
        return PhaseFactorNumeric(path, q, values, context);
    }

    if (isSubunit(myself)) {
        const auto key = std::make_pair(r1, r2);
        if (context) {
            const auto found = context->phases.find(key);
            if (found != context->phases.end()) return found->second;
        }
        const double result = getSubunit(getName(r1))->NumericPhaseFactor(
            getReference(r1), getReference(r2), q, values);
        if (context && context->phases.size() < 4096)
            context->phases.emplace(key, result);
        return result;
    }

    throw SEBException(
        "Unknown object in numerical reference-to-reference traversal",
        "World::GenerateRefToRefNumeric()"
    );
}

double World::GenerateRefToAllNumeric(
    refPoint ref,
    double q,
    const ParameterList& values,
    NumericEvaluationContext* context)
{
    if (!std::isfinite(q)) {
        throw SEBException("q must be finite", "World::GenerateRefToAllNumeric()");
    }

    const string myself = prefix(ref);
    if (isStructure(myself)) {
        double amplitude = 0.0;
        const GraphID gid = getStructure(myself)->getGraphID();

        for (auto child = subgraph_cbegin(gid); child != subgraph_cend(gid); ++child) {
            ReferencePointList path = CachedPath(ref, myself + ":" + *child, true, context);
            double term = 0.0;
            if (path.empty()) {
                term = GenerateRefToAllNumeric(postfix(ref), q, values, context);
            } else {
                term = PhaseFactorNumeric(path, q, values, context) *
                       GenerateRefToAllNumeric(path.back(), q, values, context);
            }
            amplitude += term;
        }
        return amplitude;
    }

    if (isSubunit(myself)) {
        if (context) {
            const auto found = context->amplitudes.find(ref);
            if (found != context->amplitudes.end()) return found->second;
        }
        const double result = getSubunit(myself)->NumericFormFactorAmplitudeUnnormalized(
            getReference(ref), q, values);
        if (context && context->amplitudes.size() < 4096)
            context->amplitudes.emplace(ref, result);
        return result;
    }

    throw SEBException(
        "Unknown object in numerical reference-to-all traversal",
        "World::GenerateRefToAllNumeric()"
    );
}

double World::GenerateAllToAllNumeric(
    string name,
    double q,
    const ParameterList& values,
    NumericEvaluationContext* context)
{
    if (!std::isfinite(q)) {
        throw SEBException("q must be finite", "World::GenerateAllToAllNumeric()");
    }

    if (isSubunit(name)) {
        return getSubunit(name)->NumericFormFactorUnnormalized(q, values);
    }

    if (isStructure(name)) {
        double formFactor = 0.0;
        const GraphID gid = getStructure(name)->getGraphID();

        for (auto child1 = subgraph_cbegin(gid); child1 != subgraph_cend(gid); ++child1) {
            for (auto child2 = child1; child2 != subgraph_cend(gid); ++child2) {
                if (child1 == child2) {
                    formFactor += GenerateAllToAllNumeric(*child1, q, values, context);
                    continue;
                }

                ReferencePointList path = CachedPath(
                    name + ":" + *child1,
                    name + ":" + *child2,
                    false, context
                );
                const double amplitude1 = GenerateRefToAllNumeric(path.front(), q, values, context);
                const double phase = PhaseFactorNumeric(path, q, values, context);
                const double amplitude2 = GenerateRefToAllNumeric(path.back(), q, values, context);
                formFactor += 2.0 * amplitude1 * phase * amplitude2;
            }
        }
        return formFactor;
    }

    throw SEBException(
        "Unknown object in numerical all-to-all traversal",
        "World::GenerateAllToAllNumeric()"
    );
}

double World::PhaseFactorNumeric(
    ReferencePointList& path,
    double q,
    const ParameterList& values,
    NumericEvaluationContext* context)
{
    if (path.size() < 2) {
        return 1.0;
    }

    double phase = 1.0;
    auto current = path.begin();
    auto next = current;
    ++next;
    while (next != path.end()) {
        if (!LinkedForEvaluation(*current, *next, context)) {
            phase *= GenerateRefToRefNumeric(*current, *next, q, values, context);
        }
        current = next;
        ++next;
    }
    return phase;
}

World::GuinierTerm World::AddGuinier(
    const GuinierTerm& left,
    const GuinierTerm& right)
{
    return GuinierTerm(
        left.constant + right.constant,
        left.q2Coefficient + right.q2Coefficient
    );
}

World::GuinierTerm World::MultiplyGuinier(
    const GuinierTerm& left,
    const GuinierTerm& right)
{
    return GuinierTerm(
        left.constant * right.constant,
        left.q2Coefficient * right.constant +
            left.constant * right.q2Coefficient
    );
}

World::GuinierTerm World::ScaleGuinier(
    double scale,
    const GuinierTerm& term)
{
    return GuinierTerm(
        scale * term.constant,
        scale * term.q2Coefficient
    );
}

World::GuinierTerm World::GenerateRefToRefGuinier(
    refPoint r1,
    refPoint r2,
    const ParameterList& values)
{
    const string myself = prefix(r1);
    if (myself != prefix(r2)) {
        throw SEBException(
            "Both reference points must start at the same structure level",
            "World::GenerateRefToRefGuinier()"
        );
    }

    if (isStructure(myself)) {
        if (isLinked(r1, r2) || r1 == r2) {
            return GuinierTerm(1.0, 0.0);
        }

        if (prefix(postfix(r1)) == prefix(postfix(r2))) {
            return GenerateRefToRefGuinier(postfix(r1), postfix(r2), values);
        }

        ReferencePointList path = findpath(r1, r2, false);
        return PhaseFactorGuinier(path, values);
    }

    if (isSubunit(myself)) {
        SubUnit* subunit = getSubunit(getName(r1));
        const double sigma = subunit->NumericSigmaMSDRef2Ref(
            getReference(r1),
            getReference(r2),
            values
        );
        return GuinierTerm(1.0, sigma / 6.0);
    }

    throw SEBException(
        "Unknown object in Guinier reference-to-reference traversal",
        "World::GenerateRefToRefGuinier()"
    );
}

World::GuinierTerm World::GenerateRefToAllGuinier(
    refPoint ref,
    const ParameterList& values)
{
    const string myself = prefix(ref);
    if (isStructure(myself)) {
        GuinierTerm amplitude;
        const GraphID gid = getStructure(myself)->getGraphID();

        for (auto child = subgraph_cbegin(gid); child != subgraph_cend(gid); ++child) {
            ReferencePointList path = findpath(ref, myself + ":" + *child);
            GuinierTerm term;
            if (path.empty()) {
                term = GenerateRefToAllGuinier(postfix(ref), values);
            } else {
                term = MultiplyGuinier(
                    PhaseFactorGuinier(path, values),
                    GenerateRefToAllGuinier(path.back(), values)
                );
            }
            amplitude = AddGuinier(amplitude, term);
        }
        return amplitude;
    }

    if (isSubunit(myself)) {
        SubUnit* subunit = getSubunit(myself);
        const refPoint reference = getReference(ref);
        const double beta = subunit->NumericTotalBeta(values);
        const double sigma = subunit->NumericSigmaMSDRef2Scat(reference, values);
        return GuinierTerm(beta, beta * sigma / 6.0);
    }

    throw SEBException(
        "Unknown object in Guinier reference-to-all traversal",
        "World::GenerateRefToAllGuinier()"
    );
}

World::GuinierTerm World::GenerateAllToAllGuinier(
    string name,
    const ParameterList& values)
{
    if (isSubunit(name)) {
        SubUnit* subunit = getSubunit(name);
        const double beta = subunit->NumericTotalBeta(values);
        const double rg2 = subunit->NumericRadiusOfGyration2(values);
        return GuinierTerm(beta * beta, beta * beta * rg2 / 3.0);
    }

    if (isStructure(name)) {
        GuinierTerm formFactor;
        const GraphID gid = getStructure(name)->getGraphID();

        for (auto child1 = subgraph_cbegin(gid); child1 != subgraph_cend(gid); ++child1) {
            for (auto child2 = child1; child2 != subgraph_cend(gid); ++child2) {
                GuinierTerm term;
                if (child1 == child2) {
                    term = GenerateAllToAllGuinier(*child1, values);
                } else {
                    ReferencePointList path = findpath(
                        name + ":" + *child1,
                        name + ":" + *child2,
                        false
                    );
                    term = ScaleGuinier(
                        2.0,
                        MultiplyGuinier(
                            MultiplyGuinier(
                                GenerateRefToAllGuinier(path.front(), values),
                                PhaseFactorGuinier(path, values)
                            ),
                            GenerateRefToAllGuinier(path.back(), values)
                        )
                    );
                }
                formFactor = AddGuinier(formFactor, term);
            }
        }
        return formFactor;
    }

    throw SEBException(
        "Unknown object in Guinier all-to-all traversal",
        "World::GenerateAllToAllGuinier()"
    );
}

World::GuinierTerm World::PhaseFactorGuinier(
    ReferencePointList& path,
    const ParameterList& values)
{
    GuinierTerm phase(1.0, 0.0);
    if (path.size() < 2) {
        return phase;
    }

    auto current = path.begin();
    auto next = current;
    ++next;
    while (next != path.end()) {
        if (!isLinked(*current, *next)) {
            phase = MultiplyGuinier(
                phase,
                GenerateRefToRefGuinier(*current, *next, values)
            );
        }
        current = next;
        ++next;
    }
    return phase;
}
