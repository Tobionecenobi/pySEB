// Symbolic scattering expressions and their recursive builders.
#include "World.hpp"
#include "ExpressionFunctions.hpp"

#include <utility>
#include <vector>

using sebsym::Expression;


namespace {

// Keep sums/products shallow without changing their mathematical content.
// This also avoids a linear-depth expression tree when many pairs contribute.
Expression combineExpressions(std::vector<Expression> terms, bool product)
{
    if (terms.empty()) return constant(product ? 1 : 0);
    while (terms.size() > 1) {
        std::vector<Expression> next;
        next.reserve((terms.size()+1)/2);
        for (std::size_t i=0; i<terms.size(); i+=2) {
            if (i+1 == terms.size()) next.push_back(terms[i]);
            else next.push_back(product ? terms[i]*terms[i+1] : terms[i]+terms[i+1]);
        }
        terms.swap(next);
    }
    return terms.front();
}

template<class Cache, class Key, class Builder>
Expression memoizedExpression(Cache& cache, const Key& key, Builder build)
{
    const auto found = cache.find(key);
    if (found != cache.end()) return found->second;
    Expression result = build();
    cache.emplace(key, result);
    return result;
}

} // namespace

/* User front end for generating Phase factor between two given reference points. Both r1 and r2 should start at the same
   level in the structure, and thus share prefix.
*/
Expression World::PhaseFactor( refPoint r1, refPoint r2, int depth, int varForm)
{
try{
    string myself=prefix(r1);
    if(myself != prefix(r2)) throw SEBException("Both reference points should be in the same structure / reference point.");
    if (depth<0)             throw SEBException("Depth can not be negative.");
    if (!hasName(myself))    throw SEBException("Unknown name "+myself+" in world.");

    if (isStructure(myself))
      {
        testPathSyntax(r1);
        testPathSyntax(r2);
      }

    betas.clear(); params.clear();

    SymbolicEvaluationContext context;
    return GenerateRefToRef(r1, r2, depth, varForm, &context);
}
catch (SEBException& e)
{
   e.PushCallStack("World::PhaseFactor(refPoint r1=\""+r1+"\", refPoint r2=\""+r2+"\", int depth="+to_string(depth)+"..)");
   throw;
}
}


// User front end for generating form factor amplitudes for a specified reference point. Returns normalized scattering expression
Expression World::FormFactorAmplitude( refPoint ref, int depth , int varForm  )
{
try{
   string myself=prefix(ref);
   if (depth<0)           throw SEBException("Depth can not be negative.");
   if (!hasName(myself))  throw SEBException("Unknown name "+myself+" in world.");
   if (!hasAPeriod(ref))  throw SEBException("String \""+ref+"\" does not specify a reference point.");
   if(isStructure(myself)) testPathSyntax(ref);

   betas.clear(); params.clear();

   SymbolicEvaluationContext context;
   return GenerateRefToAll(ref, depth, varForm, &context)/GenerateRefToAll(ref, depth, BETA, &context);
}
catch (SEBException& e)
{
   e.PushCallStack("World::FormFactorAmplitude(refPoint ref=\""+ref+"\",int depth="+to_string(depth)+", int varForm="+to_string(varForm)+")");
   throw;
}
}


// User front end for generating form factor amplitudes for a specified reference point. Returns beta weighted but does not divide by sum beta to normalize
Expression World::FormFactorAmplitude_Unnormalized( refPoint ref, int depth, int varForm)
{
try{
   string myself=prefix(ref);
   if (depth<0)            throw SEBException("Depth can not be negative.");
   if (!hasName(myself))   throw SEBException("Unknown name "+myself+" in world.");
   if (!hasAPeriod(ref))   throw SEBException("String \""+ref+"\" does not specify a reference point.");
   if(isStructure(myself)) testPathSyntax(ref);
   betas.clear(); params.clear();

   SymbolicEvaluationContext context;
   return GenerateRefToAll(ref, depth, varForm, &context);
}
catch (SEBException& e)
{
   e.PushCallStack("World::FormFactorAmplitude_Unnormalized( refPoint ref=\""+ref+"\",int depth="+to_string(depth)+", int varForm="+to_string(varForm)+",..)");
   throw;
}
}

// User front end for generating form factors everything is weighed by beta factors and normalized by (sum beta)^2
Expression World::FormFactor( structName myself, int depth, int varForm )
{
try{
   if (hasAColon(myself))  throw SEBException("Expected structure/sub-unit name got : in "+myself );
   if (hasAPeriod(myself)) throw SEBException("Expected structure/sub-unit name got . in "+myself);
   if (!hasName(myself))   throw SEBException("Unknown structure/sub-unit "+myself);
   if (depth<0)            throw SEBException("Depth can not be negative.");
   betas.clear(); params.clear();

   SymbolicEvaluationContext context;
   return GenerateAllToAll(myself, depth, varForm, &context)/GenerateAllToAll(myself, depth, BETA, &context);
}
catch (SEBException& e)
{
   e.PushCallStack("Expression World::FormFactor( structName \""+myself+"\", depth="+to_string(depth)+",..)");
   throw;
}
}

// User front end for generating form factors everything is weighed by beta factors but not normalized
Expression World::FormFactor_Unnormalized( string myself, int depth, int varForm)
{
try{
   if (hasAColon(myself))  throw SEBException("Expected structure/sub-unit name got : in "+myself);
   if (hasAPeriod(myself)) throw SEBException("Expected structure/sub-unit name got . in "+myself);
   if (!hasName(myself))   throw SEBException("Unknown structure/sub-unit "+myself);
   if (depth<0)            throw SEBException("Depth can not be negative.");
   betas.clear(); params.clear();

   SymbolicEvaluationContext context;
   return GenerateAllToAll(myself, depth, varForm, &context);
}
catch (SEBException& e)
{
   e.PushCallStack("Expression World::FormFactor_Unnormalized( structName \""+myself+"\", depth="+to_string(depth)+",..)");
   throw;
}
}

Expression World::PhaseFactorX( refPoint r1, refPoint r2, int depth)
{
    return PhaseFactor( r1, r2, depth, XVAR );
}

Expression World::PhaseFactorGeneric( refPoint r1, refPoint r2, int depth)
{
    return PhaseFactor(r1, r2, depth, GENERIC);
};




Expression World::FormFactorAmplitudeX       ( refPoint ref, int depth)
{
   return FormFactorAmplitude(ref, depth, XVAR);
}

Expression World::FormFactorAmplitudeGeneric( refPoint ref, int depth)
{
   return FormFactorAmplitude(ref, depth, GENERIC);
}

Expression World::FormFactorAmplitudeX_Unnormalized( refPoint ref, int depth)
{
   return FormFactorAmplitude_Unnormalized(ref, depth, XVAR);
}

Expression World::FormFactorAmplitudeGeneric_Unnormalized( refPoint ref, int depth)
{
   return FormFactorAmplitude_Unnormalized(ref, depth, GENERIC);
}

Expression World::FormFactorAmplitude_Normalization( refPoint ref, int depth )
{
   return FormFactorAmplitude_Unnormalized(ref, depth, BETA);
}




Expression World::FormFactorX( string name, int depth)
{
    return FormFactor(name, depth, XVAR);
}

Expression World::FormFactorGeneric(string name, int depth)
{
    return FormFactor(name, depth, GENERIC);
}

Expression World::FormFactorX_Unnormalized( string name, int depth)
{
    return FormFactor_Unnormalized(name, depth, XVAR);
}

Expression World::FormFactorGeneric_Unnormalized( string name, int depth)
{
    return FormFactor_Unnormalized(name, depth, GENERIC);
}

Expression World::FormFactor_Normalization( string name, int depth)
{
    return FormFactor_Unnormalized(name, depth, BETA);
}

Expression World::RadiusOfGyration2( string name, int depth )
{
    Expression q = GLEX->getSymbol("q");
    Expression F = FormFactor_Unnormalized( name, depth, GUINIER );

// When F terms are weigted by betas, but not normalized, the expansion is
// F = (sum_beta)^2 - q^2 Rg2/3 + ..
// Where Rg2 contributions are weighted by beta terms but not normalized
// Hence the proper normalized Rg2 is obtained by dividing the 2nd. coefficient by the zero coefficient.

    // Note: We can't use expand() and coeff() directly with our Expression class
    // This is a simplification that assumes the GUINIER form is already in the right form
    Expression Rg2 = GLEX->getSymbol("Rg2", name);
    return Rg2;
}


Expression World::SMSD_ref2scat( refPoint ref, int depth)
{
    Expression q = GLEX->getSymbol("q");
    Expression A = FormFactorAmplitude_Unnormalized( ref, depth, GUINIER );

// When A terms is weigted by betas but not normalized, the expansion is
// A = (sum_beta) - q^2 sigma <R2_ref>/6 + ..
// Where <R^2> contributions are weighted by beta terms but not normalized
// Hence the proper normalized <R2> is obtained by dividing by the zero coefficient.

    // Note: We can't use expand() and coeff() directly with our Expression class
    // This is a simplification that assumes the GUINIER form is already in the right form
    Expression sigmaMSD = GLEX->getSymbol("sigmaRrs2", prefix(ref), postfix(ref));
    return sigmaMSD;
}

Expression World::SMSD_ref2ref( refPoint r1, refPoint r2, int depth  )
{
    Expression q = GLEX->getSymbol("q");

    // Note: We can't use expand() and coeff() directly with our Expression class
    // This is a simplification that assumes the GUINIER form is already in the right form
    Expression sigmaMSD = GLEX->getSymbol("sigmaRrr2", prefix(r1), postfix(r1), postfix(r2));

    return sigmaMSD;
}


Expression World::Count(refPoint ref, int depth)
{
    return FormFactorAmplitude_Unnormalized(ref, depth, ONE).evalf();
}

Expression World::CountPairs( string name, int depth)
{
    return FormFactor_Unnormalized(name, depth, ONE).evalf();
}

// --------------------------------------------------------------------------------------------------------------------------------------------------
// Private helper functions below

/*  Returns the phase factor at a given level of nested structures between two reference points

     r1 and r2 should have the same prefix, thus start at the same structure or be on the same sub-unit.

     We can not calculate phase factors inside graphs, only if they ware wrapped in a structure.

     We should have a front-end function that tests input, and then calls a recursive helper,
     since there is no need to keep test input at every recursion.

*/
Expression World::GenerateRefToRef(refPoint r1, refPoint r2, int depth, int varForm,
                                         SymbolicEvaluationContext* context)
{
    const auto key = std::make_tuple(r1, r2, depth, varForm);
    return memoizedExpression(context->phases, key, [&]() -> Expression {
        const string myself = prefix(r1);
        if (isStructure(myself)) {
            if (LinkedForEvaluation(r1, r2, &context->topology) || r1 == r2) return constant(1);
            if (depth == 0) return getPsi(myself, postfix(r1), postfix(r2), varForm);
            if (prefix(postfix(r1)) == prefix(postfix(r2)))
                return GenerateRefToRef(postfix(r1), postfix(r2), depth-1, varForm, context);
            ReferencePointList path = CachedPath(r1, r2, false, &context->topology);
            return PhaseFactor(path, depth-1, myself, varForm, true, context);
        }
        if (isSubunit(myself))
            return getSubunit(getName(r1))->PhaseFactor(
                getReference(r1), getReference(r2), betas, params, varForm);
        throw SEBException("Unknown object in symbolic reference-to-reference traversal",
                           "World::GenerateRefToRef()");
    });
}





/* Returns the form factor amplitude at a certain level of the nested structures with respect to a given reference point

   For example  (using examples from DiBlockCopolymerStar)

   polyA.end1        => beta_polyA*A_polyA:end1  (GENERIC)

   Which is the scattering length weighted contribution of polyA relative to end1.

   diblock1:polyA.end1  => beta_polyA*A_polyA:end1 + beta_polyB*Psi_polyA:end1,end2*A_polyB:end1

   This time we want the form factor amplitude of diblock1 relative to polyA.end1 inside the structure.
   Each sub-unit contributes a scattering length weighted term, polyA.end1 is on polyA, hence there
   are no phase shifts and hence no Psi factors, to get to polyB, we first have jump accross polyA
   from end1 to end2 until we reach end1 on polyB. Hence the Psi term, and A_polyB:end1.


   star1:diblock1:polyA.end1 =>   (depth 1)
                                     beta_diblock1*A_diblock1:polyA.end1
                                    +beta_diblock2*A_diblock2:polyA.end1
                                    +beta_diblock3*A_diblock3:polyA.end1
                                    +beta_diblock4*A_diblock4:polyA.end1

   star1 contains 4 diblocks, and since they are all tethered by diblock*:polyA.end1 <-> polyA.end1 there are no phase factors.
   If we progress one level deeper (depth 2) the result is

                            =>     4*beta_polyA*A_polyA:end1
                                  +4*beta_polyB*Psi_polyA:end1,end2*A_polyB:end1

   Each diblock copolymer has an A block close to the centre, and a B block away from the centre. Hence there are no phase
   factors for the A_polyA:end1 which are directly connected to the reference point, whereas the to get from the reference
   point to a polyB we have to jump across an A block first, which encurs a Psi_polyA:end1,end2 phase factor.

   Same structure and level but a different reference point this time:

   star1:diblock1:polyB.end2   =>       beta_polyB*A_polyB:end2
                                       +beta_polyA*Psi_polyB:end1,end2*A_polyA:end2
                                      +3*Psi_polyB:end1,end2*Psi_polyA:end1,end2 * beta_polyA*A_polyA:end1
                                      +3*Psi_polyB:end1,end2*Psi_polyA:end1,end2 * Psi_polyA:end1,end2* A_polyB:end1*beta_polyB)

   Here diblock1:polyB is directly connected to the reference point producing the first term.
   and diblock1:polyA can be reached by jumping across polyB from end1 to end2, thus producing the second term.
   furtermore, there are 3 polyA sub-units in diblock[2-4] but to reach these, we first have to jump across polyB and polyA in diblock1, thus producing the third term.
   finally, there are 3 polyB sub-units in diblock[2-4] but to reach these, we have to jump across two polyA and one polyB thus producing the final term.




*/
Expression World::GenerateRefToAll(refPoint ref, int depth, int varForm,
                                         SymbolicEvaluationContext* context)
{
    const auto key = std::make_tuple(ref, depth, varForm);
    return memoizedExpression(context->amplitudes, key, [&]() -> Expression {
        const string myself = prefix(ref);
        if (isStructure(myself)) {
            if (depth == 0) return getFFA(myself, postfix(ref), varForm);
            std::vector<Expression> terms;
            const GraphID gid = getStructure(myself)->getGraphID();
            for (auto child=subgraph_cbegin(gid); child!=subgraph_cend(gid); ++child) {
                ReferencePointList path = CachedPath(ref, myself+":"+*child, true, &context->topology);
                if (path.empty()) terms.push_back(GenerateRefToAll(postfix(ref), depth-1, varForm, context));
                else terms.push_back(PhaseFactor(path, depth-1, myself, varForm, true, context) *
                                     GenerateRefToAll(path.back(), depth-1, varForm, context));
            }
            return combineExpressions(std::move(terms), false);
        }
        if (isSubunit(myself))
            return getSubunit(myself)->FormFactorAmplitude(getReference(ref), betas, params, varForm);
        throw SEBException("Unknown object in symbolic reference-to-all traversal",
                           "World::GenerateRefToAll()");
    });
}


/*    Calculates form factor of a structure or sub-unit within the world.

*/
Expression World::GenerateAllToAll(string myself, int depth, int varForm,
                                         SymbolicEvaluationContext* context)
{
    const auto key = std::make_tuple(myself, depth, varForm);
    return memoizedExpression(context->intensities, key, [&]() -> Expression {
        if (isStructure(myself)) {
            if (depth == 0) return getFF(myself, varForm);
            std::vector<Expression> terms;
            const GraphID gid = getStructure(myself)->getGraphID();
            for (auto first=subgraph_cbegin(gid); first!=subgraph_cend(gid); ++first) {
                for (auto second=first; second!=subgraph_cend(gid); ++second) {
                    if (first == second) {
                        terms.push_back(GenerateAllToAll(*first, depth-1, varForm, context));
                    } else {
                        ReferencePointList path = CachedPath(
                            myself+":"+*first, myself+":"+*second, false, &context->topology);
                        Expression a1 = GenerateRefToAll(path.front(), depth-1, varForm, context);
                        Expression phase = PhaseFactor(path, depth-1, myself, varForm, true, context);
                        Expression a2 = GenerateRefToAll(path.back(), depth-1, varForm, context);
                        terms.push_back(constant(2)*a1*phase*a2);
                    }
                }
            }
            return combineExpressions(std::move(terms), false);
        }
        if (isSubunit(myself)) return getSubunit(myself)->FormFactor(betas, params, varForm);
        throw SEBException("Unknown object in symbolic all-to-all traversal",
                           "World::GenerateAllToAll()");
    });
}

Expression World::PhaseFactor(ReferencePointList& path, int depth, string myself,
                              int varForm, bool doCheck, SymbolicEvaluationContext* context)
{
    if (depth < 0) throw SEBException("Depth can not be negative.", "World::PhaseFactor(path)");
    // Also support existing callers of this helper without a supplied context.
    if (!context) {
        SymbolicEvaluationContext local;
        return PhaseFactor(path, depth, myself, varForm, doCheck, &local);
    }
    std::vector<Expression> factors;
    if (path.empty()) return constant(1);
    auto current=path.begin(), next=current;
    ++next;
    for (; next!=path.end(); ++current, ++next) {
        if (LinkedForEvaluation(*current, *next, &context->topology)) continue;
        if (depth == 0) {
            if (hasAColon(*current))
                factors.push_back(getPsi(prefix(*current), postfix(*current), postfix(*next), varForm));
            else factors.push_back(getPsi(getName(*current), getReference(*current), getReference(*next), varForm));
        } else factors.push_back(GenerateRefToRef(*current, *next, depth, varForm, context));
    }
    return combineExpressions(std::move(factors), true);
}

// Code for generating a GENERIC Psi for a structure depending on various varforms
Expression World::getPsi(string myself, refPoint r1, refPoint r2, int varform)
{
     if      (varform == BETA)     return constant(1);
     else if (varform == GUINIER)
                                   {
                                      Expression q = GLEX->getSymbol("q");
                                      Expression R2 = GLEX->getSymbol("sigmaRrr2", myself, r1, r2);
                                      params[R2.to_string()]=0;
                                      return constant(1) - q * q * R2 / constant(6);
                                   }
     else if (varform == ONE)      return constant(1);
     else
                                   {
                                      Expression psi = GLEX->getSymbol("Psi", myself, r1, r2);
                                      params[psi.to_string()]=0;
                                      return psi;
                                   }
}

// Code for generating a GENERIC Form Factor Amplitude for a structure depending on various varforms
Expression World::getFFA(string myself, refPoint r, int varform)
{
     if      (varform == BETA)
                                   {
                                      Expression beta = GLEX->getSymbol("beta", myself);
                                      betas[beta.to_string()]=0;
                                      return beta;
                                   }
     else if (varform == GUINIER)  {
                                      Expression beta = GLEX->getSymbol("beta", myself);
                                      Expression q = GLEX->getSymbol("q");
                                      Expression R2 = GLEX->getSymbol("sigmaRrs2", myself, r);

                                      params[R2.to_string()]=0;
                                      betas[beta.to_string()]=0;

                                      return beta * (constant(1) - q * q * R2 / constant(6));
                                   }
     else if (varform == ONE)      return constant(1);
     else                          {
                                      Expression beta = GLEX->getSymbol("beta", myself);
                                      Expression A = GLEX->getSymbol("A", myself, r );

                                      betas[beta.to_string()]=0;
                                      params[A.to_string()]=0;

                                      return beta * A;
                                   }
}

// Code for generating a GENERIC Form Factor for a structure depending on various varforms
Expression World::getFF(string myself, int varform)
{
     if      (varform == BETA)     {
                                      Expression beta = GLEX->getSymbol("beta", myself);
                                      betas[beta.to_string()]=0;
                                      return beta * beta;
                                   }
     else if (varform == GUINIER)  {
                                      Expression beta = GLEX->getSymbol("beta", myself);
                                      Expression q = GLEX->getSymbol("q");
                                      Expression Rg2 = GLEX->getSymbol("Rg2", myself);

                                      betas[beta.to_string()]=0;
                                      params[Rg2.to_string()]=0;

                                      return beta * beta * (constant(1) - q * q * Rg2 / constant(3));
                                   }
     else if (varform == ONE)      return constant(1);
     else                          {
                                      Expression beta = GLEX->getSymbol("beta", myself);
                                      Expression F = GLEX->getSymbol("F", myself );

                                      betas[beta.to_string()]=0;
                                      params[F.to_string()]=0;

                                      return beta * beta * F;
                                   }
}
