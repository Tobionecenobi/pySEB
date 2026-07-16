//===========================================================================
// Included guards
#ifndef INCLUDE_GUARD_WORLD
#define INCLUDE_GUARD_WORLD

//===========================================================================
// included dependencies
#include <iostream>
#include <fstream>
#include <cmath>

#include "Types.hpp"
#include "Constants.hpp"
#include "SymbolInterface.hpp"
#include "SpecialFunctions.hpp"

#include "Structure.hpp"
#include "Subunit.hpp"
#include "Subunits/Subunits.hpp"


//===========================================================================
// used namespaces
using namespace std;

class World
{
private:

    /* Name of world, not really used anywhere. */
    string worldId;

    /* links between reference points tying sub-units to sub-units or structures to structures.

           links are either
                pairs of subunit_name.referencepoint when sub-units are directly connected,
             or pairs of paths of the form topstructurename:structurename:..:subunitname.referencepoint

           the latter explains how topstructurename is connected to something else, by specifying a path
           through  nested structures until a sub-unit, and finally a reference point is met.

           Thus LHS:RHS can be thought of as a defining LHS as containing a RHS inside it.
    */
    list<LinkPair> links;

    /* storage of all sub-units/structures known by world by their unique name. Here names not paths, so no : anywhere. */
    map<string, ABSSubUnit*> nameCatalog;

    /* Each connected graph of sub-units or structures has a unique graphID. When a new sub-unit is linked to an
    existing sub-unit, it automagically shares its graphid, similar when a structure is linked to a structure */
    int totalNumberofGraphs = 0;

    /* Map listing all the sub-unit / structure names inside each sub-graph. */
    map<GraphID, list<string>> subGraphs;

    /*   sub-unit name -> type map.  WHY THIS WHEN sub-units knows their own type? */
    map<subName, int> typeCatalog;

    // ?? for optimization?
    map<pair<string, string>, list<refPoint>> alreadyKnownPaths;
    int numberOfPathsFound = 0; //For print to the user.

    /* Owned symbol interface for this World instance */
    SymbolInterface *GLEX = new SymbolInterface();

    // Parameter lists
    ParameterList betas, params;

public:
    /*World Constructer with a given id as a string*/
    World(string id="World")
    {
        worldId = id;
#ifdef SPECIALFUNCTIONSERIES
        cout << "WARNING: This was compiled with -DSPECIALFUNCTIONSERIES, which should only be used when validating the sub-unit scattering expressions, NOT during production.\n";
#endif
    };

    ~World(){
        // release allocated sub-units / structures.
        for (auto it = nameCatalog.begin(); it != nameCatalog.end(); ++it) delete it->second;

        // release this World's symbol interface
        delete GLEX;
    };

    /* Get graphID from a sub-unit or structure name / sub-unit pointer */
    GraphID getGraphID(string);              // NB. fails if given a path as structure_name:subunit_name, since that's two different graphIDs!
    GraphID getGraphID(SubUnit*);            // Supplied for user convenience. NOT USED INTERNALLY
    GraphID getGraphID(Structure*);          // Supplied for user convenience. NOT USED INTERNALLY

    // Get pointers from strings
    SubUnit* getSubunit(string);             // return sub-unit pointer to specified sub-unit.
    Structure* getStructure(string);         // return structure pointer to specified structure

    bool isSubunit(string);                  // Queries whether a string is a sub-unit name. Throws if string does not not exist.
    bool isStructure(string);                // Queries whether a string is a structure name.  Throws if string does not not exist.

    /* Is name already known to this world. */
    bool hasName(string);

    // Specified sub-graph gid, contains name. Throws if gid is not known.
    bool doesSubgraphContainName(int, string);

    // Test that parent structures contains the specified sibling. (searches one level down, does not recursively search if sibling is somewhere inside parents siblings)
    bool doesStructureContainName(string parent,string sibling);

    // returns true if r1 <-> r2.
    bool isLinked(refPoint r1, refPoint r2);

    // print structure as folder hierarchy.
    void folderprint(string, string  ="", bool =false);

    // Print links known currently
    void printLinks();

    // Print links known currently, potentially specify prefix and postfix strings.
    void printPath(ReferencePointList, string ="", string = "\n");

    // print out graphIDs
    void printGraphIDs();

    // print paths.
    void printList(ReferencePointList, string ="", string="\n");

    // These methods are used to build structures   ---------------------------------------------------------------------------------

    /*Adds a new sub unit to the world by creating a new graph for it */
    GraphID Add(SubUnit *sub, subName name, string tag = "");
    GraphID Add(string, subName name, string tag = "");

    /* Adds a new sub-unit and links it with an existing sub-unit by a the specified reference point. This grows the graph */
    GraphID Link(SubUnit *sub2, refPoint r2, refPoint r1, string tag = "");
    GraphID Link(string, refPoint r2, refPoint r1, string tag = "");

    /* Wraps a given graph in a structure name. The graph now is a node that we can build with. */
    GraphID Add(GraphID gid, structName name);

    /*Adds a new structure wrapping an existing structure/sub-unit. */
//    GraphID Add(string, string);

    /* Adds the structure r2, by linking it to structure r1, thus growing the graph. */
    GraphID Link(GraphID gid, refPoint r1, refPoint r2);

    // The following methods provides access to querying and setting parameters  -----------------------------------------------------

    // Return a list of all parameters for LAST evaluated expression.
    ParameterList getParams();

    // Returns a copy of all parameters, involving q.
    ParameterList getParamsq();

    // Returns a parameter list containing only q, and potentially set its value.
    ParameterList getq(double =0);
    ParameterList q(double =0);

    // Sets value of string in a parameter list.
    double setParameter(ParameterList& pl, string str, double x);

    // The following methods allows easy evaluation of expressions   -----------------------------------------------------

    // Methods to generate vectors of q values either linear distributed, log distributed, or log distributed with specific base.
    DoubleVector linspace(double q1,double q2,int n);
    DoubleVector logspace(double q1,double q2,int n);
    DoubleVector logspace(double q1,double q2,double base);

    // Evaluate expression for given parameters returning a value.
    double Evaluate(sebsym::Expression, ParameterList&);

    // Evaluate expression for given parameters, at a single q value and returns a double.
    double Evaluate(sebsym::Expression, ParameterList&, double);

    // Evaluate expression for given parameters for a vector of q values, returns a vector of values
    DoubleVector Evaluate(sebsym::Expression, ParameterList&, DoubleVector& );

    // This does the same as the function above, but saves the result to a file.
    // The first optional argument is a user text, the second the character denoting a comment.
    DoubleVector Evaluate(sebsym::Expression e, ParameterList& pl, DoubleVector& q, string, string ="", string = "#");

    // Direct numerical evaluation of analytic and numerical sub-units in mixed structures.
    double EvaluateFormFactor(string name, const ParameterList& values, double q);
    DoubleVector EvaluateFormFactor(string name, const ParameterList& values, const DoubleVector& q);
    double EvaluateFormFactorUnnormalized(string name, const ParameterList& values, double q);
    DoubleVector EvaluateFormFactorUnnormalized(string name, const ParameterList& values, const DoubleVector& q);

    double EvaluateFormFactorAmplitude(refPoint ref, const ParameterList& values, double q);
    DoubleVector EvaluateFormFactorAmplitude(refPoint ref, const ParameterList& values, const DoubleVector& q);
    double EvaluateFormFactorAmplitudeUnnormalized(refPoint ref, const ParameterList& values, double q);
    DoubleVector EvaluateFormFactorAmplitudeUnnormalized(refPoint ref, const ParameterList& values, const DoubleVector& q);

    double EvaluatePhaseFactor(refPoint r1, refPoint r2, const ParameterList& values, double q);
    DoubleVector EvaluatePhaseFactor(refPoint r1, refPoint r2, const ParameterList& values, const DoubleVector& q);

    double EvaluateRadiusOfGyration2(string name, const ParameterList& values);
    double EvaluateSMSDRef2Scat(refPoint ref, const ParameterList& values);
    double EvaluateSMSDRef2Ref(refPoint r1, refPoint r2, const ParameterList& values);

    // These methods are used to produce analytic expressions for scattering tems   --------------------------------------

    // Methods for getting phase factors, normalization is never an issue for phase factors.
    sebsym::Expression PhaseFactor        ( refPoint r1, refPoint r2, int depth = WORLDMAXDEPTH, int varform=QVAR);
    sebsym::Expression PhaseFactorX       ( refPoint r1, refPoint r2, int depth = WORLDMAXDEPTH);
    sebsym::Expression PhaseFactorGeneric ( refPoint r1, refPoint r2, int depth = WORLDMAXDEPTH);

    // Methods for getting F->1 for q->0 normalized form factors and form factor amplitudes
    sebsym::Expression FormFactorAmplitude        ( refPoint ref, int depth = WORLDMAXDEPTH, int varForm = QVAR );
    sebsym::Expression FormFactorAmplitudeX       ( refPoint ref, int depth = WORLDMAXDEPTH);
    sebsym::Expression FormFactorAmplitudeGeneric ( refPoint ref, int depth = WORLDMAXDEPTH);
    sebsym::Expression FormFactor        ( string name, int depth = WORLDMAXDEPTH, int varForm = QVAR );
    sebsym::Expression FormFactorX       ( string name, int depth = WORLDMAXDEPTH);
    sebsym::Expression FormFactorGeneric ( string name, int depth = WORLDMAXDEPTH);

    // Methods for getting unnormalized form factors and - amplitudes.
    // For q->0 they converge to sum beta
    sebsym::Expression FormFactorAmplitude_Unnormalized( refPoint ref, int depth = WORLDMAXDEPTH, int varForm = QVAR );
    sebsym::Expression FormFactorAmplitudeX_Unnormalized       ( refPoint ref, int depth = WORLDMAXDEPTH);
    sebsym::Expression FormFactorAmplitudeGeneric_Unnormalized( refPoint ref, int depth = WORLDMAXDEPTH);

    // For q->0 they converge to (sum beta)^2
    sebsym::Expression FormFactor_Unnormalized        ( string name, int depth = WORLDMAXDEPTH, int varForm = QVAR );
    sebsym::Expression FormFactorX_Unnormalized       ( string name, int depth = WORLDMAXDEPTH);
    sebsym::Expression FormFactorGeneric_Unnormalized( string name, int depth = WORLDMAXDEPTH);

    // Normalization constants
    sebsym::Expression FormFactorAmplitude_Normalization( refPoint ref, int depth = WORLDMAXDEPTH );    // = sum beta
    sebsym::Expression FormFactor_Normalization       ( string name, int depth = WORLDMAXDEPTH);        // = (sum beta)^2

    // These methods are used to provide analytic expressions for Radius of gyration etc. ---------------------------------

    // Measure apparent sizes. Weighed by beta terms and normalized by (sum beta)^2 and (sum beta), respectively.
    sebsym::Expression RadiusOfGyration2( string name, int depth = WORLDMAXDEPTH );
    sebsym::Expression SMSD_ref2scat( refPoint ref, int depth = WORLDMAXDEPTH );

    // Measure absolute mean-square distance between reference point pairs.
    sebsym::Expression SMSD_ref2ref( refPoint r1, refPoint r2, int depth = WORLDMAXDEPTH);

    // Miscellanious methods. --------------------------------------------------------------------------------------------

    // Count number of structures/sub-units connected to ref at the specifie depth
    sebsym::Expression Count(refPoint ref, int depth = WORLDMAXDEPTH);

    // Count number of pairs of structures/sub-units within structure at specific depth. Count / Countpairs returns N vs. N^2 elements.
    sebsym::Expression CountPairs( string name, int depth = WORLDMAXDEPTH);

    // Find a path of reference points connecting the two reference points, recursing fown to the specified level.  (NOT used by SEB)
    ReferencePointList  Path( refPoint r1, refPoint r2, int depth=WORLDMAXDEPTH, bool =true);

    // Expose symbols to the user.
    SymbolInterface* GetSymbolInterface() { return GLEX; }

    // finds horizontal path between two reference points / structures / subunits.          Default check userspecified arguments.
    ReferencePointList findpath(string name1, string name2, bool =true);

  private:
    struct GuinierTerm {
        double constant;
        double q2Coefficient;

        GuinierTerm(double c0 = 0.0, double c2 = 0.0)
            : constant(c0), q2Coefficient(c2) {}
    };

    // These three methods does all the traversal used in all scattering expressions above.
    sebsym::Expression GenerateRefToRef( refPoint r1, refPoint r2, int depth, int varForm );
    sebsym::Expression GenerateRefToAll( refPoint r,               int depth, int varForm );
    sebsym::Expression GenerateAllToAll( string name,              int depth, int varForm );

    // Parallel numerical traversal.
    double GenerateTotalBetaNumeric(string name, const ParameterList& values);
    double GenerateRefToRefNumeric(refPoint r1, refPoint r2, double q, const ParameterList& values);
    double GenerateRefToAllNumeric(refPoint ref, double q, const ParameterList& values);
    double GenerateAllToAllNumeric(string name, double q, const ParameterList& values);
    double PhaseFactorNumeric(ReferencePointList& path, double q, const ParameterList& values);

    // Guinier traversal using c0 - q^2 c2 + O(q^4).
    GuinierTerm GenerateRefToRefGuinier(refPoint r1, refPoint r2, const ParameterList& values);
    GuinierTerm GenerateRefToAllGuinier(refPoint ref, const ParameterList& values);
    GuinierTerm GenerateAllToAllGuinier(string name, const ParameterList& values);
    GuinierTerm PhaseFactorGuinier(ReferencePointList& path, const ParameterList& values);
    static GuinierTerm AddGuinier(const GuinierTerm& left, const GuinierTerm& right);
    static GuinierTerm MultiplyGuinier(const GuinierTerm& left, const GuinierTerm& right);
    static GuinierTerm ScaleGuinier(double scale, const GuinierTerm& term);

    // Makes a list of all neighbors, that is link partners, and reference points inside the same structure / sub-unit
    ReferencePointList getNeighbors( refPoint last, ReferencePointList& VisitedAlready);

    // Returns a list of all reference points contained in structure / subunit.
    ReferencePointList getReferencePoints( refPoint structure );

    // Generate a list of all reference points inside a structure / sub-unit
    void addStructureReferences( refPoint r, ReferencePointList& ret );

    // is r in VisitedAlready?
    bool isVisited( refPoint& r, ReferencePointList& VisitedAlready);

    // Generate a link between the two reference points
    LinkPair generateLink(refPoint r1, refPoint r2);

    // Helper to generate product of phase factors along a given path.
    sebsym::Expression PhaseFactor(ReferencePointList& path, int depth, string, int varForm, bool doCheck=true);

    // Const iterators over subgraphs
    list<string>::const_iterator subgraph_cbegin(GraphID);
    list<string>::const_iterator subgraph_cend(GraphID);

    // Calculates depth of a structure   depth = number of : in string
    int  Depth(string);
    int  Depth(string,string);

    /* checks whether strings contains certain symbols. */
    bool hasAHash(refPoint r);
    bool hasAColon(refPoint r);
    bool hasAPeriod(refPoint r);

    // Useful string operations for paths in structures.
    string prefix(string name);      // Returns the sub string that is before the first occurence of ":", or before ".", or everything. Throws if string starts with : or .
    string postfix(string name);     // Returns the sub string that is after the first occurance of ":" . Throws if no :, or : is the last character.
    string getPath(string);          // Returns the sub string from start to .  Throws if . not in string, or if string starts with .

    string getName(string);          // Returns the sub-unitname, that is from start or right most :, until end or .    Throws if :. or .xxx: found.
    string getReference(string);     // Returns "ref#pr" from "structure..:name.ref#pt". That is everything after . Throws if . not in string.

    string getReferenceBase(string);       // Returns "refbase" from "subunitname.redbase#mypoint". That is everything between . or #/end
    string getReferenceBaseHash(string);   // Identical to getReferenceBase except it throws unless noth . and # are in string.
    string getReferenceAfterHash(string);  // Returns everything after #,   throws if no # or # last character

    // test syntax of a path to a reference point if valid.
    void testPathSyntax(string, bool =false);

    // tests if a supplied pointer points to a sub-unit.
    bool testSubunitPointer(SubUnit*);

    // tests if a supplied graphID is valid.
    bool testGraphID(int);

    // builds a human-readable list of valid reference points for a sub-unit, e.g. "end1, end2, contour (distributed)"
    string subunitReferenceList(SubUnit* s);


    // Code for generating a symbolic Psi depending on various varforms
    sebsym::Expression getPsi(string myself, refPoint r1, refPoint r2, int varform);

    // Code for generating a symbolic Form Factor Amplitude depending on various varforms
    sebsym::Expression getFFA(string myself, refPoint r, int varform);

    // Code for generating a symbolic Form Factor depending on various varforms
    sebsym::Expression getFF(string myself, int varform);
};


#endif // INCLUDE_GUARD_WORLD
