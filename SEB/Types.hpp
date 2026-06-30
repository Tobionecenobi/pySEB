#include <string>
#include <set>
#include <vector>
#include <utility>
#include <list>
#include <map>
#include "Expression.hpp"
using namespace std;

//===========================================================================
// type definintions

typedef string refPoint;
typedef string distRefPoint;

typedef int id;
typedef int GraphID;

typedef string structName;
typedef string subName;
typedef pair<refPoint, refPoint> LinkPair;

typedef set<refPoint> ReferencePointSet;     // Set of reference points, fast for querying.
typedef list<refPoint> ReferencePointList;   // Paths are lists of reference point strings
typedef std::map<std::string, double> ParameterList;  // Lists of parameters with their associated values.
typedef sebsym::Expression ex;                          // Symbolic expression type
typedef std::map<sebsym::Expression, sebsym::Expression> ExpressionMap;  // Map of expressions to expressions for substitution

typedef vector<double> DoubleVector;         // Vectors of values.
