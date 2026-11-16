// Graph queries, path searches, and reference-name handling.
#include "World.hpp"

#include <algorithm>
#include <set>
#include <utility>
#include <vector>


/*  Find a path inside a given structure.   (here using DiBlockStarChain.cpp as an example)

    Example:  chain:star1:diblock1:polyB.end2 and    chain:star5:diblock3:polyB.end2

    The path is between two fully specified sub-units and potentially reference points hence the lengthy expressions.
    That is a path from a common structural parent (chain) down to the specific sub-unit (and reference point).

    These searches are always horizontal in a structure in the sense that a path is defined on the immediate
    children of the parent structure. Here star1-star5 are the immediate children of chain, thus the path is
    defined as the star level, and never involves steps at the diblock level.

    Paths always exist and are unique because by construction. They exist because we require the
    searches to be rooted at the same level in the structure, and everything is connected at every
    level. Secondly, paths are unique because with addAndLink we grow graphs and structures by
    adding leafs to an existing graph, thus we can never form a loop by construction.


      CASE 1: Path between two reference points:

      Example:
         name1 is star1:diblock1:polyB.end2
         name2 is star1:diblock3:polyB.end2

      The result is diblock1:polyB.end2 -> diblock1:polyA.end1 <-> diblock3:polyA.end1 -> diblock3:polyB.end2


      Example:
        name1 is chain:star1:diblock1:polyB.end2
        name2 is chain:star5:diblock3:polyB.end2

      The result is star1:diblock1:polyB.end2 -> star1:diblock3:polyB.end2 <-> star2:diblock1:polyB.end2 -> star2:diblock3:polyB.end2 <->
                    star3:diblock1:polyB.end2 -> star3:diblock3:polyB.end2 <-> star4:diblock1:polyB.end2 -> star4:diblock3:polyB.end2 <->
                    star5:diblock1:polyB.end2 -> star5:diblock3:polyB.end2

      Paths are always alternating sequences of jumps across sub-units/structures (->) connecting two reference points on the same
      sub-unit/structure and jumps across links (<->).

      Example: The lowest level we can search is inside a structure that directly resolves to sub-units
        name1 is  diblock1:polyA.end1
        name1 is  diblock1:polyB.end2

      The result is   polyA.end1 -> polyA.end2 <-> polyB.end1 -> polyB.end2



      CASE 2: path between a reference point and a structure / sub-unit. Here we flip so name1 is always the refence point.

        name1 is chain:star1:diblock1:polyB.end2
        name2 is chain:star5:diblock3:polyB

       The result is star1:diblock1:polyB.end2 -> star1:diblock3:polyB.end2 <-> star2:diblock1:polyB.end2 -> star2:diblock3:polyB.end2 <->
                     star3:diblock1:polyB.end2 -> star3:diblock3:polyB.end2 <-> star4:diblock1:polyB.end2 -> star4:diblock3:polyB.end2 <->
                     star5:diblock1:polyB.end2

       The search terminates because chain:star5:diblock3:polyB is found inside star5. So we have found the structure
       containing the target. Because searches are horizontal we do not search at the diblock level towards the target.

       But to find out how to get from star5:diblock1:polyB.end2 to chain:star5:diblock3:polyB we would have again search
       for that path, this is a search at the diblock level one below star5. The result is
                   diblock1:polyB.end2 -> diblock1:polyA.end1 <-> diblock3:polyA.end1

       Again the search terminates at diblock3:polyA.end1 because diblock3:polyB is found inside diblock3. To find the
       path from the reference point to the unit, we have to search diblock3:polyA.end1 to diblock3:polyB

       The result is polyA.end1 -> polyA.end2 <-> polyB.end1  and we have found which reference point
       chain:star5:diblock3:polyB.end1 which is closest to chain:star1:diblock1:polyB.end2



      CASE 3: path between a pair of structures / sub-unit2.

         name1 is chain:star1:diblock1:polyB
         name2 is chain:star5:diblock3:polyB

         In this case we seed the search with start1:*.refs and then search outwards to other structures, and
         terminate when we hit a star5:* reference point. Hence the path above is identical to the search.

         name1 is chain:star1
         name2 is chain:star5


         The path found is

                  star1:diblock3:polyB.end2 <-> star2:diblock1:polyB.end2 -> star2:diblock3:polyB.end2
                                            <-> star3:diblock1:polyB.end2 -> star3:diblock3:polyB.end2
                                            <-> star4:diblock1:polyB.end2 -> star4:diblock3:polyB.end2
                                            <-> star5:diblock1:polyB.end2

         The paths starts at star1:diblock3:polyB.end2   since star1:diblock1:polyB is inside star1.
         The path ends at    star5:diblock1:polyB.end2   since star5:diblock3:polyB is inside star5

         Thus the corresponding scattering expression would be

           Form factor amplitude of star1:diblock1:polyB relative to star1:diblock3:polyB.end2
         x Phase factor of star2:diblock1:polyB.end2 -> star2:diblock3:polyB.end2
         x Phase factor of star3:diblock1:polyB.end2 -> star3:diblock3:polyB.end2
         x Phase factor of star4:diblock1:polyB.end2 -> star4:diblock3:polyB.end2
         x Form factor amplitude of star5:diblock3:polyB relative to star5:diblock1:polyB.end2

         Additional path searches between these are required to resolve the structural
         scattering terms in terms or their substructures / subunits.

      Example:

          name1 is star1:diblock1:polyB
          name2 is star1:diblock3:polyB

          The path found is diblock1:polyA.end1 <-> diblock3:polyA.end1

          The path starts at diblock1:polyA.end1  since diblock1:polyB is inside diblock1
          The path ends at   diblock3:polyA.end1  since diblock3:polyB is inside diblock3


      Example:

          name1 is diblock1:polyB
          name2 is diblock1:polyA.end1

          The path is polyB.end1 <-> polyA.end2 -> polyA.end1


      Example:

          name1 is diblock3:polyA.end1
          name2 is diblock3:polyB

          The path is polyA.end1 -> polyA.end2 <-> polyB.end1



      The breadth-first search queues references with predecessor indices and
      reconstructs the successful path. Neighbors retain their original order.

      CASE1: We seed the search with the starting reference point and terminate when we have found the termination point.

      CASE2: We seed the search with the starting/termination reference point, and terminate when we have found a
             structure containing the target/starting sub-unit.

      CASE3: We seed the search with all reference points within a structure, and terminate when we have found a
             structure containing the target sub-unit.

*/

ReferencePointList World::findpath(string name1, string name2, bool check)
{
    return findPath(std::move(name1), std::move(name2), check, nullptr);
}

ReferencePointList World::findPath(string name1, string name2, bool check,
                                  NumericEvaluationContext* context)
{
try{

    if(prefix(name1) != prefix(name2))
        throw SEBException("Both search paths should start in the same structure. name1='" + name1 + "' (prefix='" + prefix(name1) + "'), name2='" + name2 + "' (prefix='" + prefix(name2) + "').");
    if (check) testPathSyntax(name1);
    if (check) testPathSyntax(name2);

    name1=postfix(name1);
    name2=postfix(name2);

    if (prefix(name1)==prefix(name2))       // Returns an empty path.
       {
             ReferencePointList Path;
             return Path;
       }

    bool reversesearch=false;
    if (!hasAPeriod(name1) && hasAPeriod(name2))
         {
            swap(name1,name2);                                 // Optimization: Search from referece point to sub-unit, rather than from
            reversesearch=true;                                // sub-unit to reference point, since there are many potential starting points.
         }

    // A generation asks for many targets from the same source. Reuse the
    // ordered BFS predecessor tree, not a separate search for every pair.
    // Scalar/public path queries retain the original early-exit search below.
    if (context) {
        auto cached = context->searches.find(name1);
        if (cached == context->searches.end() && context->searches.size() < 64 &&
            context->referenceEntries < 262144) {
            cached = context->searches.emplace(name1, PathSearch{}).first;
            ++context->referenceEntries;
        }
        // Full searches pay off only for repeated sources. Leave one-off
        // phase queries and small blocks on the early-exit algorithm.
        if (cached != context->searches.end() && !cached->second.attempted &&
            ++cached->second.uses >= 3) {
            cached->second.attempted = true;
            const std::size_t noParent = static_cast<std::size_t>(-1);
            std::vector<PathSearchNode> tree;
            std::unordered_set<refPoint> visited;
            bool complete = true;
            const auto append = [&](const refPoint& reference, std::size_t parent) {
                if (tree.size() >= 4096) return false;
                tree.push_back({reference, prefix(reference), parent});
                visited.insert(reference);
                return true;
            };
            try {
                if (hasAPeriod(name1)) complete = append(name1, noParent);
                else for (const auto& ref : getReferencePoints(prefix(name1)))
                    if (!append(ref, noParent)) { complete = false; break; }
                for (std::size_t head = 0; complete && head < tree.size(); ++head) {
                    // Snapshot filtering and insertion order match the
                    // original BFS, including any duplicate neighbors.
                    const auto neighbors = getNeighbors(tree[head].reference, visited, context);
                    for (const auto& next : neighbors)
                        if (!append(next, head)) { complete = false; break; }
                }
            } catch (const SEBException&) {
                // An optional full search must not introduce an error beyond
                // where the original target-specific search would terminate.
                complete = false;
            }
            const std::size_t entries = 2*tree.size();
            if (complete && entries <= 262144-context->referenceEntries) {
                context->referenceEntries += entries;
                cached->second.nodes = std::move(tree);
            }
        }
        if (cached != context->searches.end() && !cached->second.nodes.empty()) {
            const bool targetIsReference = hasAPeriod(name2);
            const auto targetObject = prefix(name2);
            const auto& tree = cached->second.nodes;
            const std::size_t noParent = static_cast<std::size_t>(-1);
            for (std::size_t i = 0; i < tree.size(); ++i) {
                if (tree[i].parent != noParent &&
                    (targetIsReference ? tree[i].reference == name2 : tree[i].object == targetObject)) {
                    ReferencePointList path;
                    for (std::size_t node = i; node != noParent; node = tree[node].parent)
                        path.push_front(tree[node].reference);
                    if (reversesearch) path.reverse();
                    return path;
                }
            }
            // Missing targets use the original search/error reporting.
        }
    }

    // Store one predecessor per queued reference, rather than copying a
    // complete candidate path whenever BFS visits or branches to a neighbor.
    struct SearchNode { refPoint reference; std::size_t parent; };
    const std::size_t noParent = static_cast<std::size_t>(-1);
    std::vector<SearchNode> queue;
    std::unordered_set<refPoint> visited;
    const bool targetIsReference = hasAPeriod(name2);
    const auto seed = [&](const refPoint& reference) {
        queue.push_back({reference, noParent});
        visited.insert(reference);
    };
    if (hasAPeriod(name1)) seed(name1);
    else for (const auto& reference : getReferencePoints(prefix(name1))) seed(reference);

    // Preserve seed and neighbor order, including the existing reverse-search
    // convention. Reconstruct only the successful path from predecessor indices.
    for (std::size_t head=0; head<queue.size(); ++head) {
        const ReferencePointList neighbors = getNeighbors(queue[head].reference, visited, context);
        for (const auto& next : neighbors) {
            queue.push_back({next, head});
            if ((targetIsReference && next == name2) ||
                (!targetIsReference && prefix(next) == prefix(name2))) {
                ReferencePointList path;
                for (std::size_t node=queue.size()-1; node!=noParent; node=queue[node].parent)
                    path.push_front(queue[node].reference);
                if (reversesearch) path.reverse();
                return path;
            }
            visited.insert(next);
        }
    }

    throw SEBException("Didn't find a path between "+name1+"  "+name2+"    This should not happen!");
}
catch (SEBException& e)
{
   e.PushCallStack("ReferencePointList World::findpath(string name1=\""+name1+"\", string name2=\""+name2+"\")" );
   throw;
}
}


/*  Not used by SEB, but provides paths between pairs of reference points at a user specified depth. Thus path is specified at the
    same height as the two starting reference points.  */
ReferencePointList  World::Path( refPoint r1, refPoint r2, int depth, bool doCheck)
{
try{
    string myself=prefix(r1);
    if(myself != prefix(r2))   throw SEBException("Both reference points should be in the same structure / reference point.");
    if (depth<0)               throw SEBException("Depth can not be negative.");
    if (!hasName(myself))      throw SEBException(string("Unknown name ")+myself+" in world.");


    // If we have two paths to a reference point inside a structure:
    if (isStructure(myself))
      {
        // These checks are expensive, and disabled in subsequent recursions.
        if(doCheck) testPathSyntax(r1);
        if(doCheck) testPathSyntax(r2);

        // Trivial cases
        ReferencePointList rpl;
        if ( isLinked(r1, r2) || r1 == r2 || depth==0)
            {
              rpl.push_back(r1);
              rpl.push_back(r2);
              return rpl;
            }

        // depth>0   hence we could have XX:YY... XX:ZZ...  or  XX:YY... XX:YY...  in the second case the yy-yo-yy path would be empty, so we have to handle thus case
        if (prefix(postfix(r1)) == prefix(postfix(r2)))   // XX:YY.. and XX:yy..
           {
               rpl = Path( postfix(r1) , postfix(r2), depth - 1, false);
           }
         else
           {                                                                                    // The two points have different second prefix, hence the path will not be empty.
               ReferencePointList path = findpath(r1, r2, false);                              // Find path, don't check arguments.

               auto it=path.begin();
               auto its=it; its++;
               while (its != path.end())
                {
                    if(!isLinked(*it, *its))                                                                                 // steps across links contribute no phase
                             rpl.splice( rpl.end(), Path(*it, *its, depth-1, false) );

                    it=its; its++;
                }
           }

        for (auto& p : rpl) p=myself+":"+p;   // Prepent myself, since findpath removes the common root structure.
        return rpl;
    }
    else
     if (isSubunit(myself))     // r1,r2 must have the form   r1=subunit.ref1  r2=subunitname.ref2
          {
              ReferencePointList rpl { r1, r2};
              return rpl;
         }

    throw SEBException("Internal error.");

}
catch (SEBException& e)
{
   e.PushCallStack("ReferencePointList  World::Path(refPoint "+r1+", refPoint "+r2+", int "+to_string(depth)+"..)");
   throw;
}
}

/* Is specified name in the nameCatalog? */
bool World::hasName(string n)
{
   return nameCatalog.find(n) != nameCatalog.end();
}

/* Is specified name in the specified subgraph? */
bool World::doesSubgraphContainName(int gid, string n)
{
   if (!testGraphID(gid))    throw SEBException("Bad graphid "+to_string(gid), "World::doesSubgraphContainName(GraphID ="+to_string(gid)+", string=\""+n+"\")");

   for ( auto item : subGraphs.at(gid))
      if (item == n) return true;
   return false;
}


/* Checks that parent structure contains subling (structure / subunit name). Only checks one level down. */
bool World::doesStructureContainName(string parent, string sibling)
{
   if (!isStructure(parent)) throw SEBException("String does not refer to a structure.", "World::doesStructureContainName(string parent=\""+parent+"\", string sibling=\""+sibling+"\")");

   GraphID gid=getStructure(parent)->getGraphID();
   return doesSubgraphContainName(gid, sibling);
}


/* Returns a GraphID from a unique name of a sub unit or a subgraph

   Should be e.g. structure_name or subunit_name, will fail if given structure:sub-unit, since then its unclear what graphID referred to, that of the structure, of the sub-unit.
*/
GraphID World::getGraphID(string name)
{
    if (hasAColon(name))  throw SEBException("The specified name contains a :", "World::getGraphID(string "+name+")");

    for (auto itr : subGraphs)
        for (auto itr2 : itr.second)
            if (itr2 == name) return itr.first;

    throw SEBException("The specified name was not found in any sub-graph", "World::getGraphID(string "+name+")");
}

/*Returns the graphID id given a pointer to a sub-unit or structure */
GraphID World::getGraphID(SubUnit *s)    {  return getGraphID( s->getName() ); }
GraphID World::getGraphID(Structure *s)  {  return getGraphID( s->getName() ); }

void World::SetGraphLabel(GraphID gid, const string& label)
{
    if (!testGraphID(gid)) throw SEBException("Bad graphid:" + to_string(gid), "World::SetGraphLabel()");
    if (label.empty()) throw SEBException("Graph label must not be empty", "World::SetGraphLabel()");
    for (const auto& item : graphLabels) {
        if (item.first != gid && item.second == label) {
            throw SEBException("Graph label " + label + " is already in use", "World::SetGraphLabel()");
        }
    }
    graphLabels[gid] = label;
}

// return sub-unit pointer to specified sub-unit. Throws if name unknown, or pointer is not to a sub-unit
SubUnit* World::getSubunit(string n)
{
   if (!hasName(n)) throw SEBException("The specified subunit is unknown", "World::getSubunit(string "+n+")");
   SubUnit* sptr = reinterpret_cast<SubUnit*>( nameCatalog.find(n)->second );

   if (sptr->getType() != SUBUNIT && sptr->getType() != SUBUNITCHILD)
       throw SEBException("The specified name does not point to a sub-unit!", "World::getSubunit(string "+n+")");

   return sptr;
}

// return structure pointer to specified structure.  Throws if name unknown, or pointer is not to a structue
Structure* World::getStructure(string n)
{
   if (!hasName(n)) throw SEBException("The specified structure is unknown", "World::getStructure(string "+n+")");
   Structure* sptr = reinterpret_cast<Structure*>( nameCatalog.find(n)->second );

   if (sptr->getType() != STRUCTURE)
       throw SEBException("The specified name does not point to a structure!", "World::getStructure(string "+n+")");

   return sptr;
}

// Helper returning whether "name" is a sub-unit. Throws if name does not exist.
bool World::isSubunit(string n)
{
   if (!hasName(n)) throw SEBException("The specified name is unknown", "World::isSubunit(string "+n+")");

   ABSSubUnit* aptr=nameCatalog.find(n)->second;
   int typ=aptr->getType();
   return  typ== SUBUNIT || typ == SUBUNITCHILD;
}

// Helper returning whether "name" is a structure. Throws if name does not exist.
bool World::isStructure(string n)
{
   if (!hasName(n)) throw SEBException("The specified name is unknown", "World::isStructure(string "+n+")");

   ABSSubUnit* aptr=nameCatalog.find(n)->second;
   return aptr->getType() == STRUCTURE;
}


/* Checks if to reference point are linked */
bool World::isLinked(refPoint r1, refPoint r2)
{
    if (r2<r1) swap(r1,r2);   // We store links alphabetically sorted.
    for(const auto& link : links)
        if(link.first == r1 && link.second == r2)  return true;

    return false;
}


// Helper function that gets all reference points inside a structure, that is all paths
// connecting structure to a reference point.
ReferencePointList World::getReferencePoints( refPoint structure)
{
// test structure string is good.

    ReferencePointList references;
    addStructureReferences(structure, references );
    return references;
}

/* Sorts the link in ascii order so the smallest ascii number comes first in the link where a link = pair<refpoint, refpoint>
   The sanity of the strings should have been tested before here*/
LinkPair World::generateLink(refPoint r1, refPoint r2)
{
    if (r1 == r2) throw SEBException("A reference point cannot be linked with itself", "World::generateLink(refPoint r1=" + r1 + ", refPoint r2=" + r2 + ")");
    if (r2 < r1) swap(r1, r2);   // Store alphabetically so isLinked can find them.
    return make_pair(r1, r2);
}



// Const iterator over subgraphs
list<string>::const_iterator World::subgraph_cbegin(GraphID id)
{
   return (subGraphs.find(id)->second).cbegin();
}

// const iterator over subgraphs termination point.
list<string>::const_iterator World::subgraph_cend(GraphID id)
{
   return (subGraphs.find(id)->second).cend();
}



/* Checks for various characters */
bool World::hasAHash(refPoint r)  {  return r.find("#") != string::npos; }
bool World::hasAColon(refPoint r) {  return r.find(":") != string::npos; }
bool World::hasAPeriod(refPoint r){  return r.find(".") != string::npos; }


/* Checks if a sub-unit pointer is good. This could still coredump if pointer points to a bad place in memory. */
bool World::testSubunitPointer(SubUnit* s)
{
   if (s==nullptr) return false;
   if (s->getType() != SUBUNIT && s->getType() != SUBUNITCHILD ) return false;

   return true;
}



/* Recursively checks a path through a structure.

      structure1:structure2:...:subunitname [.refname]

      Checks that all structure names are known,
             that the subunit name is known.
             that the sub-unit has the basename of the reference point (if defined).
             and that structure1 has a structure2 in it and so on.

    In case we are bulding a new structure, the prefix structure does not exist yet.   (skipfirst=true in this case)
*/

void World::testPathSyntax(string r, bool skipfirst)  // default false, that is check everything
{
    if (hasAColon(r))
           {           // CASE: structure1:X

               if (!skipfirst)                                       // If skipfirst is true then we skip testing structure1, since it hasn't been created yet.
                   {
                      if (!hasName(prefix(r)))
                               throw SEBException("Structure "+prefix(r)+" is not known");

                      if (!doesStructureContainName(prefix(r), prefix(postfix(r))))
                               throw SEBException("Structure "+ prefix(r)+" does not contain a "+prefix(postfix(r))+" inside it");
                   }

               return testPathSyntax(postfix(r));                    // Its known, so continue to test X.

           }
      else
           {              // r=subunit.ref
               string name=getName(r);
               if (!hasName(name))
                      throw SEBException("Subunit "+name+" is not known");

               // Does the sub-unit recognize the base reference point?
               if (hasAPeriod(r) && !getSubunit(name)->hasReference( getReferenceBase(r) ))
               {
                      SubUnit* s = getSubunit(name);
                      throw SEBException("Subunit '"+name+"' does not have reference point '"+getReferenceBase(r)+"'. Valid reference points: "+subunitReferenceList(s));
               }
           }
}


bool World::testGraphID(int g)
{
   return g>=1 && g<=totalNumberofGraphs;
}

// Builds a human-readable list of all valid reference points for a sub-unit.
string World::subunitReferenceList(SubUnit* s)
{
    string result;
    auto append = [&](const string& name) {
        if (!result.empty()) result += ", ";
        result += name;
    };
    for (auto it = s->SpecificReferencepoints_cbegin(); it != s->SpecificReferencepoints_cend(); ++it)
        append(*it);
    for (auto it = s->DistributedReferencepoints_cbegin(); it != s->DistributedReferencepoints_cend(); ++it)
        append(*it + " (distributed)");
    if (result.empty()) result = "(none)";
    return result;
}


/*  Returns the prefix before ":" or if not found before "." or the entire string.

    If called with
        structure1:structure2:subunit.ref     returns structure1
        subunit.ref                           returns subunit
        subunit                               returns subunit

    Throws if called with strings :xxx  or .xxx that has an empty prefix.

*/
string World::prefix(string n)
{
    int start=n.find(":");

    if (start==0) throw SEBException("String of form :xx has empty prefix", "World::prefix(\""+n+"\")");
    if (start != string::npos)     return n.substr(0, start);

    start=n.find(".");
    if (start==0) throw SEBException("String of form .xx has empty prefix", "World::prefix(\""+n+"\")");
    if (start != string::npos)     return n.substr(0, start);

    return n;
}

/*
Need a helper to count to occurances of a character to make the logic below more sharp.
*/

// Returns everything what is after the first :
// Throws if : does not exist, or : is at the end of the string.
string World::postfix(string n)
{
    int from=n.find(":");
    if (from == string::npos)     throw SEBException(": not found in string", "World::postfix(\""+n+"\")");
    from++;
    if (from == n.length())       throw SEBException("Strings of form xx: has empty postfix", "World::postfix(\""+n+"\")");

    return n.substr(from);      // return the postfix
}


// Returns the path part of a string like "structure:structure:subunit.ref" by discarding the reference part.
//  Throws if no .
string World::getPath(string n)
{
    int end=n.find(".");
    if (end == string::npos)     throw SEBException(". not found in string", "World::getPath(\""+n+"\")");
    if (end == 0)                throw SEBException("String of form .xx has no path", "World::getPath(\""+n+"\")");

    return n.substr(0, end);   // return the path
}

//  Returns the subname from a string of the form  [structures:]subunitname[.referencepoint.]  Thus everything from start or rigtht most : to . or end.
string World::getName(string n)
{
    int start=0;
    if (n.rfind(":") != string::npos)  start=n.rfind(":")+1;

    int end=n.length();
    if (n.find(".") != string::npos)  end=n.find(".");

    if (end<start)          throw SEBException("String of form .xxx: has no name", "World::getName(\""+n+"\")");
    if (start == end)       throw SEBException("String of form :. has no name", "World::getName(\""+n+"\")");

    return n.substr(start, end-start);   // return name part
}


/*  Returns the reference point part from a string, that is everything after .
    Hence referencepoint#mypoint from e.g. structurenamae:..;subunitname.referencepoint[#mypoint]
    Throws if the string does not contain .
*/
string World::getReference(string n)
{
    int start=n.rfind(".");
    if (start == -1)               start=n.find(".");
    if (start == string::npos)     throw SEBException("String without . has no reference found", "World::getReference(\""+n+"\")");
    if (start == n.length()-1)     throw SEBException("String of form xx. has no reference found", "World::getReference(\""+n+"\")");

    return n.substr(start+1);   // return reference part
}


//  Returns the refbase from a string of the form  [structures:]subunitname.refbase#mypoint.
//   Thus everything after . until #/end. Throws if no . found or if # is before .
string World::getReferenceBase(string n)
{
    int start=n.find(".");
    if (start == string::npos)     throw SEBException("String without .", "World::getReferenceBase(\""+n+"\")");
    start++;

    int end=n.find("#");
    if (end == string::npos)       return  n.substr(start);                  // CASE: .refbase
    if (start==end)                throw SEBException("String of form .# has no refpoint", "World::getReferenceBase(\""+n+"\")");
    if (end<start)                 throw SEBException("String of form #xx. has no refpoint", "World::getReferenceBase(\""+n+"\")");

    return  n.substr(start, end-start);       // CASE: .refbase#something
}


/*  Returns the refbase from a distributed reference point. The return value is identical to the getReferenceBase above, except that this
    method throws an exception unless both . and # are found.
*/
string World::getReferenceBaseHash(string n)
{
    int start=n.find(".");
    if (start == string::npos)     throw SEBException("String without .", "World::getReferenceBaseHash(\""+n+"\")");
    start++;

    int end=n.find("#");
    if (end == string::npos)       throw SEBException("String without #",  "World::getReferenceBaseHash(\""+n+"\")");
    if (end<start)                 throw SEBException("String of form #xx.", "World::getReferenceBaseHash(\""+n+"\")");
    if (start==end)                throw SEBException("String of form .# has no refpoint", "World::getReferenceBaseHash(\""+n+"\")");

    return n.substr(start, end-start);   // return reference type partm
}


/*  Returns everything from # to the end of the string.
    Throws if # or not found, or at the end of the string.
*/
string World::getReferenceAfterHash(string n)
{
    int start=n.find("#");
    if (start == string::npos)       throw SEBException("String without #",  "World::getReferenceAfterHash(\""+n+"\")");
    if (start+1 == n.length())       throw SEBException("String of form xx#", "World::getReferenceAfterHash(\""+n+"\")");

    return n.substr(start+1);   // return reference type partm
}



/*
   Returns all the non-visited neighbours of the last reference point, that is

      Last is structure1:structure2:....:subunitname.reference

      then potentially we can jump through a link to another structure,
      or we can jump to any reference point inside this structure.
*/

ReferencePointList World::getNeighbors(const refPoint& last,
                                       const std::unordered_set<refPoint>& visited,
                                       NumericEvaluationContext* context)
{
    ReferencePointList uncached;
    const ReferencePointList* candidates = nullptr;
    if (context) {
        const auto found = context->neighbors.find(last);
        if (found != context->neighbors.end()) candidates = &found->second;
    }
    if (!candidates) {
        // Adjacency depends only on topology. Filtering by this search's visited
        // references happens afterwards, so the list can serve other pair searches.
        for (const auto& link : links) {
            refPoint next;
            if (link.first == last) next = link.second;
            if (link.second == last) next = link.first;
            if (!next.empty()) uncached.push_back(next);
        }
        for (const auto& next : getReferencePoints(prefix(last)))
            if (next != last) uncached.push_back(next);
        if (context && context->neighbors.size() < 4096 &&
            uncached.size() <= 262144 - context->referenceEntries) {
            context->referenceEntries += uncached.size();
            candidates = &context->neighbors.emplace(last, std::move(uncached)).first->second;
        } else candidates = &uncached;
    }
    ReferencePointList neighbors;
    for (const auto& next : *candidates)
        if (visited.count(next) == 0) neighbors.push_back(next);
    return neighbors;
}

// Recursive funktion that does the search for structures / sub-units into a structure.
void World::addStructureReferences( refPoint r, ReferencePointList& ret )
{
    refPoint name=getName(r);                                                          // get right most structure/sub-unit in path.
    if (!hasName(name)) throw SEBException("Unknown name "+name,"addStructureReferences(\""+r+"\", ..)");

    if (isSubunit(name)) getSubunit(name)->addSubunitReferences(r, ret);
    else
    if (isStructure(name))
            {
                GraphID gid=getStructure(name)->getGraphID();                           // get gid inside "name".
                for (auto it = subgraph_cbegin(gid); it != subgraph_cend(gid); it++)    // traverse that subgraph
                     addStructureReferences( r+":"+*it , ret);                          // and recursive call function to add below.
            }
    else throw SEBException("Should not happen!","addStructureReferences()");
}
