// World construction, model registration, and structural descriptions.
#include "World.hpp"
#include "Subunits/CreateSubunit.hpp"
#include "Subunits/FileDefined.hpp"

#include <algorithm>
#include <memory>
#include <set>


/* Adds a new sub unit that is not connected to anything. Thus it receives its own graphID.

   Arguments:
           sub a pointer to an uninitialized sub-unit object class, Init initializes it.
           name  a string for naming the sub-unit.
           tag   a string for grouping identical sub-units together.

   Returns graphID.
*/

GraphID World::Add(SubUnit *sub, subName name, string tag)
{
try{
    if (!GLEX->testnamestring(name))                 throw SEBException("Bad symbol in subunit name:"+name+". "+GLEX->invalid_name_message(name));
    if (!tag.empty() && !GLEX->testnamestring(tag))  throw SEBException("Bad symbol in subunit tag:"+tag+". "+GLEX->invalid_name_message(tag));
    if (hasName(name))                               throw SEBException("Name "+name+" already exists in the world");
    if (!testSubunitPointer(sub))                    throw SEBException("Supplied subUnit pointer does not point a sub-unit.");

    sub->Init(name, tag, GLEX);

    totalNumberofGraphs++;                // Every add creates a new graph since the sub-unit is not connected to anything.
    list<subName> subList;
    subList.push_back(name);
    subGraphs.insert(pair<GraphID, list<subName>>(totalNumberofGraphs, subList));
    nameCatalog.insert(pair<subName, SubUnit *>(name, sub));
    typeCatalog.insert(pair<subName, int>(name, sub->getType()));

    return totalNumberofGraphs;
}
catch (SEBException& e)
{
   e.PushCallStack("GraphID World::Add(SubUnit*, subName "+name+", string "+tag+")" );
   throw;
}
}

// Helper
GraphID World::Add(string subtype, subName name, string tag)
{
   std::unique_ptr<SubUnit> subunit(subunitRegistry.Create(subtype));
   GraphID result = Add(subunit.get(), name, tag);
   string canonical = subtype;
   for (const auto& model : subunitRegistry.List()) {
       if (model.id == subtype || model.apiName == subtype) {
           canonical = model.id;
           break;
       }
   }
   modelCatalog[name] = canonical;
   subunit.release();
   return result;
}

GraphID World::AddFile(string path, subName name, string tag)
{
   pyseb::SubunitDefinition definition = pyseb::LoadSubunitDefinitionFile(path);
   std::unique_ptr<SubUnit> subunit(new FileDefinedSubunit(definition));
   GraphID result = Add(subunit.get(), name, tag);
   modelCatalog[name] = definition.id;
   subunit.release();
   return result;
}


/*
    Adds a new structure to the world, think of this as a bag that contains a lot of connected sub-units, but gives them a new name.
    Using the structure we can continue to build structures by linking additional structures into it, analogously to linking
    sub-units to sub-units.

    It is not possible to addAndLink a sub-unit to a structure, but the single sub-unit can be wrapped in a structure and
    linked that way.

    Arguments: graphID existing graph in the world
               name a new unique name for this structure

    Returns graphID.
*/

GraphID World::Add(GraphID gid, structName name)
{
try{
    if (!GLEX->testnamestring(name))      throw SEBException("Bad symbol in structure name:"+name+". "+GLEX->invalid_name_message(name));
    if (!testGraphID(gid))                throw SEBException("Bad graphid:"+to_string(gid) );
    if (hasName(name))                    throw SEBException("Name "+name+" already exists in the world");

    Structure *struc = new Structure(name, gid, GLEX);
    totalNumberofGraphs++;

    list<structName> struccontainer;
    struccontainer.push_back(name);
    subGraphs.insert(pair<GraphID, list<structName>>(totalNumberofGraphs, struccontainer));

    nameCatalog.insert(pair<string, ABSSubUnit *>(name, struc));
    typeCatalog.insert(pair<string, int>(name, struc->getType()));

    return totalNumberofGraphs;   // return new graphID
}
catch (SEBException& e)
{
   e.PushCallStack("GraphID World::Add(GraphID "+to_string(gid)+", structName "+name+")");
   throw;
}

}

/* Adds a structure based on a name, which can be the name of an existing sub-unit or structure. */
/*GraphID World::Add(string oldname, string newname)
{
try{
     return Add(getGraphID(oldname), newname);
   }
catch (SEBException& e)
{
   e.PushCallStack("GraphID World::Add(string "+oldname+", string "+newname+")");
   throw;
}
}
*/


/* addsAndLink a new sub-unit to an existing sub-unit, thus growing the graph by one sub-unit.

   Arguments:
     sub  is the pointer to an uninitialized sub-unit object.
     newr has the form name.ref     name is given to the new sub-unit, ref has to be a well defined reference point on that type of sub-unit given by *sub.
     oldr has the form name.ref     name has to exist already in the world, and ref has to be a well defined reference point on that sub-unit.
     tag is an optional string.     Used for defining groups of identical sub-units.

     The new sub-unit is then linked to the existing (old) sub-unit by the two reference points.


   Returns graphID which is the same for the new and old sub-units.
*/

GraphID World::Link(SubUnit *sub, refPoint newr, refPoint oldr, string tag)
{
try{
    string newname = getName(newr); // extract name from  name.ref
    string oldname = getName(oldr);

    if (!GLEX->testnamestring(newname))                  throw SEBException("Bad symbol in sub-unit name:"+newname+". "+GLEX->invalid_name_message(newname));
    if (!tag.empty() && !GLEX->testnamestring(tag))      throw SEBException("Bad symbol in sub-unit tag:"+tag+". "+GLEX->invalid_name_message(tag));
    if (hasName(newname) )                               throw SEBException("Name "+newname+" already exists in the world");
    if (!hasName(oldname))                               throw SEBException("Name "+oldname+" does not exist in the world");

// Initialize sub-unit, this also sets up the known reference points for this type.
    sub->Init(newname, tag,  GLEX);

// Test that both reference points are good.
    string newref = getReferenceBase(newr); // extract refbase   from e.g. name.refbase#mypoint
    string oldref = getReferenceBase(oldr);

// test new reference name.
    if (!sub->hasReference(newref))      throw SEBException("Subunit named "+newname+" does not have reference point '"+newref+"'. Valid reference points: "+subunitReferenceList(sub));
    if (sub->hasDistributedReference(newref) && !hasAHash(newr))
                                         throw SEBException("Missing #... for distributed reference point: "+newref);

// rest old reference name
    SubUnit* olds = getSubunit(oldname);
    if (!olds->hasReference(oldref))     throw SEBException("Subunit '"+oldname+"' does not have reference point '"+oldref+"'. Valid reference points: "+subunitReferenceList(olds));
    if (olds->hasDistributedReference(oldref) && !hasAHash(oldr))
                                         throw SEBException("Missing #... for distributed reference point: "+oldref);
// Handle .distributedref#myname
    if(hasAHash(newr))
       {
         refPoint random=getReferenceAfterHash(newr);                // get mypoint
         sub->addRandomDistributedReferencePoint(newref, random);    // add if not already exists.
       }

    if(hasAHash(oldr))
       {
         refPoint random=getReferenceAfterHash(oldr);
         olds->addRandomDistributedReferencePoint(oldref, random);
       }

    GraphID gId = getGraphID(oldname);

    nameCatalog.insert(pair<string, ABSSubUnit *>(newname, sub));
    typeCatalog.insert(pair<string, int>(newname, sub->getType()));
    subGraphs.find(gId)->second.push_back(newname);
    links.push_back(generateLink(oldr, newr));

    return gId;
}
catch (SEBException& e)
{
   e.PushCallStack("GraphID World::Link(SubUnit*, refPoint "+newr+", refPoint "+oldr+", string tag"+tag+")"  );
   throw;
}

}

// helper
GraphID World::Link(string subtype, refPoint newr, refPoint oldr, string tag)
{
   std::unique_ptr<SubUnit> subunit(subunitRegistry.Create(subtype));
   GraphID result = Link(subunit.get(), newr, oldr, tag);
   string canonical = subtype;
   for (const auto& model : subunitRegistry.List()) {
       if (model.id == subtype || model.apiName == subtype) {
           canonical = model.id;
           break;
       }
   }
   modelCatalog[getName(newr)] = canonical;
   subunit.release();
   return result;
}

GraphID World::LinkFile(string path, refPoint newr, refPoint oldr, string tag)
{
   pyseb::SubunitDefinition definition = pyseb::LoadSubunitDefinitionFile(path);
   std::unique_ptr<SubUnit> subunit(new FileDefinedSubunit(definition));
   GraphID result = Link(subunit.get(), newr, oldr, tag);
   modelCatalog[getName(newr)] = definition.id;
   subunit.release();
   return result;
}

pyseb::SubunitModelInfo World::RegisterSubunitFile(string path)
{
   return subunitRegistry.RegisterFile(path);
}

vector<pyseb::SubunitModelInfo> World::RegisterSubunitDirectory(string path)
{
   return subunitRegistry.RegisterDirectory(path);
}

vector<pyseb::SubunitModelInfo> World::ListSubunitModels() const
{
   return subunitRegistry.List();
}

/*
    Grows a structure by linking the structure in gid to it.

    Instead of a subunit pointer, the method takes a graphID.
    The first reference point r1 has structure     mystructure:..:sub-unit.refpoint   such that the new structure is named mystructure.
    The second reference point r2  is of the form structure:..:sub-unit.refpoint      and should already exist in a previously defined structure.

    Returns graphID which is the same for the new and old structures.

*/
GraphID World::Link(GraphID gid, refPoint r1, refPoint r2)
{
try{
    if (!testGraphID(gid))  throw SEBException("Bad graphid:"+to_string(gid) );

// Check names of structures, at the highest level
    structName newname = prefix(r1);
    structName oldname = prefix(r2);

    if (hasName(newname))  throw SEBException("Name "+newname+" already exists in the world" );
    if (!hasName(oldname)) throw SEBException("Name "+oldname+" does not exist in the world" );

// Test entire path down to sub-unit and reference point base are good. true=> skip test of first structure in r1, since it has not yet been created.
    testPathSyntax(r1,true);
    testPathSyntax(r2);

// could be a new named distributed reference point  subunit.distributedreftype#myref  where we have to handle the unknown myref
    string newsubname = getName(r1);
    string oldsubname = getName(r2);
    string newref = getReference(r1);
    string oldref = getReference(r2);
    string newrefbase = getReferenceBase(r1);
    string oldrefbase = getReferenceBase(r2);

// testPathSyntax ensures that the path is correct and that the reference point exists.
// BUT we could still have a case where the reference point base is distributed, but the user has failed to specify a #.
    if (getSubunit(newsubname)->hasDistributedReference(newrefbase) && !hasAHash(newref))
        throw SEBException("Missing #... for distributed reference point: "+newref);
    if (getSubunit(oldsubname)->hasDistributedReference(oldrefbase) && !hasAHash(oldref))
        throw SEBException("Missing #... for distributed reference point: "+oldref);

    if(hasAHash(newref))                                                        // e.g. contour #mypoint
        {
            refPoint refbase=getReferenceBase(r1);                          // contour
            refPoint random=getReferenceAfterHash(r1);                      // mypoint
            getSubunit(newsubname)->addRandomDistributedReferencePoint(refbase, random);    // add if not already exists. Throw if refbase does not exist
        }

    if(hasAHash(oldref))                                                        // e.g. contour #mypoint
        {
            refPoint refbase=getReferenceBase(r2);                          // contour
            refPoint random=getReferenceAfterHash(r2);                      // mypoint
            getSubunit(oldsubname)->addRandomDistributedReferencePoint(refbase, random);    // add if not already exists.;
        }

    Structure *struc = new Structure(newname, gid, GLEX);

    GraphID oldgraphID = getGraphID(oldname);
    subGraphs.find(oldgraphID)->second.push_back(newname);

    nameCatalog.insert(pair<string, Structure *>(newname, struc));
    typeCatalog.insert(pair<string, int>(newname, struc->getType()));
    links.push_back(generateLink(r1, r2));

    GraphID newgraphID = getGraphID( prefix(r1) );
    if (newgraphID!= oldgraphID)
        throw SEBException("Internal error: mismatch in graphids");

    return newgraphID;
}
catch (SEBException& e)
{
   e.PushCallStack("GraphID World::Link(GraphID "+to_string(gid)+", refPoint "+r1+", refPoint "+r2+")" );
   throw;
}

}

pyseb::WorldDefinition World::Describe() const
{
    pyseb::WorldDefinition result;
    result.schemaVersion = 1;
    result.worldId = worldId;

    map<GraphID, string> labels;
    set<string> usedLabels;
    for (const auto& graph : subGraphs) {
        const auto known = graphLabels.find(graph.first);
        if (known != graphLabels.end()) {
            labels.emplace(graph.first, known->second);
            usedLabels.insert(known->second);
        }
    }
    int generatedLabel = 1;
    for (const auto& graph : subGraphs) {
        if (labels.count(graph.first)) continue;
        string candidate;
        do {
            candidate = string("graph-") + to_string(generatedLabel++);
        } while (usedLabels.count(candidate));
        labels.emplace(graph.first, candidate);
        usedLabels.insert(candidate);
    }
    auto graphLabel = [&](GraphID gid) { return labels.at(gid); };
    auto endpointObject = [](const string& endpoint) {
        const size_t colon = endpoint.find(':');
        const size_t period = endpoint.find('.');
        size_t cut = string::npos;
        if (colon != string::npos) cut = colon;
        if (period != string::npos) cut = min(cut, period);
        return cut == string::npos ? endpoint : endpoint.substr(0, cut);
    };
    auto graphForObject = [&](const string& name) {
        for (const auto& graph : subGraphs) {
            for (const auto& member : graph.second) {
                if (member == name) return graph.first;
            }
        }
        throw SEBException("The specified name was not found in any sub-graph", "World::Describe()");
    };

    for (const auto& item : nameCatalog) {
        if (item.second->getType() == STRUCTURE) {
            Structure* structure = reinterpret_cast<Structure*>(item.second);
            pyseb::WorldStructureDefinition value;
            value.name = item.first;
            value.graph = graphLabel(structure->getGraphID());
            result.structures.push_back(value);
            continue;
        }

        SubUnit* subunit = reinterpret_cast<SubUnit*>(item.second);
        if (subunit->getSubunitType() != FILEDEFINEDSUBUNIT) {
            throw SEBException(
                "Subunit '" + item.first + "' is not representable in world schema version 1",
                "World::Describe()");
        }
        pyseb::WorldSubunitDefinition value;
        value.name = item.first;
        const auto origin = modelCatalog.find(item.first);
        if (origin != modelCatalog.end()) {
            value.model = origin->second;
        } else {
            FileDefinedSubunit* fileDefined = dynamic_cast<FileDefinedSubunit*>(subunit);
            if (fileDefined == nullptr) {
                throw SEBException(
                    "Subunit '" + item.first + "' is not representable in world schema version 1",
                    "World::Describe()");
            }
            value.model = fileDefined->getDefinition().id;
        }
        value.tag = subunit->getTag();
        result.subunits.push_back(value);
    }

    map<GraphID, std::size_t> graphIndices;
    for (const auto& graph : subGraphs) {
        pyseb::WorldGraphDefinition value;
        value.id = graphLabel(graph.first);
        for (const auto& member : graph.second) value.members.push_back(member);
        graphIndices.emplace(graph.first, result.graphs.size());
        result.graphs.push_back(value);
    }

    for (const auto& link : links) {
        const GraphID firstGraph = graphForObject(endpointObject(link.first));
        const GraphID secondGraph = graphForObject(endpointObject(link.second));
        if (firstGraph != secondGraph) {
            throw SEBException("Link endpoints belong to different graphs", "World::Describe()");
        }
        result.graphs.at(graphIndices.at(firstGraph)).links.push_back(link);
    }

    set<GraphID> referencedGraphs;
    for (const auto& item : nameCatalog) {
        if (item.second->getType() == STRUCTURE) {
            referencedGraphs.insert(reinterpret_cast<Structure*>(item.second)->getGraphID());
        }
    }
    for (const auto& graph : subGraphs) {
        if (!referencedGraphs.count(graph.first)) result.roots.push_back(graphLabel(graph.first));
    }
    return result;
}

/*
Prints out all nested structures in a structure using a directory format.

Here casting is used to get a structure pointer from namecatalog that contains ABSSubunits pointers
 https://stackoverflow.com/questions/36802354/print-binary-tree-in-a-pretty-way-using-c
 Prints the the sub units with a print inspired by folder structure
 */

void World::folderprint(string name, string StructureString, bool isLast)
{
    if (isStructure(name))
       {
            cout << StructureString;
            cout << (isLast ? "└──" : "├──");
            cout << name + "\n";

            GraphID gid=getStructure(name)->getGraphID();                           // get gid inside "name".
            for (auto it = subgraph_cbegin(gid); it != subgraph_cend(gid); it++)    // traverse that subgraph
                {
                    string n=*it;
                    auto itt=it; itt++;  // Next element in list.
                    folderprint(n, StructureString + (isLast ? "    " : "│   "), itt==subgraph_cend(gid) );
                }
       }
     else if (isSubunit(name))
       {
            cout << StructureString;
            cout << (isLast ? "└──" : "├──");
            cout << name << "\n";
       }
     else
       throw SEBException("Should never happen!","void World::print(folderprintname=\""+name+"\", StructureString=\""+StructureString+"\", isLast="+to_string(isLast)+")" );
}


// Prints out which sub-unit / structure names belong to which sub-graphs
void World::printGraphIDs()
{
   for (auto g : subGraphs)
     {
        cout << g.first << " -> ";
        for (auto n : g.second )     // g is pair<GraphID,list<string>>
           {
              cout << n;
              if (isStructure(n)) cout << "(" << getStructure(n)->getGraphID() << ")";
              cout << " ";
           }
        cout << endl;
     }

}


void World::printPath(ReferencePointList path, string prefix, string postfix)
{
   for (auto p : path)
         cout << prefix << p << postfix;
}


/*Prints all links*/
void World::printLinks()
{
    for (auto it : links)
    {
        cout.width(50);
        cout << it.first << right << "<--->" + it.second + "\n";
    }
}


void World::printList(ReferencePointList path, string prefix, string postfix)
{
   for (auto it=path.begin() ; it!= path.end(); ++it )
     {
        cout << prefix << *it;
        if (*it != path.back()) cout << postfix;
     }
}
