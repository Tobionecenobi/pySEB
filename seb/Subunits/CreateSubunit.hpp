/*
    Helper method for creating different subunit types

*/


#ifndef INCLUDE_GUARD_CREATESUBUNIT
#define INCLUDE_GUARD_CREATESUBUNIT

inline SubUnit* CreateSubunit(string subtype)
{
        if (subtype == "SymbolicSubunit")     return new SymbolicSubunit();
   else 
      throw SEBException("Unknown sub-unit type "+subtype,"Create");
}

#endif
