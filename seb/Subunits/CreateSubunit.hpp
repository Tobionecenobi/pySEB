/*
    Helper method for creating different subunit types

*/


#ifndef INCLUDE_GUARD_CREATESUBUNIT
#define INCLUDE_GUARD_CREATESUBUNIT

inline SubUnit* CreateSubunit(string subtype)
{
        if (subtype == "GaussianLoop")        return new GaussianLoop();
   else if (subtype == "ThinRod")             return new ThinRod();
   else if (subtype == "ThinCircle")          return new ThinCircle();
   else if (subtype == "ThinDisk")            return new ThinDisk();
   else if (subtype == "ThinSphericalShell")  return new ThinSphericalShell();
   else if (subtype == "SolidSphere")         return new SolidSphere();
   else if (subtype == "SolidSphericalShell") return new SolidSphericalShell();
   else if (subtype == "SolidCylinder")       return new SolidCylinder();
   else if (subtype == "SymbolicSubunit")     return new SymbolicSubunit();
   else 
      throw SEBException("Unknown sub-unit type "+subtype,"Create");
}

#endif
