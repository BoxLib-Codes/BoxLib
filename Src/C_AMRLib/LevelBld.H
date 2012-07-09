
#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

#include <REAL.H>
#include <Box.H>
#include <Geometry.H>
#include <Amr.H>
#include <AmrRegion.H>

extern LevelBld* getLevelBld ();

//
// Builds problem-specific AmrRegions.
//
// Abstract base class specifying an interface for building problem-specific
// AmrRegions.  
//

class LevelBld
{
public:

    virtual ~LevelBld () {}
    //
    // Perform any problem-dependent setup such as physical
    // boundary condition and derived quantities.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void variableSetUp () = 0;
    //
    // Perform any problem-dependent cleanup.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual void variableCleanUp () = 0;
    //
    // This is a virtual constructor for types derived
    // from AmrRegion.  The derived type is initialized
    // with the default constructor.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual AmrRegion* operator() () = 0;
    //
    // This is a virtual constructor for types derived
    // from AmrRegion.  The derived type is initialized
    // with the five specified variables.
    // This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual AmrRegion* operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& geom_lev,
                                  const BoxArray& ba,
                                  Real            time) = 0;
};

#endif /*_LEVELBLD_H_*/
