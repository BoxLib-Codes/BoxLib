
#include "LevelBld.H"
#include "ADR.H"

class ADRBld
    :
    public LevelBld
{
    virtual void variable_setup();
    virtual void variable_cleanup();

    // hack copies for BoxLib overriding
    virtual void variableSetUp();
    virtual void variableCleanUp();

    virtual AmrRegion *operator() ();
    virtual AmrRegion *operator() (Amr& papa, ID id,
                                  const Geometry& level_geom,
                                  const BoxArray& ba, Real time);
};

ADRBld ADR_bld;

LevelBld*
get_level_bld ()
{
    return &ADR_bld;
}

void
ADRBld::variable_setup ()
{
    ADR::variable_setup();
}

void
ADRBld::variable_cleanup ()
{
    ADR::variable_cleanup();
}

AmrRegion*
ADRBld::operator() ()
{
    return new ADR;
}

AmrRegion*
ADRBld::operator() (Amr&            papa,
                    ID              id,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    Real            time)
{
    return new ADR(papa, id, level_geom, ba, time);
}

// override hacks, copies of above
LevelBld*
getLevelBld()
{
    return &ADR_bld;
}

void ADRBld::variableSetUp()
{
    ADR::variable_setup();
}
void ADRBld::variableCleanUp()
{
    ADR::variable_cleanup();
}

