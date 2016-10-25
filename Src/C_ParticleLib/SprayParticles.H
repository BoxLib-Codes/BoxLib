#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include <Particles.H>
#include <Geometry.H>

// Use a particle container storing 2x the data
// 0:NR-1 particle state
// NR:2*NR-1 particle state source terms computed at src_time
class SprayParticleContainer
    : public ParticleContainer<2*SPRAY_COMPONENTS,0,std::vector<Particle<2*SPRAY_COMPONENTS,0> > >
{
public:

    SprayParticleContainer (ParGDBBase* gdb)
	: ParticleContainer<2*SPRAY_COMPONENTS,0,std::vector<Particle<2*SPRAY_COMPONENTS,0> > >(gdb)
	{ nstate = SPRAY_COMPONENTS; src_time = -1.0; }

    SprayParticleContainer (const Geometry            & geom, 
			     const DistributionMapping & dmap,
			     const BoxArray            & ba)
	: ParticleContainer<2*SPRAY_COMPONENTS,0,std::vector<Particle<2*SPRAY_COMPONENTS,0> > >(geom,dmap,ba)
	{ nstate = SPRAY_COMPONENTS; }

    ~SprayParticleContainer () {}

    void ComputeParticleSource (const MultiFab& S, MultiFab& S_src, 
            int level, Real dt, bool calc_field_src = false);

    void AdvectWithUcc (const MultiFab& ucc, int level, Real dt);

    int nState() { return nstate; src_time = -1.0;}

	// using ParticleContainer<2*NR,0,std::vector<Particle<2*NR,0> > >::PBox;
	// using ParticleContainer<2*NR,0,std::vector<Particle<2*NR,0> > >::PMap;
	// using ParticleContainer<2*NR,0,std::vector<Particle<2*NR,0> > >::ParticleType;

    void deposit(ParticleBase& prt, const Geometry& geom,
            FArrayBox& fld_src_fab, Real fld_src, int icomp);

    void test_fcn();

private:
    int nstate;
    Real src_time;
};


#endif

