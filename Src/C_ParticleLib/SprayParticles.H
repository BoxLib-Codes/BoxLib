#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include <Particles.H>
#include <Geometry.H>

// Use a particle container storing 2x the data
// 0:NR-1 particle state
// NR:2*NR-1 particle state source terms computed at src_time
template <int NR>
class SprayParticleContainer
    : public ParticleContainer<BL_SPACEDIM +2*NR,0,std::vector<Particle<BL_SPACEDIM +2*NR,0> > >
{
public:

    SprayParticleContainer<NR> (ParGDBBase* gdb)
	: ParticleContainer<BL_SPACEDIM +2*NR,0,std::vector<Particle<BL_SPACEDIM +2*NR,0> > >(gdb)
	{ nstate = NR; src_time = -1.0; }

    SprayParticleContainer<NR> (const Geometry            & geom, 
			     const DistributionMapping & dmap,
			     const BoxArray            & ba)
	: ParticleContainer<BL_SPACEDIM +2*NR,0,std::vector<Particle<BL_SPACEDIM +2*NR,0> > >(geom,dmap,ba)
	{ nstate = NR; }

    ~SprayParticleContainer () {}

    void ComputeParticleSource (const MultiFab& S, MultiFab& S_src, 
            int level, Real dt, bool calc_field_src = false);

    int nState() { return nstate; src_time = -1.0;}

	using ParticleContainer<BL_SPACEDIM +2*NR,0,std::vector<Particle<BL_SPACEDIM +2*NR,0> > >::PBox;
	using ParticleContainer<BL_SPACEDIM +2*NR,0,std::vector<Particle<BL_SPACEDIM +2*NR,0> > >::PMap;
	using ParticleContainer<BL_SPACEDIM +2*NR,0,std::vector<Particle<BL_SPACEDIM +2*NR,0> > >::ParticleType;

    void deposit(ParticleBase& prt, const Geometry& geom,
            FArrayBox& fld_src_fab, Real fld_src, int icomp);

private:
    int nstate;
    Real src_time;
};


#endif

