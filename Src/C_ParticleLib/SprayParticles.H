#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include <Particles.H>
#include <Geometry.H>

#define NREG 2 // Number of copies of state + source terms necessary for time advance

// Use a particle container storing 2x the data
// 0:NR-1 particle state
// NR:2*NR-1 particle state source terms computed at src_time
class SprayParticleContainer
    : public ParticleContainer<2*NREG*SPRAY_COMPONENTS,0,std::vector<Particle<2*NREG*SPRAY_COMPONENTS,0> > >
{
public:

    SprayParticleContainer (ParGDBBase* gdb)
	: ParticleContainer<2*NREG*SPRAY_COMPONENTS,0,std::vector<Particle<2*NREG*SPRAY_COMPONENTS,0> > >(gdb)
	{ nstate = SPRAY_COMPONENTS; 
        for (int i=0; i<NREG; ++i){
            src_time[i] = -1.0; } }

    SprayParticleContainer (const Geometry            & geom, 
			     const DistributionMapping & dmap,
			     const BoxArray            & ba)
	: ParticleContainer<2*NREG*SPRAY_COMPONENTS,0,std::vector<Particle<2*NREG*SPRAY_COMPONENTS,0> > >(geom,dmap,ba)
	{ nstate = SPRAY_COMPONENTS; }

    ~SprayParticleContainer () {}

    void ComputeParticleSource (const MultiFab& S, MultiFab& S_src, 
            int level, Real dt, int ireg, bool calc_field_src = false);

    void AdvectWithUcc (const MultiFab& ucc, int level, Real dt);

    int nState() { return nstate;}

    void deposit(ParticleBase& prt, const Geometry& geom,
            FArrayBox& fld_src_fab, Real fld_src, int icomp);

    void test_fcn();

    void set_src_time(int ireg, Real time){ src_time[ireg] = time; }

    int state_start_idx (int ireg) { return (ireg*2*SPRAY_COMPONENTS); }
    int src_start_idx (int ireg) { return ((ireg*2+1)*SPRAY_COMPONENTS); }

    void PosToState(int ireg, int lev);

    void StateToPos(int ireg, int lev);

    void Update( int ireg_src_terms, int ireg_start_state, int ireg_dest_state, int lev, Real time, Real dt);

private:
    int nstate;
    Real src_time[NREG];
};


#endif

