
#include <SprayParticles.H>
#include <Particles_F.H>

//
// Computes particle state source terms as well as field source terms from particles
//
void
SprayParticleContainer::ComputeParticleSource (const MultiFab& S, 
        MultiFab& S_src, int lev, Real dt, int ireg,  bool calc_field_src )
{
    // BL_PROFILE("SprayParticleContainer::ComputeParticleSource()");
    BL_ASSERT(Ucc.nGrow() > 0);
    BL_ASSERT(OK(true, lev, S.nGrow()-1));
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    BL_ASSERT(!Ucc.contains_nan());


    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    BL_ASSERT(OnSameGrids(lev,Ucc));

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};


    // A level of particles is stored in a map indexed by the grid number.
    // typedef typename std::map<int,PBox> PMap;
    // A PBox is the storage for particles - something like
    // std::vector<Particle<BL_SPACEDIM +3,0> >
    PMap& pmap = m_particles[lev];

    for (auto& kv : pmap)
    {
        const int grid = kv.first;
        PBox&     pbox = kv.second;
        const int n    = pbox.size();

        const FArrayBox& field_fab = S[grid];
        FArrayBox& fld_src_fab = S_src[grid];

        fld_src_fab.setVal(0.0); // Here or outside and just add to it?

        // Number of field components
        int ncomp = field_fab.nComp();

        // Scratch arrays that will be passed to FORTRAN routine

        Box bx(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(n-1,0,0)));

        // to hold particle state
        FArrayBox pstate_SOA(bx,nstate);

        // to hold field variables interpolated to particle positions
        std::vector<Real> state_p(ncomp); // 1 particle
        FArrayBox fld_at_part_SOA(bx,ncomp); // Several particles

        // to hold particle source terms on return from FORTRAN routine
        FArrayBox pstate_src_SOA(bx,nstate);

        // to hold field source terms at particle positions
        FArrayBox fld_src_at_part_SOA(bx,ncomp);

        int state_start = state_start_idx(ireg);
        int src_start = src_start_idx(ireg);

        std::cout << "offset (state, source) stored: " << state_start << ", " << src_start << std::endl;
#ifdef _OPENMP
#pragma omp parallel 
#endif
        // This is the loop over particles to populate arrays
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            // Interpolate field data, pack into SoA
            // This interpolates the first 2-3 fields in the field_fab onto the particle locations
            // in p->m_pos

            // std::cout << i << " old x " << p.m_pos[0] << "; new x " << p.m_data[state_start] << std::endl;
            D_TERM( p.m_pos[0] = p.m_data[state_start];,
                    p.m_pos[1] = p.m_data[state_start+1];,
                    p.m_pos[2] = p.m_data[state_start+2];);

            // std::cout << "BL_SPACEDIM = " << BL_SPACEDIM << " ncomp = " << ncomp << std::endl;

            ParticleBase::Interp(p, geom, field_fab, idx, 
                    state_p.data(), BL_SPACEDIM);

            // std::cout << " old x " << p.m_pos[0] << "; new x " << p.m_data[state_start] << std::endl;
            // std::cout << p.m_id << ", " << i << " Interp state_p : " << state_p[0] << " ," << state_p[1] << std::endl;

            for (int icomp = 0; icomp < ncomp; icomp++)
            {
                fld_at_part_SOA(IntVect(D_DECL(i,0,0)),icomp) = state_p[icomp];
            }

            // This just straight up AoS to SoA conversion 
            for (int istate = 0; istate < nstate; istate++)
            {
                pstate_SOA(IntVect(D_DECL(i,0,0)),istate)
                    = p.m_data[state_start+istate];
            }

        }

#if 0
        // Next call a bunch of FORTRAN to compute RHS for spray particle state update
        spray_part_src ( pstate_SOA.dataPtr(), fld_at_part_SOA.dataPtr(),
			 pstate_src_SOA.dataPtr(),
			 &n, &nstate, &ncomp);

        spray_part_src_to_fld_src ( pstate_SOA.dataPtr(), fld_at_part_SOA.dataPtr(),
				    pstate_src_SOA.dataPtr(), fld_src_at_part_SOA.dataPtr(), 
				    &n, &nstate, &ncomp);
#endif


        std::cout << "Not using fortran, hack" << std::endl;
        // Unpack data for this group of particles back into AoS and deposit field src onto grid
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            // This just straight up SoA to AoS conversion 
            // Store particle state source terms in NR:2*NR-1 
            // TODO(RG) : Fix this
            for (int istate = 0; istate < nstate; istate++)
            {
                p.m_data[src_start+istate] = 
                    fld_at_part_SOA(IntVect(D_DECL(i,0,0)), istate);
                //    = pstate_src_SOA(IntVect(D_DECL(i,0,0)),istate);
            }

            // std::cout << " Storing src starting at index: " << src_start << std::endl;
            // std::cout << src_start << " ; " <<  p.m_id << ", " << i << " src data: " << p.m_data[src_start] << ", " << p.m_data[src_start+1] << std::endl;


            // Depost field source at particle location onto S_Src
            for (int icomp = 0; icomp < ncomp; icomp++)
            {
                Real fld_src = fld_src_at_part_SOA(IntVect(D_DECL(i,0,0)),icomp);
                deposit(p, geom, fld_src_fab, fld_src, icomp);
            }
        }
    } // End loop over all the PBoxes for this level
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
                ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

                if (ParallelDescriptor::IOProcessor())
                {
                std::cout << "SprayParticleContainer::ComputeParticleSource() time: " << stoptime << '\n';
                }
#ifdef BL_LAZY
                });
#endif
    }
}


void
SprayParticleContainer::deposit(ParticleBase& prt, const Geometry& geom, 
        FArrayBox& fld_src_fab, Real fld_src, int icomp)
{
    return;

}

//
// Uses midpoint method to advance particles using cell-centered velocity
//
void
SprayParticleContainer::AdvectWithUcc (const MultiFab& Ucc, int lev, Real dt)
{
    BL_ASSERT(Ucc.nGrow() > 0);
    BL_ASSERT(OK(true, lev, Ucc.nGrow()-1));
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    BL_ASSERT(!Ucc.contains_nan());
    BL_ASSERT(SPRAY_COMPONENTS >= BL_SPACEDIM);

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    BL_ASSERT(OnSameGrids(lev,Ucc));

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};

    for (int ipass = 0; ipass < 2; ipass++)
    {
        PMap& pmap = m_particles[lev];

        for (auto& kv : pmap)
        {
            const int grid = kv.first;
            PBox&     pbox = kv.second;
            const int n    = pbox.size();

	    const FArrayBox& fab = Ucc[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                
                if (p.m_id <= 0) continue;

                BL_ASSERT(p.m_grid == grid);

		Real v[BL_SPACEDIM];

		ParticleBase::Interp(p, geom, fab, idx, v, BL_SPACEDIM);

		if (ipass == 0) {
		    //
		    // Save old position and the vel & predict location at dt/2.
		    //
		    for (int d = 0; d < BL_SPACEDIM; d++)
		    {
			p.m_data[d] = p.m_pos[d];
                        p.m_pos[d] += 0.5*dt*v[d];
                    }
		} else {
		    //
		    // Update to final time using the orig position and the vel at dt/2.
		    //
		    for (int d = 0; d < BL_SPACEDIM; d++)
		    {
                        p.m_pos[d]  = p.m_data[d] + dt*v[d];
                        // Save the velocity for use in Timestamp().
			p.m_data[d] = v[d];
                    }
                }
                
                ParticleBase::RestrictedWhere(p,m_gdb, Ucc.nGrow()); 
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "TracerParticleContainer::AdvectWithUcc() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}


//
// Set the first BL_SPACEDIM elements in the particle data to the position
//
void
SprayParticleContainer::PosToState (int ireg, int lev)
{
    BL_ASSERT(lev >= 0 && lev < m_particles.size());
    BL_ASSERT(SPRAY_COMPONENTS >= BL_SPACEDIM);

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};

    PMap& pmap = m_particles[lev];
    int state_start = state_start_idx(ireg);

    for (auto& kv : pmap)
    {
        const int grid = kv.first;
        PBox&     pbox = kv.second;
        const int n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                p.m_data[state_start + d] = p.m_pos[d];
            }

        }
    }
}


//
// Set particle position to the  first BL_SPACEDIM in the particle state
//
void
SprayParticleContainer::StateToPos (int ireg, int lev)
{
    BL_ASSERT(lev >= 0 && lev < m_particles.size());
    BL_ASSERT(SPRAY_COMPONENTS >= BL_SPACEDIM);

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};

    PMap& pmap = m_particles[lev];
    int state_start = state_start_idx(ireg);

    for (auto& kv : pmap)
    {
        const int grid = kv.first;
        PBox&     pbox = kv.second;
        const int n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                p.m_pos[d] = p.m_data[state_start + d];
            }

        }
    }
}

//
// Apply tentative update to particle state
//
void
SprayParticleContainer::Update( int ireg_src_terms, int ireg_start_state, int ireg_dest_state, int lev, Real time, Real dt)
{
    BL_ASSERT(lev >= 0 && lev < m_particles.size());
    BL_ASSERT(SPRAY_COMPONENTS >= BL_SPACEDIM);

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};

    PMap& pmap = m_particles[lev];
    int start_state = state_start_idx(ireg_start_state);
    int start_src = src_start_idx(ireg_src_terms);
    int dest_state = state_start_idx(ireg_dest_state);

    std::cout << "registers (state in, soruce, state dest) " << ireg_start_state << "," << ireg_src_terms<< "," << ireg_dest_state << std::endl;
    std::cout << "offset (state in, soruce, state dest) " << start_state << "," << start_src << "," << dest_state << std::endl;
    std::cout.flush();

    int nstate = nState();

    for (auto& kv : pmap)
    {
        const int grid = kv.first;
        PBox&     pbox = kv.second;
        const int n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            // std::cout << "idx: " << start_src << "; " << p.m_id << ", " << i << " update data: " << p.m_data[start_src] << ", " << p.m_data[start_src+1] << std::endl;
            // std::cout << " Update looking for src starting at index: " << start_src << std::endl;
            for (int istate = 0; istate < nstate; istate++)
            {
                p.m_data[dest_state+istate] = p.m_data[start_state+istate] 
                    + dt*p.m_data[start_src+istate];
            }
            // ParticleBase::RestrictedWhere(p,m_gdb, Ucc.nGrow()); 

        }
    }
}

