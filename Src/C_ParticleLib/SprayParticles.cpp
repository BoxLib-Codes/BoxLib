
#include <SprayParticles.H>
#include <Particles_F.H>

//
// Computes particle state source terms as well as field source terms from particles
//
void
SprayParticleContainer::ComputeParticleSource (const MultiFab& S, 
        MultiFab& S_src, int lev, Real dt, bool calc_field_src )
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
            ParticleBase::Interp(p, geom, field_fab, idx, 
                    state_p.data(), BL_SPACEDIM);

            for (int icomp = 0; icomp < ncomp; icomp++)
            {
                fld_at_part_SOA(IntVect(D_DECL(i,0,0)),icomp) = state_p[icomp];
            }

            // This just straight up AoS to SoA conversion 
            for (int istate = 0; istate < nstate; istate++)
            {
                pstate_SOA(IntVect(D_DECL(i,0,0)),istate)
                    = p.m_data[istate];
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


        // Unpack data for this group of particles back into AoS and deposit field src onto grid
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);
            // This just straight up SoA to AoS conversion 
            // Store particle state source terms in NR:2*NR-1 
            for (int istate = 0; istate < nstate; istate++)
            {
                p.m_data[istate+nstate]
                    = pstate_src_SOA(IntVect(D_DECL(i,0,0)),istate);
            }


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
