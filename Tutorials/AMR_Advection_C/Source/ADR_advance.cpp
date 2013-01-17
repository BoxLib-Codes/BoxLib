#include <winstd.H>

#include "ADR.H"
#include "ADR_F.H"

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

using std::string;

Real
ADR::advance (Real time,
              Real dt,
              int  iteration,
              int  ncycle)
{
    const Real strt = ParallelDescriptor::second();

    const int finest_level = master->finestLevel();
    int finest_level_to_advance;
    bool nosub = !master->subCycle();
    
    ExecutionTree advance_tree(m_id);
    if (nosub)
    {
        if (level > 0)
            return dt;
        advance_tree.buildFromStructure(master->getRegions().getStructure(m_id));
        finest_level_to_advance = finest_level;
        RegionIterator it = master->getRegionIterator(m_id);
        for ( ; !it.isFinished(); ++it)
        {
            ID id = it.getID();
            advance_tree.setData(id, 1);
        }
    }
    else
    {
        if (strict_subcycling)
        {
            finest_level_to_advance = level;
            //Leave advance tree as is.
        }
        else
        {
            // Check if the parent region advanced this.
            if (level > 0 && master->dtRegion(m_id) == master->dtRegion(parent_region->getID()))
                return dt;
            // Find all the subregions with the same timestep.
            finest_level_to_advance = level;
            advance_tree.buildFromStructure(master->getRegions().getStructure(m_id));
            RegionIterator it = master->getRegionIterator(m_id);
            for( ; !it.isFinished(); ++it)
            {
                ID id = it.getID();
                if (master->dtRegion(id) == dt)
                {
                    advance_tree.setData(id,1);
                    if (it.getLevel() > finest_level_to_advance)
                        finest_level_to_advance = it.getLevel();
                }
                else
                {
                    advance_tree.setData(id,0);
                }
            }
            advance_tree.prune();
        }
    }

    //
    // Move current data to previous, clear current.
    // Don't do this if a coarser level has done this already.
    //
    if (level == 0 || iteration > 1)
    {
        for (int lev = level; lev <= finest_level; lev++)
        {
            RegionIterator it = master->getRegionIterator(m_id, lev);
            for ( ; !it.isFinished(); ++it)
            {
                ID id = it.getID();
                Real dt_reg = master->dtRegion(id);
                for (int k = 0; k < NUM_STATE_TYPE; k++)
                {
                    get_region(id).state[k].allocOldData();
                    get_region(id).state[k].swapTimeLevels(dt_reg);
                }
            }
        }
    }

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();

    //
    // Call the advance at each level to be advanced
    //
    for (int lev = level; lev <= finest_level_to_advance; lev++)
    {
        ExecutionTreeIterator it = advance_tree.getIterator(m_id, lev);
        for ( ; !it.isFinished(); ++it)
        {
            get_region(it.getID()).advance_level(time, dt);
        }
    }
    
    //
    // We must reflux here
    //
    if (do_reflux)
    {
        for (int lev = level; lev < finest_level_to_advance; lev++)
        {
            ExecutionTreeIterator it = advance_tree.getIterator(m_id, lev);
            for ( ; !it.isFinished(); ++it)
                get_region(it.getID()).reflux(1);
        }
    }
    
    // Always average down from finer to coarser.
    for (int lev = finest_level_to_advance - 1; lev >= level; lev--)
    {
        ExecutionTreeIterator it = advance_tree.getIterator(m_id, lev);
        for ( ; !it.isFinished(); ++it)
            get_region(it.getID()).average_down(true);
    }

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Time  at end of routine " << end << '\n';
    }
    
    return dt;
}

Real
ADR::advance_level (Real time,
                    Real dt)
{
    const Real prev_time    = state[State_Type].prevTime();
    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = master->finestLevel();
    MultiFab&  S_old        = get_old_data(State_Type);
    MultiFab&  S_new        = get_new_data(State_Type);

    if (std::abs(time-prev_time) > (1.e-10*cur_time) )
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "advance_level:  prev_time = " << prev_time << std::endl;
            std::cout << "advance_level:       time = " <<      time << std::endl;
        }
        BoxLib::Abort("time should equal prev_time in advance_level!");
    }

#ifndef NDEBUG
    if (S_old.contains_nan(Density, S_old.nComp(), 0))
    {
        for (int i = 0; i < S_old.nComp(); i++)
        {
           if (ParallelDescriptor::IOProcessor())
                std::cout << "advance_level: testing component " << i << " for NaNs" << std::endl;
            if (S_old.contains_nan(Density+i,1,0))
                BoxLib::Abort("S_old has NaNs in this component");
        }
    }
#endif

    // It's possible for interpolation to create very small negative values for
    // species so we make sure here that all species are non-negative after this
    // point
    enforce_nonnegative_species(S_old);

    if (do_reflux && level < finest_level)
    {
        //
        // Set reflux registers to zero.
        //
        PList<AmrRegion> children;
        master->getRegions().getChildrenOfNode(m_id, children);
        for (PList<AmrRegion>::iterator it = children.begin(); it != children.end(); ++it)
        {
            ID c_id = (*it)->getID();
            get_flux_reg(c_id).setVal(0);
        }
    }

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    /// We've updated fine to be a tested bool since there is no longer a guarantee of
    /// a single fine flux reg
    //FluxRegister* fine    = 0;
    bool fine = false;
    FluxRegister* current = 0;

    if (do_reflux && level < finest_level)
        fine = true;
        //fine = &get_flux_reg(level+1);
    if (do_reflux && level > 0)
        current = &get_flux_reg(m_id);

    const Real* dx     = geom.CellSize();
    Real        courno = -1.0e+200;

#ifdef REACTIONS
    react_first_half_dt(S_old,time,dt);
#endif

    FArrayBox flux[BL_SPACEDIM], u_gdnv[BL_SPACEDIM];

    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux && fine)
    {
        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            BoxArray ba = S_new.boxArray();
            ba.surroundingNodes(j);
            fluxes[j].define(ba, NUM_STATE, 0, Fab_allocate);
        }
    }

    MultiFab ext_src_old(grids,NUM_STATE,1,Fab_allocate);
    ext_src_old.setVal(0.0);

    if (add_ext_src)
       getOldSource(prev_time,dt,ext_src_old);

#ifdef DIFFUSION
    MultiFab OldDiffTerm(grids,NUM_STATE,1);
    OldDiffTerm.setVal(0.);
    if (diffuse_spec == 1) {
       for (int ispec = 0; ispec < NumSpec; ispec++)
          add_diffusion_to_old_source(ext_src_old,OldDiffTerm,prev_time,FirstSpec+ispec);
    }
#endif
    ext_src_old.FillBoundary();
        
    for (FillPatchIterator fpi(*this, S_new, NUM_GROW, time, State_Type, 0, NUM_STATE);
         fpi.isValid(); ++fpi)
    {
        int mfiindex = fpi.index();

        Box bx(fpi.UngrownBox());

        // Create FAB for extended grid values (including boundaries) and fill.
        FArrayBox &state = fpi();
        FArrayBox &stateout = S_new[fpi];
        
        // Allocate fabs for fluxes.
        for (int i = 0; i < BL_SPACEDIM ; i++)  
            flux[i].resize(BoxLib::surroundingNodes(bx,i),NUM_STATE);

        // And for Godunov.
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            Box ebx = bx;
            ebx.surroundingNodes(i);
            ebx.grow(1);
            u_gdnv[i].resize(ebx, 1);
            u_gdnv[i].setVal(1.e200);
        }

        const int*  domain_lo = geom.Domain().loVect();
        const int*  domain_hi = geom.Domain().hiVect();

        Real cflLoc = -1.0e+20;

        BL_FORT_PROC_CALL(ADVECT,advect)
                (&time,
                 bx.loVect(), bx.hiVect(),
                 domain_lo, domain_hi,
                 BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
		 BL_TO_FORTRAN(u_gdnv[0]),
		 BL_TO_FORTRAN(u_gdnv[1]),
#if (BL_SPACEDIM == 3)
		 BL_TO_FORTRAN(u_gdnv[2]),
#endif
                 BL_TO_FORTRAN(ext_src_old[fpi]),
                 D_DECL(BL_TO_FORTRAN(flux[0]), 
                        BL_TO_FORTRAN(flux[1]), 
                        BL_TO_FORTRAN(flux[2])), 
                 dx, &dt, verbose);

            if (do_reflux)
            {
                if (fine)
                {
                    for (int i = 0; i < BL_SPACEDIM ; i++)
                        fluxes[i][mfiindex].copy(flux[i]);
                }
                if (current)
                {
                    for (int i = 0; i < BL_SPACEDIM ; i++)
                        current->FineAdd(flux[i],i,mfiindex,0,0,NUM_STATE,1);
                }
            }

            courno = std::max(courno, cflLoc);
        }

        if (add_ext_src)  
        {
           getOldSource(prev_time,dt,ext_src_old);
           ext_src_old.mult(-0.5*dt);

           // Compute source at new time (no ghost cells needed)
           MultiFab ext_src_new(grids,NUM_STATE,0,Fab_allocate);
           ext_src_new.setVal(0.0);

           getNewSource(prev_time,cur_time,dt,ext_src_new);
           ext_src_new.mult(0.5*dt);

           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
           MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
        }

    for (int i = 0; i < BL_SPACEDIM ; i++)
    {
        flux[i].clear();
        u_gdnv[i].clear();
    }

    if (do_reflux && fine)
    {
        PList<AmrRegion> children;
        master->getRegions().getChildrenOfNode(m_id, children);
        for (PList<AmrRegion>::iterator it = children.begin(); it != children.end(); ++it)
        {
            ID c_id = (*it)->getID();
            for (int i = 0; i < BL_SPACEDIM ; i++)
            {
                get_flux_reg(c_id).CrseInit(fluxes[i], i, 0, 0, NUM_STATE, -1);
            }
        }
        for (int i = 0; i < BL_SPACEDIM ; i++)
        {
            fluxes[i].clear();
        }
    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "OOPS -- EFFECTIVE CFL AT THIS LEVEL " << level
                      << " IS " << courno << '\n';

        BoxLib::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

#ifndef NDEBUG
    if (S_new.contains_nan(Density, S_new.nComp(), 0))
    {
        for (int i = 0; i < S_new.nComp(); i++)
        {
            if (S_new.contains_nan(Density + i, 1, 0))
            {
                BoxLib::Abort("S_new has NaNs after the update");
            }
        }
    }
#endif


#ifdef DIFFUSION
     if (diffuse_spec == 1) {
        for (int ispec = 0; ispec < NumSpec; ispec++)
           time_center_diffusion(S_new, OldDiffTerm, cur_time, dt, FirstSpec+ispec);
     }
#endif

#ifdef REACTIONS
    react_second_half_dt(S_new,cur_time,dt);
#endif

    // Copy old velocity into new velocity -- assumes velocity unchanging.
    MultiFab::Copy(S_new,S_old,Xvel,Xvel,BL_SPACEDIM,0);

    return dt;
}
