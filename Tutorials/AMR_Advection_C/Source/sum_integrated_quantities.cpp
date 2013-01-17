#include <iomanip>

#include <ADR.H>
#include <ADR_F.H>

void
ADR::sum_integrated_quantities ()
{
    if (verbose <= 0) return;

    int finest_level = master->finestLevel();
    Real dt_crse     = master->dtLevel(0);
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real xvel        = 0.0;
    Real yvel        = 0.0;
#if (BL_SPACEDIM==3)
    Real zvel        = 0.0;
#endif

    RegionIterator it = master->getRegionIterator();
    for ( ; !it.isFinished(); ++it)
    {
        ADR& adr_reg = get_region(it.getID());
        int lev = adr_reg.Level();

        mass     += adr_reg.volWgtSum("density", time);
        xvel     += adr_reg.volWgtSum("xvel", time);
        yvel     += adr_reg.volWgtSum("yvel", time);
#if (BL_SPACEDIM==3)
        zvel     += adr_reg.volWgtSum("zvel", time);
#endif
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << '\n';
        std::cout << "TIME= " << time << " MASS        = "   << mass  << '\n';
        std::cout << "TIME= " << time << " XVEL        = "   << xvel     << '\n';
        std::cout << "TIME= " << time << " YVEL        = "   << yvel     << '\n';
#if (BL_SPACEDIM==3)
        std::cout << "TIME= " << time << " ZVEL        = "   << zvel     << '\n';
#endif
	std::cout<<'\n';
    }
}

Real
ADR::volWgtSum (const std::string& name,
                   Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < master->finestLevel())
    {
        baf = master->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < master->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
	BL_FORT_PROC_CALL(SUMMASS,summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
