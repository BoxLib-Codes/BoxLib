#include <winstd.H>

#include <algorithm>
#include <cstdio>
#include <vector>
#include <iostream>
#include <string>
using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <Utility.H>
#include <CONSTANTS.H>
#include <ADR.H>
#include <ADR_F.H>
#include <Derive_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <ParmParse.H>

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

static int  sum_interval = -1;
static Real fixed_dt     = -1.0;
static Real initial_dt   = -1.0;
static Real dt_cutoff    = 0.0;

int          ADR::strict_subcycling = 0;

bool         ADR::dump_old      = false;

int          ADR::verbose       = 0;
int          ADR::show_timings  = 0;
Real         ADR::cfl           = 0.8;
Real         ADR::init_shrink   = 1.0;
Real         ADR::change_max    = 1.1;
ErrorList    ADR::err_list;
int          ADR::radius_grow   = 1;
BCRec        ADR::phys_bc;
int          ADR::NUM_STATE     = -1;
int          ADR::do_reflux     = 1;
int          ADR::NUM_GROW      = -1;

int          ADR::Density       = -1;
int          ADR::Xvel          = -1;
int          ADR::Yvel          = -1;
int          ADR::Zvel          = -1;


int          ADR::NumSpec       = 0;
int          ADR::FirstSpec     = -1;
int          ADR::LastSpec      = -1;

int          ADR::NumAdv        = 0;
int          ADR::FirstAdv      = -1;
int          ADR::LastAdv       = -1;

int          ADR::do_react = -1;
int          ADR::add_ext_src = 0;

#ifdef DIFFUSION
Diffusion*    ADR::diffusion  = 0;
int           ADR::diffuse_spec = 0;
#endif

int          ADR::allow_untagging = 0;
int          ADR::normalize_species = 0;

// Note: ADR::variableSetUp is in ADR_setup.cpp

void
ADR::variable_cleanup () 
{
#ifdef DIFFUSION
  if (diffusion != 0) {
    if (verbose && ParallelDescriptor::IOProcessor()) {
      cout << "Deleting diffusion in variable_cleanup..." << '\n';
    }
    delete diffusion;
    diffusion = 0;
  }
#endif

    desc_lst.clear();
}

void
ADR::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("act");   

    pp.query("v",verbose);
    verbose = (verbose ? 1 : 0);
    pp.get("init_shrink",init_shrink);
    pp.get("cfl",cfl);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("initial_dt",initial_dt);
    pp.query("sum_interval",sum_interval);
    pp.query("do_reflux",do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp.get("dt_cutoff",dt_cutoff);

    pp.query("strict_subcycling",strict_subcycling);

    pp.query("dump_old",dump_old);

    // Get boundary conditions
    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir<BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "ADR::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    BoxLib::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "ADR::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    BoxLib::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir=0; dir<BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "ADR::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                BoxLib::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "ADR::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                BoxLib::Error();
            }
        }
    }

    if ( Geometry::IsRZ() && (lo_bc[0] != Symmetry) ) {
        std::cerr << "ERROR:ADR::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
        BoxLib::Error();
    }

#if (BL_SPACEDIM == 1)
    if ( Geometry::IsSPHERICAL() )
    {
      if ( (lo_bc[0] != Symmetry) && (Geometry::ProbLo(0) == 0.0) ) 
      {
        std::cerr << "ERROR:ADR::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
        BoxLib::Error();
      }
    }
#elif (BL_SPACEDIM == 2)
        if ( Geometry::IsSPHERICAL() )
        {
          BoxLib::Abort("We don't support spherical coordinate systems in 2D");
        }
#elif (BL_SPACEDIM == 3)
        if ( Geometry::IsRZ() )
        {
          BoxLib::Abort("We don't support cylindrical coordinate systems in 3D"); 
        }
        else if ( Geometry::IsSPHERICAL() )
        {
          BoxLib::Abort("We don't support spherical coordinate systems in 3D");
        }
#endif

    pp.get("do_react",do_react);
    pp.query("add_ext_src",add_ext_src);

#ifdef DIFFUSION
    pp.query("diffuse_spec",diffuse_spec);
#endif

    pp.query("allow_untagging",allow_untagging);
    pp.query("normalize_species",normalize_species);
}

ADR::ADR ()
{
    flux_reg = 0;
}

ADR::ADR (Amr&            papa,
          ID              id,
          const Geometry& level_geom,
          const BoxArray& bl,
          Real            time)
    :
    AmrRegion(papa,id,level_geom,bl,time) 
{
    buildMetrics();

    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);

#ifdef DIFFUSION
      // diffusion is a static object, only alloc if not already there
      if (diffusion == 0) 
	diffusion = new Diffusion(master, master->finestLevel(), &phys_bc);

      std::cout << "CALLING INSTALL REGION AT ID " << m_id << std::endl;
      diffusion->install_region(m_id,this,volume,area);
#endif

}

ADR::~ADR () 
{
    delete flux_reg;
}

void
ADR::restart (Amr&     papa,
                 istream& is,
                 bool     bReadSpecial)
{
    AmrRegion::restart(papa,is,bReadSpecial);

    buildMetrics();

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);

#ifdef DIFFUSION
    if (level == 0) {
       BL_ASSERT(diffusion == 0);
       diffusion = new Diffusion(master,master->finestLevel(),&phys_bc);
    }
#endif
}

void
ADR::checkPoint(const std::string& dir,
                   std::ostream&  os,
                   VisMF::How     how,
                   bool dump_old_default)
{
  AmrRegion::checkPoint(dir, os, how, dump_old);
}

std::string
ADR::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
ADR::setPlotVariables ()
{
  AmrRegion::setPlotVariables();

  ParmParse pp("act");

  bool plot_X;

  if (pp.query("plot_X",plot_X))
  {
      if (plot_X)
      {
          //
	  // Get the number of species from the network model.
          //
	  BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);
          //
	  // Get the species names from the network model.
          //
	  for (int i = 0; i < NumSpec; i++)
          {
              int len = 20;
              Array<int> int_spec_names(len);
              //
              // This call return the actual length of each string in "len"
              //
              BL_FORT_PROC_CALL(GET_SPEC_NAMES, get_spec_names)
                  (int_spec_names.dataPtr(),&i,&len);
              char* spec_name = new char[len+1];
              for (int j = 0; j < len; j++) 
                  spec_name[j] = int_spec_names[j];
              spec_name[len] = '\0';
	      string spec_string = "X(";
              spec_string += spec_name;
              spec_string += ')';
	      master->addDerivePlotVar(spec_string);
              delete [] spec_name;
	  }
      }
  }
}

void
ADR::writePlotFile (const std::string& dir,
                       ostream&       os,
                       VisMF::How     how)
{
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (master->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if (master->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive++;
	}
    }

    int n_data_items = plot_var_map.size() + num_derive;

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for ( std::list<std::string>::iterator it = derive_names.begin();
	      it != derive_names.end(); ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
            os << rec->variableName(0) << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << master->cumTime() << '\n';
        int f_lev = master->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << master->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << master->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << master->regionSteps(ROOT_ID) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << master->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << master->regionSteps(ROOT_ID) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
	cnt++;
    }
    //
    // Cull data from derived variables.
    // 
    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::iterator it = derive_names.begin();
	     it != derive_names.end(); ++it) 
	{
	    MultiFab* derive_dat = derive(*it,cur_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,1,nGrow);
	    delete derive_dat;
	    cnt++;
	}
    }

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);
}

void
ADR::buildMetrics ()
{
    const int ngrd = grids.size();

    radius.resize(ngrd);

    const Real* dx = geom.CellSize();

    for (int i = 0; i < ngrd; i++)
    {
        const Box& b = grids[i];
        int ilo      = b.smallEnd(0)-radius_grow;
        int ihi      = b.bigEnd(0)+radius_grow;
        int len      = ihi - ilo + 1;

        radius[i].resize(len);

        Real* rad = radius[i].dataPtr();

        if (Geometry::IsCartesian())
        {
            for (int j = 0; j < len; j++)
            {
                rad[j] = 1.0;
            }
        }
        else
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

            const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dx[0];

            for (int j = 0; j < len; j++)
            {
                rad[j] = xlo + j*dx[0];
            }
        }
    }
    //
    // Build volume, face area and dLogArea arrays.
    // volume is not PArrayManaged, must manually delete.
    //
    volume.clear();
    //
    // area is not PArrayManaged, must manually delete.
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clear();
    }
    dLogArea[0].clear();
    geom.GetVolume(volume,grids,NUM_GROW);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        geom.GetFaceArea(area[dir],grids,dir,NUM_GROW);
    }
#if (BL_SPACEDIM <= 2)
    geom.GetDLogA(dLogArea[0],grids,0,NUM_GROW);
#endif
}

void
ADR::setTimeLevel (Real time,
                      Real dt_old,
                      Real dt_new)
{
    AmrRegion::setTimeLevel(time,dt_old,dt_new);
}

void
ADR::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    int ns          = NUM_STATE;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Initializing the data at level " << level << std::endl;

    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        RealBox    gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        BL_FORT_PROC_CALL(INITDATA,initdata)
	  (level, cur_time, lo, hi, ns,
	   BL_TO_FORTRAN(S_new[mfi]), dx,
	   gridloc.lo(), gridloc.hi());
    }

    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Done initializing the level " << level << " data " << std::endl;
}

void
ADR::init (AmrRegion &old)
{
    ADR* oldlev = (ADR*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = master->dtRegion(m_id);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    
    for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
          fpi.isValid();
          ++fpi)
    {
        S_new[fpi].copy(fpi());
    }
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
ADR::init ()
{
    Real dt        = master->dtRegion(m_id);
    Real cur_time  = parent_region->get_state_data(State_Type).curTime();
    Real prev_time = parent_region->get_state_data(State_Type).prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)master->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

Real
ADR::initial_time_step ()
{
    Real dummy_dt = 0.0;
    return (initial_dt > 0.0) ? initial_dt : init_shrink*est_time_step(dummy_dt);
}

Real
ADR::est_time_step (Real dt_old)
{
    if (fixed_dt > 0.0)
        return fixed_dt;

    // This is just a dummy value to start with 
    Real estdt  = 1.0e+20;

    const Real* dx    = geom.CellSize();
    const MultiFab& stateMF = get_new_data(State_Type);

    for (MFIter mfi(stateMF); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        Real dt = estdt;
	   BL_FORT_PROC_CALL(ESTDT,estdt)
            (BL_TO_FORTRAN(stateMF[mfi]),
             box.loVect(),box.hiVect(),dx,&dt);
  
        estdt = std::min(estdt,dt);
    }
    ParallelDescriptor::ReduceRealMin(estdt);
    estdt *= cfl;

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "ADR::est_time_step in region " << m_id << ":  estdt = " << estdt << '\n';

    return estdt;
}
void
ADR::computeNewDt (int                   finest_level,
                   int                   sub_cycle,
                   Tree<int>&            n_cycle,
                   const Array<IntVect>& ref_ratio,
                   Tree<Real>&           dt_reg_min,
                   Tree<Real>&           dt_region,
                   Real                  stop_time,
                   int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    ID id;
    TreeIterator<int> nc_it = n_cycle.getIterator();
    TreeIterator<Real> dt_it = dt_reg_min.getIterator(Prefix);
    for (; !dt_it.isFinished(); ++dt_it)
    {
        id = dt_it.getID();
        ADR& adv_level = get_region(id);
        dt_reg_min.setData(id,adv_level.est_time_step(dt_region.getData(id)));
    }

    if (fixed_dt <= 0.0)
    {
        bool sub_unchanged=true;
        if ((master->maxLevel() > 0) && (level == 0) &&
            (master->subcyclingMode() == "Optimal") && 
            (master->okToRegrid(m_id) || master->regionSteps(ROOT_ID) == 0) )
        {
            Tree<int> new_cycle(n_cycle);
            // The max allowable dt
            Tree<Real> dt_max(dt_reg_min);
            // find the maximum number of cycles allowed.
            std::list<int> structure = new_cycle.getStructure(ROOT_ID);
            Tree<int> cycle_max;
            Tree<Real> est_work;
            cycle_max.setRoot(1);
            est_work.setRoot(master->coarseRegion().estimateWork());
            cycle_max.buildFromStructure(ROOT_ID, structure);
            est_work.buildFromStructure(ROOT_ID, structure);
            // set non-root elements of est_work and cycle_max
            nc_it = cycle_max.getIterator(Prefix);
            for ( ++nc_it; !nc_it.isFinished(); ++nc_it)
            {
                id = nc_it.getID();
                int level = id.level();
                cycle_max.setData(id, master->MaxRefRatio(level-1));
                est_work.setData(id, master->getRegion(id).estimateWork());
            }

            // this value will be used only if the subcycling pattern is changed.
            dt_0 = master->computeOptimalSubcycling(new_cycle, dt_max, est_work, cycle_max);
            nc_it = cycle_max.getIterator(Prefix);
            for ( ; !nc_it.isFinished(); ++nc_it)
            {
                id = nc_it.getID();
                if (n_cycle.getData(id) != new_cycle.getData(id))
                {
                    sub_unchanged = false;
                    n_cycle.setData(id,new_cycle.getData(id));
                }
            }
            
        }
        
        if (sub_unchanged)
    //
        // Limit dt's by change_max * old dt
        //
        {
            dt_it = dt_reg_min.getIterator(Prefix);
            for (; !dt_it.isFinished(); ++dt_it)
            {
                id = dt_it.getID();
                if (verbose && ParallelDescriptor::IOProcessor())
                {
                    if (dt_reg_min.getData(id) > change_max*dt_region.getData(id))
                    {
                        cout << "ADR::compute_new_dt : limiting dt at level "
                             << id.size()-1 << " in region " << id << '\n';
                        cout << " ... new dt computed: " << dt_reg_min.getData(id)
                             << '\n';
                        cout << " ... but limiting to: "
                             << change_max * dt_region.getData(id) << " = " << change_max
                             << " * " << dt_region.getData(id) << '\n';
                    }
                }

                dt_reg_min.setData(id,std::min(dt_reg_min.getData(id), change_max * dt_region.getData(id)));
            }
            //
            // Find the max allowed Dt by finding the minmax over all levels
            //
            master->FindMaxDt(dt_0, n_cycle, dt_reg_min);
        }
        else
        {
            if (verbose && ParallelDescriptor::IOProcessor())
            {
                cout << "ADR: Changing subcycling pattern. New pattern:\n";
                cout << n_cycle.toString();
            }
        }
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    nc_it = n_cycle.getIterator(Prefix);
    dt_region.setRoot(dt_0);
    // Apply figure out all dt's from dt_0
    nc_it = n_cycle.getIterator(Prefix);
    for (++nc_it; !nc_it.isFinished(); ++nc_it)
    {
        id = nc_it.getID();
        dt_region.setData(id, dt_region.getData(nc_it.getParentID()) / n_cycle.getData(id));
    }
 
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "ADR: New dt:\n";
        cout << dt_region.toString();
    }
}

void
ADR::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Tree<int>&            n_cycle,
                       const Array<IntVect>& ref_ratio,
                       Tree<Real>&           dt_region,
                       Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    ID id;
    if (master->subcyclingMode() == "Optimal")
    {
        Tree<int> new_cycle(n_cycle);
        // find the maximum number of cycles allowed.
        std::list<int> structure = new_cycle.getStructure(ROOT_ID);
        Tree<int> cycle_max;
        Tree<Real> est_work;
        cycle_max.setRoot(1);
        est_work.setRoot(master->coarseRegion().estimateWork());
        cycle_max.buildFromStructure(ROOT_ID, structure);
        est_work.buildFromStructure(ROOT_ID, structure);
        // set non-root elements of est_work and cycle_max
        TreeIterator<int> nc_it = cycle_max.getIterator(Prefix);
        for ( ++nc_it; !nc_it.isFinished(); ++nc_it)
        {
            id = nc_it.getID();
            int level = id.level();
            cycle_max.setData(id, master->MaxRefRatio(level-1));
            est_work.setData(id, master->getRegion(id).estimateWork());
        }
        // The max allowable dt
        Tree<Real> dt_max;
        dt_max.buildFromStructure(ROOT_ID, structure);
        TreeIterator<Real> dt_it = dt_max.getIterator(Prefix);
        for (; !dt_it.isFinished(); ++dt_it)
        {
            id = dt_it.getID();
            dt_max.setData(id, get_region(id).initial_time_step());
        }
        dt_0 = master->computeOptimalSubcycling(new_cycle, dt_max, est_work, cycle_max);
        nc_it = n_cycle.getIterator(Prefix);
        for ( ; !nc_it.isFinished(); ++nc_it)
        {
            id = nc_it.getID();
            n_cycle.setData(id,new_cycle.getData(id));
            BL_ASSERT(n_cycle.getData(id) <= cycle_max.getData(id));
        }

        if (verbose && ParallelDescriptor::IOProcessor())
        {
            cout << "ADR: Initial subcycling pattern:\n";
            cout << n_cycle.toString();
        }
    }
    else
    {
        TreeIterator<Real> dt_it = dt_region.getIterator();
        for (; !dt_it.isFinished(); ++dt_it)
        {
            id = dt_it.getID();
            dt_region.setData(id,get_region(id).initial_time_step());
        }
        //
        // Find max coarse dt by minmaxing over levels
        //
        master->FindMaxDt(dt_0, n_cycle, dt_region);
    }
    
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }
    n_factor = 1;
    dt_region.setRoot(dt_0);
    // Apply figure out all dt's from dt_0
    TreeIterator<int> nc_it = n_cycle.getIterator(Prefix);
    for (++nc_it; !nc_it.isFinished(); ++nc_it)
    {
        id = nc_it.getID();
        dt_region.setData(id, dt_region.getData(nc_it.getParentID()) / n_cycle.getData(id));
    }
    
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        cout << "ADR: Initial dt:\n";
        cout << dt_region.toString();
    }
}

void
ADR::computeRestrictedDt (ID  base_region,
                          Tree<int>&  n_cycle,
                          Tree<Real>& dt_region)
{
    std::list<int> structure = n_cycle.getStructure(ROOT_ID);
    Tree<int> cycle_max;
    Tree<Real> dt_max;
    cycle_max.buildFromStructure(structure);
    dt_max.buildFromStructure(structure);
    TreeIterator<int> it = cycle_max.getIterator(base_region);
    ID id;
    // Setup cycle_max and dt_max
    for ( ; !it.isFinished(); ++it)
    {
        id = it.getID();
        if (id.size() == 1)
            cycle_max.setData(id, 1);
        else
            cycle_max.setData(id, master->MaxRefRatio(id.size()-2));
        dt_max.setData(id, get_region(id).est_time_step(dt_region.getData(id)));
        std::cout << "SETTING DT_MAX " << id << " " << dt_max.getData(id) << std::endl;
    }
    master->setRestrictedSubcycling(base_region, n_cycle, dt_region, dt_max, cycle_max);
}

void
ADR::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = master->finestLevel();

    if (do_reflux && level < finest_level) {

        MultiFab& S_new_crse = get_new_data(State_Type);

        //We must reflux if a region at the next finer level is subcycled relative to this level;
        //   otherwise the reflux was done as part of the multilevel advance
        reflux(-1);

        // We need to do this before anything else because refluxing changes the values of coarse cells
        //    underneath fine grids with the assumption they'll be over-written by averaging down
        if (level < finest_level)
           average_down();

        // This needs to be done after any changes to the state from refluxing.
        enforce_nonnegative_species(S_new_crse);

    }

    if (level < finest_level)
        average_down();

    if (level == 0)
    {
        int nstep = master->regionSteps(ROOT_ID);

        if ((sum_interval > 0) && (nstep%sum_interval == 0) )
            sum_integrated_quantities();
    }
}
void
ADR::post_restart ()
{
#ifdef DIFFUSION
      // diffusion is a static object, only alloc if not already there

      if (diffusion == 0)
	diffusion = new Diffusion(master, master->finestLevel(), &phys_bc);

      if (level == 0) 
      {
           RegionIterator it = master->getRegionIterator();
           for ( ; !it.isFinished(); ++it)
           {
               ADR& fine_region = get_region(it.getID());
               diffusion->install_region(it.getID(), *it, 
                                     fine_region.Volume(),fine_region.Area());
           }
      }
#endif
}

void
ADR::postCoarseTimeStep (Real cumtime)
{
    //
    // Only level 0 calls this routine.
    //
}

void
ADR::post_regrid (ID base_region, int new_finest)
{
}

void
ADR::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = master->finestLevel();

    for (int k = finest_level - 1; k >= 0; k--)
    {
        RegionIterator it = master->getRegionIterator(k);
        for ( ; !it.isFinished(); ++it)
        {
            get_region(it.getID()).average_down();
        }
    }

    if ( (sum_interval > 0) && (master->regionSteps(ROOT_ID)%sum_interval == 0) )
        sum_integrated_quantities();
}

int
ADR::okToContinue ()
{
    if (level > 0)
        return 1;

    int test = 1;
    if (master->dtRegion(ROOT_ID) < dt_cutoff)
        test = 0;

    return test;
}

void
ADR::getOldSource (Real old_time, Real dt, MultiFab&  ext_src)
{
   const Real* dx = geom.CellSize();

   MultiFab& S_old = get_old_data(State_Type);
   const int ncomp = S_old.nComp();

   ext_src.setVal(0.0,ext_src.nGrow());

   MultiFab levelArea[BL_SPACEDIM];
   for (int i = 0; i < BL_SPACEDIM ; i++)
       geom.GetFaceArea(levelArea[i],grids,i,NUM_GROW);

   for (FillPatchIterator Old_fpi(*this,S_old,4,old_time,State_Type,Density,ncomp);
                          Old_fpi.isValid();++Old_fpi)
   {
        const Box& bx = grids[Old_fpi.index()];
        BL_FORT_PROC_CALL(EXT_SRC,ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             dx,&old_time,&dt);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);
}

void
ADR::getNewSource (Real old_time, Real new_time, Real dt, MultiFab& ext_src)
{
   const Real* dx = geom.CellSize();

   MultiFab& S_old = get_old_data(State_Type);
   const int ncomp = S_old.nComp();

   ext_src.setVal(0.0,ext_src.nGrow());

   for (FillPatchIterator Old_fpi(*this,S_old,4,old_time,State_Type,Density,ncomp),
                          New_fpi(*this,S_old,4,new_time,State_Type,Density,ncomp);
                          Old_fpi.isValid() && New_fpi.isValid();++Old_fpi,++New_fpi)
   {
        const Box& bx = grids[Old_fpi.index()];
        BL_FORT_PROC_CALL(EXT_SRC,ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(  New_fpi()),
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             dx,&old_time,&dt);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);
}

void
ADR::reflux (int check_children)
{
    BL_ASSERT(level<master->finestLevel());

    const Real strt = ParallelDescriptor::second();

    BL_ASSERT(level<master->finestLevel());
 
    if (check_children != -1 && check_children != 0 && check_children != 1)
        BoxLib::Abort("Must have check_children = {-1,0,1} in reflux()");
 
    // Possible values for check_children
    // -1: exclude children with the same timestep (if called from post_timestep)
    // 0: no check
    // 1: exclude children with a different timestep (if called from multilevel advance)
 
    PList<AmrRegion> children;
    master->getRegions().getChildrenOfNode(m_id, children);
    for (PList<AmrRegion>::iterator it = children.begin(); it != children.end(); ++it)
    {
        ID c_id = (*it)->getID();
 
        // Note that if (strict_subcycling == 1) then every region advances
        //      independently so we reflux from every child regardless of its dt.
 
        // Skip this child if it wasn't part of the ml advance.
        if ( (check_children == 1) && (strict_subcycling == 0) &&
              master->dtRegion(m_id) != master->dtRegion(c_id))
            continue;
 
        // Skip this child if it was    part of the ml advance.
        if ( (check_children == -1) && (strict_subcycling == 0) &&
              master->dtRegion(m_id) == master->dtRegion(c_id))
            continue;

        std::cout << "Refluxing from region " << c_id << " to region " << m_id << std::endl;
        get_flux_reg(c_id).Reflux(get_new_data(State_Type),volume,1.0,0,0,NUM_STATE,geom);
    }

    if (verbose & 0)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "ADR::reflux() at level " << level << " : time = " << end << std::endl;
    }
}

void
ADR::average_down (bool check_children)
{
  if (level == master->finestLevel()) return;

  average_down(State_Type,check_children);
}

void
ADR::enforce_nonnegative_species (MultiFab& S_new)
{
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
       const Box bx = mfi.validbox();
       BL_FORT_PROC_CALL(ENFORCE_NONNEGATIVE_SPECIES,enforce_nonnegative_species)
           (BL_TO_FORTRAN(S_new[mfi]),bx.loVect(),bx.hiVect());
    }
}

void
ADR::average_down (int state_index, bool check_children)
{
    if (level == master->finestLevel()) return;
 
    MultiFab& S_crse = get_new_data(state_index);
 
    //iterate over chidlren
    RegionIterator c_it = master->getRegionIterator(m_id, level+1);
    for ( ; !c_it.isFinished(); ++c_it)
    {
        if ( check_children && strict_subcycling == 0 && 
             master->dtRegion(m_id) != master->dtRegion(c_it.getID()))
            continue;
        ADR& fine_region = get_region(c_it.getID());
        MultiFab& S_fine = fine_region.get_new_data(state_index);
        const int num_comps = S_fine.nComp();

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA(S_fine.boxArray().size());
 
        for (int i = 0; i < S_fine.boxArray().size(); ++i)
        {
            crse_S_fine_BA.set(i, BoxLib::coarsen(S_fine.boxArray()[i],
                                                  fine_ratio));
        }
 
        MultiFab crse_S_fine(crse_S_fine_BA, num_comps, 0);
        MultiFab crse_fvolume(crse_S_fine_BA,1,0);
        crse_fvolume.copy(volume);

        MultiFab& fvolume = fine_region.Volume();
        BL_ASSERT(S_crse.boxArray() == volume.boxArray());
        BL_ASSERT(fvolume.boxArray() == S_fine.boxArray());

        for (MFIter mfi(S_fine); mfi.isValid(); ++mfi)
        {
            const int        i        = mfi.index();
            const Box&       ovlp     = crse_S_fine_BA[i];
            FArrayBox&       crse_fab = crse_S_fine[i];
            const FArrayBox& crse_vol = crse_fvolume[i];
            const FArrayBox& fine_fab = S_fine[i];
            const FArrayBox& fine_vol = fvolume[i];
 
       	    BL_FORT_PROC_CALL(AVGDOWN,avgdown)
                (BL_TO_FORTRAN(crse_fab), num_comps,
                 BL_TO_FORTRAN(crse_vol),
                 BL_TO_FORTRAN(fine_fab),
                 BL_TO_FORTRAN(fine_vol),
                 ovlp.loVect(),ovlp.hiVect(),fine_ratio.getVect());
        }

        S_crse.copy(crse_S_fine);

     }
}

void
ADR::allocOldData ()
{
    for (int k = 0; k < NUM_STATE_TYPE; k++)
        state[k].allocOldData();
}

void
ADR::removeOldData()
{
    AmrRegion::removeOldData();
}

void
ADR::errorEst (TagBoxArray& tags,
                  int          clearval,
                  int          tagval,
                  Real         time,
                  int          n_error_buf,
                  int          ngrow)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();
    Array<int>  itags;

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

        for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
        {
            int         idx     = mfi.index();
            RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());
            itags               = tags[idx].tags();
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tags[idx].box().loVect();
            const int*  thi     = tags[idx].box().hiVect();
            const int*  lo      = mfi.validbox().loVect();
            const int*  hi      = mfi.validbox().hiVect();
            const Real* xlo     = gridloc.lo();
            Real*       dat     = (*mf)[mfi].dataPtr();
            const int*  dlo     = (*mf)[mfi].box().loVect();
            const int*  dhi     = (*mf)[mfi].box().hiVect();
            const int   ncomp   = (*mf)[mfi].nComp();

            err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
				  &clearval, dat, ARLIM(dlo), ARLIM(dhi),
				  lo,hi, &ncomp, domain_lo, domain_hi,
				  dx, xlo, prob_lo, &time, &level);
            //
            // Don't forget to set the tags in the TagBox.
            //
            if (allow_untagging == 1) 
            {
               tags[idx].tags_and_untags(itags);
            } else {
               tags[idx].tags(itags);
            }
        }

        delete mf;
    }
}

MultiFab*
ADR::derive (const std::string& name,
                Real           time,
                int            ngrow)
{
   return AmrRegion::derive(name,time,ngrow);
}

void
ADR::derive (const std::string& name,
                Real           time,
                MultiFab&      mf,
                int            dcomp)
{
    AmrRegion::derive(name,time,mf,dcomp);
}

Real
ADR::sumDerive (const std::string& name,
                   Real           time)
{
    Real sum     = 0.0;
    MultiFab* mf = derive(name, time, 0);

    BL_ASSERT(!(mf == 0));

    BoxArray ba_fine;

    if (level < master->finestLevel())
    {
        ba_fine = master->boxArray(m_id, level + 1);
        ba_fine.coarsen(fine_ratio);
    }

    std::vector< std::pair<int,Box> > isects;

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < master->finestLevel())
        {
            ba_fine.intersections(grids[mfi.index()],isects);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0, isects[ii].second, 0, fab.nComp());
            }
        }

        sum += fab.sum(0);
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

//
// Helper function for ADR::SyncInterp().
//

static
void
set_bc_new (int*            bc_new,
            int             n,
            int             src_comp,
            const int*      clo,
            const int*      chi,
            const int*      cdomlo,
            const int*      cdomhi,
            const BoxArray& cgrids,
            int**           bc_orig_qty)
            
{
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
        bc_new[bc_index]             = INT_DIR;
        bc_new[bc_index+BL_SPACEDIM] = INT_DIR;
 
        if (clo[dir] < cdomlo[dir] || chi[dir] > cdomhi[dir])
        {
            for (int crse = 0; crse < cgrids.size(); crse++)
            {
                const int* c_lo = cgrids[crse].loVect();
                const int* c_hi = cgrids[crse].hiVect();

                if (clo[dir] < cdomlo[dir] && c_lo[dir] == cdomlo[dir])
                    bc_new[bc_index] = bc_orig_qty[crse][bc_index];
                if (chi[dir] > cdomhi[dir] && c_hi[dir] == cdomhi[dir])
                    bc_new[bc_index+BL_SPACEDIM] = bc_orig_qty[crse][bc_index+BL_SPACEDIM]; 
            }
        }
    }
}

//
// Interpolate a cell-centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//
// This routine interpolates the num_comp components of CrseSync
// (starting at src_comp) and either increments or puts the result into
// the num_comp components of FineSync (starting at dest_comp)
// The components of bc_orig_qty corespond to the quantities of CrseSync.
//

void
ADR::SyncInterp (MultiFab&      CrseSync,
                    int            c_lev,
                    MultiFab&      FineSync,
                    int            f_lev,
                    IntVect&       ratio,
                    int            src_comp,
                    int            dest_comp,
                    int            num_comp,
                    int            increment,
                    Real           dt_clev, 
                    int**          bc_orig_qty,
                    SyncInterpType which_interp,
                    int            state_comp)
{
#if 0
    BL_ASSERT(which_interp >= 0 && which_interp <= 5);

    Interpolater* interpolater = 0;

    switch (which_interp)
    {
    case PC_T:           interpolater = &pc_interp;           break;
    case CellCons_T:     interpolater = &cell_cons_interp;    break;
    case CellConsLin_T:  interpolater = &lincc_interp;        break;
    case CellConsProt_T: interpolater = &protected_interp;    break;
    default:
        BoxLib::Abort("ADR::SyncInterp(): how did this happen");
    }

    ADR&   fine_level = getLevel(f_lev);
    const BoxArray& fgrids     = fine_level.boxArray();
    const Geometry& fgeom      = master->Geom(f_lev);
    const BoxArray& cgrids     = getLevel(c_lev).boxArray();
    const Geometry& cgeom      = master->Geom(c_lev);
    const Real*     dx_crse    = cgeom.CellSize();
    Box             cdomain    = BoxLib::coarsen(fgeom.Domain(),ratio);
    const int*      cdomlo     = cdomain.loVect();
    const int*      cdomhi     = cdomain.hiVect();
    int*            bc_new     = new int[2*BL_SPACEDIM*(src_comp+num_comp)];

    BoxArray cdataBA(fgrids.size());

    for (int i = 0; i < fgrids.size(); i++)
        cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
    //
    // Note: The boxes in cdataBA may NOT be disjoint !!!
    //
    MultiFab cdataMF(cdataBA,num_comp,0);

    cdataMF.setVal(0);

    cdataMF.copy(CrseSync, src_comp, 0, num_comp);
    //
    // Set physical boundary conditions in cdataMF.
    //
    // HACK HACK HACK -- for now to get it to compile
#if 1
    for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
        int         i       = mfi.index();
        RealBox     gridloc = RealBox(fine_level.boxArray()[i],
                                      fine_level.Geom().CellSize(),
                                      fine_level.Geom().ProbLo());
        FArrayBox&  cdata   = cdataMF[mfi];
        const int*  clo     = cdata.loVect();
        const int*  chi     = cdata.hiVect();
        const Real* xlo     = gridloc.lo();

        for (int n = 0; n < num_comp; n++)
        {
            set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);

            BL_FORT_PROC_CALL(FILCC,filcc)
                (BL_TO_FORTRAN(cdata),
                 cdomlo, cdomhi, dx_crse, xlo,
                 &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }
    }
#endif
    cgeom.FillPeriodicBoundary(cdataMF, 0, num_comp);
    //
    // Interpolate from cdataMF to fdata and update FineSync.
    // Note that FineSync and cdataMF will have the same distribution
    // since the length of their BoxArrays are equal.
    //
    FArrayBox    fdata;
    Array<BCRec> bc_interp(num_comp);

    MultiFab* fine_stateMF = 0;
    if (interpolater == &protected_interp)
    {
        fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));
    }

    for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
        int        i     = mfi.index();
        FArrayBox& cdata = cdataMF[mfi];
        const int* clo   = cdata.loVect();
        const int* chi   = cdata.hiVect();

        fdata.resize(fgrids[i], num_comp);
        //
        // Set the boundary condition array for interpolation.
        //
        for (int n = 0; n < num_comp; n++)
        {
            set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);
        }

        for (int n = 0; n < num_comp; n++)
        {
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
                bc_interp[n].setLo(dir,bc_new[bc_index]);
                bc_interp[n].setHi(dir,bc_new[bc_index+BL_SPACEDIM]);
            }
        }

        interpolater->interp(cdata,0,fdata,0,num_comp,fgrids[i],ratio,
                             cgeom,fgeom,bc_interp,src_comp,State_Type);

        if (increment)
        {
            fdata.mult(dt_clev);

            if (interpolater == &protected_interp)
            {
              cdata.mult(dt_clev);
              FArrayBox& fine_state = (*fine_stateMF)[i];
              interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
                                    num_comp,fgrids[i],ratio,
                                    cgeom,fgeom,bc_interp);
              Real dt_clev_inv = 1./dt_clev;
              cdata.mult(dt_clev_inv);
            }
            
            FineSync[i].plus(fdata,0,dest_comp,num_comp);
        }
        else
        {
            FineSync[i].copy(fdata,0,dest_comp,num_comp);
        }
    }

    delete [] bc_new;
#endif
}

void
ADR::adr_network_init ()
{
   BL_FORT_PROC_CALL(F_NETWORK_INIT,f_network_init) ();
}

void
ADR::restructure(std::list<int> structure)
{
#ifdef DIFFUSION
    diffusion->restructure(m_id, structure);
#endif
}
