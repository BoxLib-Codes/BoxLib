#include <winstd.H>
#include <algorithm>
#include <cstdio>
#include <list>
#include <iostream>
#include <sstream>

#include <Geometry.H>
#include <TagBox.H>
#include <Array.H>
#include <CoordSys.H>
#include <ParmParse.H>
#include <BoxDomain.H>
#include <Cluster.H>
#include <LevelBld.H>
#include <AmrRegion.H>
#include <PROB_AMR_F.H>
#include <Amr.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <DistributionMapping.H>
#include <FabSet.H>

#ifdef BL_USE_ARRAYVIEW
#include <DatasetClient.H>
#endif
//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif
//
// Static class members.  Set defaults in Initialize()!!!
//
std::list<std::string> Amr::state_plot_vars;
std::list<std::string> Amr::derive_plot_vars;
bool                   Amr::first_plotfile;
Array<BoxArray>        Amr::initial_ba;
Array<BoxArray>        Amr::regrid_ba;

namespace
{
    const std::string CheckPointVersion("Region_CheckPointVersion_1.0");

    bool initialized = false;
}

namespace
{
    //
    // These are all ParmParse'd in.  Set defaults in Initialize()!!!
    //
    int  plot_nfiles;
    int  mffile_nstreams;
    int  probinit_natonce;
    bool plot_files_output;
    int  checkpoint_nfiles;
    int  regrid_on_restart;
    int  use_efficient_regrid;
    bool refine_grid_layout;
    int  plotfile_on_restart;
    int  checkpoint_on_restart;
    bool checkpoint_files_output;
    int  compute_new_dt_on_regrid;
}

void
Amr::Initialize ()
{
    if (initialized) return;
    //
    // Set all defaults here!!!
    //
    Amr::first_plotfile      = true;
    plot_nfiles              = 64;
    mffile_nstreams          = 1;
    probinit_natonce         = 32;
    plot_files_output        = true;
    checkpoint_nfiles        = 64;
    regrid_on_restart        = 0;
    use_efficient_regrid     = 0;
    refine_grid_layout       = true;
    plotfile_on_restart      = 0;
    checkpoint_on_restart    = 0;
    checkpoint_files_output  = true;
    compute_new_dt_on_regrid = 0;

    BoxLib::ExecOnFinalize(Amr::Finalize);

    initialized = true;
}

void
Amr::Finalize ()
{
    Amr::state_plot_vars.clear();
    Amr::derive_plot_vars.clear();
    Amr::regrid_ba.clear();
    Amr::initial_ba.clear();

    initialized = false;
}

AmrRegion&
Amr::getLevel (int lev)
{
    if (multi_region)
        BoxLib::Abort("Called getLevel on multi-region level");
    ID id(lev+1);
    for (int i = 0; i <= lev; i++)
        id[i] = 0;
    return amr_regions.getData(id);
}

int
Amr::nCycle (int lev) const
{
    if (multi_region)
        BoxLib::Abort("Called nCycle on multi-region level");
    ID id(lev+1);
    for (int i = 0; i <= lev; i++)
        id[i] = 0;
    return n_cycle.getData(id);
}

int
Amr::levelSteps (int lev) const
{
    if (multi_region)
        BoxLib::Abort("Called levelSteps on multi-region level");
    ID id(lev+1);
    for (int i = 0; i <= lev; i++)
        id[i] = 0;
    return region_steps.getData(id);
}

int
Amr::regionSteps (ID region) const
{
    return region_steps.getData(region);
}

int
Amr::levelCount (int lev) const
{
    if (multi_region)
        BoxLib::Abort("Called levelCount on multi-region level");
    ID id(lev+1);
    for (int i = 0; i <= lev; i++)
        id[i] = 0;
    return region_count.getData(id);
}

bool Amr::Plot_Files_Output () { return plot_files_output; }

std::ostream&
Amr::DataLog (int i)
{
    return datalog[i];
}

int
Amr::NumDataLogs ()
{
    return datalog.size();
}

bool
Amr::RegridOnRestart () const
{
    return regrid_on_restart;
}

BoxArray
Amr::boxArray (int lev) const
{
    return boxArray(ROOT_ID, lev);
}

BoxArray
Amr::boxArray (ID base_region, int lev) const
{
    BoxList bl;
    PTreeConstIterator<AmrRegion> it = amr_regions.getConstIterator(base_region, lev);
    for ( ; !it.isFinished(); ++it)
    {
        bl.join((*it)->boxArray().boxList());
    }
    const BoxArray ba(bl);
    return ba;
}

BoxArray
Amr::boxArray (ID base_region, int lev, ExecutionTree* exec_tree) const
{
    BoxList bl;
    ExecutionTreeIterator it = exec_tree->getIterator(base_region, lev);
    for ( ; !it.isFinished(); ++it)
    {
        bl.join(amr_regions.getData(it.getID()).boxArray().boxList());
    }
    const BoxArray ba(bl);
    return ba;
}

BoxArray
Amr::boxArray (ID region_id) const
{
    return amr_regions.getData(region_id).boxArray();
}

long
Amr::cellCount (int lev)
{
    long cnt = 0;
    RegionIterator it = amr_regions.getIterator(lev);
    for ( ; !it.isFinished(); ++it)
    {
        cnt += (*it)->countCells();
    }
    return cnt;
}

long
Amr::cellCount (ID region_id)
{
    return amr_regions.getData(region_id).countCells();
}

int
Amr::numGrids (int lev)
{
    int cnt = 0;
    RegionIterator it = amr_regions.getIterator(lev);
    for ( ; !it.isFinished(); ++it)
    {
        cnt += (*it)->numGrids();
    }
    return cnt;
}

MultiFab*
Amr::derive (const std::string& name,
             Real               time,
             int                lev,
             int                ngrow)
{
    PArray<MultiFab> mfa(PArrayManage);
    int N = 0;
    RegionIterator it = amr_regions.getIterator(lev);
    for ( ; !it.isFinished(); ++it)
    {
        mfa.resize(N+1);
        mfa.set(N,(*it)->derive(name,time,ngrow));
        N++;
    }
    // Create a temp multifab as an aggregate of the region derived mfs.
    MultiFab temp_mf(mfa);
    // Create a deep copy of that mf to return. This is slightly more costly 
    // than we'd like but should work for the moment. \\\TODO: improve.
    int ncomp = temp_mf.nComp();
    MultiFab* mf = new MultiFab(temp_mf.boxArray(),ncomp, ngrow);
    MultiFab::Copy(*mf, temp_mf, ncomp, ncomp, ncomp, ngrow);
    return mf;
}

int
Amr::MaxRefRatio (int level) const
{
    int maxval = 0;
    for (int n = 0; n<BL_SPACEDIM; n++) 
        maxval = std::max(maxval,ref_ratio[level][n]);
    return maxval;
}

Amr::Amr ()
    :
    amr_regions(PTreeManage),
    datalog(PArrayManage)
{
    Initialize();
    //
    // Setup Geometry from ParmParse file.
    // May be needed for variableSetup or even getLevelBld.
    //
    Geometry::Setup();
    //
    // Determine physics class.
    //
    levelbld = getLevelBld();
    //
    // Global function that define state variables.
    //
    levelbld->variableSetUp();
    //
    // Set default values.
    //
    grid_eff               = 0.7;
    plot_int               = -1;
    n_proper               = 1;
    max_level              = -1;
    last_plotfile          = 0;
    last_checkpoint        = 0;
    record_run_info        = false;
    record_grid_info       = false;
    file_name_digits       = 5;
    record_run_info_terse  = false;
    bUserStopRequest       = false;

    int i;
    for (i = 0; i < BL_SPACEDIM; i++)
        isPeriodic[i] = false;

    ParmParse pp("amr");
    //
    // Check for command line flags.
    //
    verbose = 0;
    pp.query("v",verbose);

    pp.query("regrid_on_restart",regrid_on_restart);
    pp.query("use_efficient_regrid",use_efficient_regrid);
    pp.query("plotfile_on_restart",plotfile_on_restart);
    pp.query("checkpoint_on_restart",checkpoint_on_restart);

    pp.query("compute_new_dt_on_regrid",compute_new_dt_on_regrid);

    pp.query("refine_grid_layout", refine_grid_layout);

    pp.query("mffile_nstreams", mffile_nstreams);
    pp.query("probinit_natonce", probinit_natonce);

    probinit_natonce = std::max(1, std::min(ParallelDescriptor::NProcs(), probinit_natonce));

    pp.query("file_name_digits", file_name_digits);

    pp.query("initial_grid_file",initial_grids_file);
    pp.query("regrid_file"      , regrid_grids_file);

    if (pp.contains("run_log"))
    {
        std::string log_file_name;
        pp.get("run_log",log_file_name);
        setRecordRunInfo(log_file_name);
    }
    if (pp.contains("run_log_terse"))
    {
        std::string log_file_name;
        pp.get("run_log_terse",log_file_name);
        setRecordRunInfoTerse(log_file_name);
    }
    if (pp.contains("grid_log"))
    {
        std::string grid_file_name;
        pp.get("grid_log",grid_file_name);
        setRecordGridInfo(grid_file_name);
    }

    if (pp.contains("data_log"))
    {
      int num_datalogs = pp.countval("data_log");
      datalog.resize(num_datalogs);
      Array<std::string> data_file_names(num_datalogs);
      pp.queryarr("data_log",data_file_names,0,num_datalogs);
      for (int i = 0; i < num_datalogs; i++) 
        setRecordDataInfo(i,data_file_names[i]);
    }

    probin_file = "probin";  // Make "probin" the default

    if (pp.contains("probin_file"))
    {
        pp.get("probin_file",probin_file);
    }
    //
    // Restart or run from scratch?
    //
    pp.query("restart", restart_file);
    //
    // Read max_level and alloc memory for container objects.
    //
    pp.get("max_level", max_level);
    int nlev     = max_level+1;
    geom.resize(nlev);
    blocking_factor.resize(nlev);
    max_grid_size.resize(nlev);
    n_error_buf.resize(nlev);
    
    initRegions();
    //
    // Set bogus values.
    //
    for (i = 0; i < nlev; i++)
    {
        dt_region.setRoot(1.e200); // Something nonzero so old & new will differ
        region_steps.setRoot(0);
        region_count.setRoot(0);
        n_cycle.setRoot(0);
        dt_reg_min.setRoot(0.0);
        n_error_buf[i] = 1;
        blocking_factor[i] = 2;
        max_grid_size[i] = (BL_SPACEDIM == 2) ? 128 : 32;
    }

    if (max_level > 0) 
    {
       regrid_int.resize(max_level);
       for (i = 0; i < max_level; i++)
           regrid_int[i]  = 0;
    }

    ref_ratio.resize(max_level);
    for (i = 0; i < max_level; i++)
        ref_ratio[i] = IntVect::TheZeroVector();
    //
    // Read other amr specific values.
    //


    pp.query("n_proper",n_proper);
    pp.query("grid_eff",grid_eff);
    pp.queryarr("n_error_buf",n_error_buf,0,max_level);
    //
    // Read in the refinement ratio IntVects as integer BL_SPACEDIM-tuples.
    //
    if (max_level > 0)
    {
        const int nratios_vect = max_level*BL_SPACEDIM;

        Array<int> ratios_vect(nratios_vect);

        int got_vect = pp.queryarr("ref_ratio_vect",ratios_vect,0,nratios_vect);

        Array<int> ratios(max_level);

        const int got_int = pp.queryarr("ref_ratio",ratios,0,max_level);
   
        if (got_int == 1 && got_vect == 1 && ParallelDescriptor::IOProcessor())
        {
            BoxLib::Warning("Only input *either* ref_ratio or ref_ratio_vect");
        }
        else if (got_vect == 1)
        {
            int k = 0;
            for (i = 0; i < max_level; i++)
            {
                for (int n = 0; n < BL_SPACEDIM; n++,k++)
                    ref_ratio[i][n] = ratios_vect[k];
            }
        }
        else if (got_int == 1)
        {
            for (i = 0; i < max_level; i++)
            {
                for (int n = 0; n < BL_SPACEDIM; n++)
                    ref_ratio[i][n] = ratios[i];
            }
        }
        else
        {
            BoxLib::Error("Must input *either* ref_ratio or ref_ratio_vect");
        }
    }
    
    //
    // Setup plot and checkpoint controls.
    //
    initPltAndChk();
    
    //
    // Setup subcycling controls.
    //
    initSubcycle();

    //
    // Read in max_grid_size.  Use defaults if not explicitly defined.
    //
    int cnt = pp.countval("max_grid_size");

    if (cnt == 1)
    {
        //
        // Set all values to the single available value.
        //
        int the_max_grid_size = 0;

        pp.get("max_grid_size",the_max_grid_size);

        for (i = 0; i <= max_level; i++)
        {
            max_grid_size[i] = the_max_grid_size;
        }
    }
    else if (cnt > 1)
    {
        //
        // Otherwise we expect a vector of max_grid_size values.
        //
        pp.getarr("max_grid_size",max_grid_size,0,max_level+1);
    }
    //
    // Read in the blocking_factors.  Use defaults if not explicitly defined.
    //
    cnt = pp.countval("blocking_factor");

    if (cnt == 1)
    {
        //
        // Set all values to the single available value.
        //
        int the_blocking_factor = 0;

        pp.get("blocking_factor",the_blocking_factor);

        for (i = 0; i <= max_level; i++)
        {
            blocking_factor[i] = the_blocking_factor;
        }
    }
    else if (cnt > 1)
    {
        //
        // Otherwise we expect a vector of blocking factors.
        //
        pp.getarr("blocking_factor",blocking_factor,0,max_level+1);
    }
    //
    // Read in the regrid interval if max_level > 0.
    //
    if (max_level > 0) 
    {
       int numvals = pp.countval("regrid_int");
       if (numvals == 1)
       {
           //
           // Set all values to the single available value.
           //
           int the_regrid_int = 0;
           pp.query("regrid_int",the_regrid_int);
           for (i = 0; i < max_level; i++)
           {
               regrid_int[i] = the_regrid_int;
           }
       }
       else if (numvals < max_level)
       {
           BoxLib::Error("You did not specify enough values of regrid_int");
       }
       else 
       {
           //
           // Otherwise we expect a vector of max_level values
           //
           pp.queryarr("regrid_int",regrid_int,0,max_level);
       }
    }
    //
    // Read computational domain and set geometry.
    //
    Array<int> n_cell(BL_SPACEDIM);
    pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
    BL_ASSERT(n_cell.size() == BL_SPACEDIM);
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();
    Box index_domain(lo,hi);
    for (i = 0; i <= max_level; i++)
    {
        geom[i].define(index_domain);
        if (i < max_level)
            index_domain.refine(ref_ratio[i]);
    }
    //
    // Now define offset for CoordSys.
    //
    Real offset[BL_SPACEDIM];
    for (i = 0; i < BL_SPACEDIM; i++)
    {
        const Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
        offset[i]        = Geometry::ProbLo(i) + delta*lo[i];
    }
    CoordSys::SetOffset(offset);

    if (max_level > 0 && !initial_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' )
        std::ifstream is(initial_grids_file.c_str(),std::ios::in);

        if (!is.good())
            BoxLib::FileOpenFailed(initial_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        initial_ba.resize(in_finest);
        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                bx.refine(ref_ratio[lev-1]);
                bl.push_back(bx);
            }
            initial_ba[lev-1].define(bl);
        }
        is.close();
#undef STRIP
    }

    if (max_level > 0 && !regrid_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' )
        std::ifstream is(regrid_grids_file.c_str(),std::ios::in);

        if (!is.good())
            BoxLib::FileOpenFailed(regrid_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        regrid_ba.resize(in_finest);
        for (int lev = 1; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                 bx.refine(ref_ratio[lev-1]);
                 if (bx.longside() > max_grid_size[lev])
                 {
                     std::cout << "Grid " << bx << " too large" << '\n';
                     BoxLib::Error();
                 }
                 bl.push_back(bx);
            }
            regrid_ba[lev-1].define(bl);
        }
        is.close();
#undef STRIP
    }
}

bool
Amr::isStatePlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = state_plot_vars.begin(), End = state_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name)
            return true;
    }
    return false;
}

void
Amr::fillStatePlotVarList ()
{
    state_plot_vars.clear();
    const DescriptorList& desc_lst = AmrRegion::get_desc_lst();
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (desc_lst[typ].getType() == IndexType::TheCellType())
                state_plot_vars.push_back(desc_lst[typ].name(comp));
}

void
Amr::clearStatePlotVarList ()
{
    state_plot_vars.clear();
}

void
Amr::addStatePlotVar (const std::string& name)
{
    if (!isStatePlotVar(name))
        state_plot_vars.push_back(name);
}

void
Amr::deleteStatePlotVar (const std::string& name)
{
    if (isStatePlotVar(name))
        state_plot_vars.remove(name);
}

bool
Amr::isDerivePlotVar (const std::string& name)
{
    for (std::list<std::string>::const_iterator li = derive_plot_vars.begin(), End = derive_plot_vars.end();
         li != End;
         ++li)
    {
        if (*li == name)
            return true;
    }

    return false;
}

void 
Amr::fillDerivePlotVarList ()
{
    derive_plot_vars.clear();
    DeriveList& derive_lst = AmrRegion::get_derive_lst();
    std::list<DeriveRec>& dlist = derive_lst.dlist();
    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), End = dlist.end();
         it != End;
         ++it)
    {
        if (it->deriveType() == IndexType::TheCellType())
        {
            derive_plot_vars.push_back(it->name());
        }
    }
}

void
Amr::clearDerivePlotVarList ()
{
    derive_plot_vars.clear();
}

void
Amr::addDerivePlotVar (const std::string& name)
{
    if (!isDerivePlotVar(name))
        derive_plot_vars.push_back(name);
}

void
Amr::deleteDerivePlotVar (const std::string& name)
{
    if (isDerivePlotVar(name))
        derive_plot_vars.remove(name);
}

Amr::~Amr ()
{
    levelbld->variableCleanUp();
    Amr::Finalize();
}

void
Amr::setRecordGridInfo (const std::string& filename)
{
    record_grid_info = true;
    if (ParallelDescriptor::IOProcessor())
    {
        gridlog.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!gridlog.good())
            BoxLib::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier();
}

void
Amr::setRecordRunInfo (const std::string& filename)
{
    record_run_info = true;
    if (ParallelDescriptor::IOProcessor())
    {
        runlog.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!runlog.good())
            BoxLib::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier();
}

void
Amr::setRecordRunInfoTerse (const std::string& filename)
{
    record_run_info_terse = true;
    if (ParallelDescriptor::IOProcessor())
    {
        runlog_terse.open(filename.c_str(),std::ios::out|std::ios::app);
        if (!runlog_terse.good())
            BoxLib::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier();
}

void
Amr::setRecordDataInfo (int i, const std::string& filename)
{
    if (ParallelDescriptor::IOProcessor())
    {
        datalog.set(i,new std::ofstream);
        datalog[i].open(filename.c_str(),std::ios::out|std::ios::app);
        if (!datalog[i].good())
            BoxLib::FileOpenFailed(filename);
    }
    ParallelDescriptor::Barrier();
}

Real
Amr::dtLevel (int level) const
{
    if (multi_region)
        BoxLib::Abort("Called dtLevel on multi-region level");
    ID id(level+1);
    for (int i = 0; i <= level; i++)
        id[i] = 0;
    return dt_region.getData(id);
}

long
Amr::cellCount ()
{
    int cnt = 0;
    RegionIterator it = amr_regions.getIterator();
    for ( ; !it.isFinished(); ++it)
    {
        cnt += (*it)->countCells();
    }
    return cnt;
}

int
Amr::numGrids ()
{
    int cnt = 0;
    RegionIterator it = amr_regions.getIterator();
    for ( ; !it.isFinished(); ++it)
    {
        cnt += (*it)->numGrids();
    }
    return cnt;
}

int 
Amr::okToContinue ()
{
    int ok = true;
    RegionIterator it = amr_regions.getIterator();
    for ( ; ok && !it.isFinished(); ++it)
        ok = ok && (*it)->okToContinue();
    if(bUserStopRequest) {
      ok = false;
    }
    return ok;
}

void
Amr::writePlotFile ()
{
    if (!Plot_Files_Output()) return;

    BL_PROFILE("Amr::writePlotFile()");

    VisMF::SetNOutFiles(plot_nfiles);

    if (first_plotfile) 
    {
        first_plotfile = false;
        coarseRegion().setPlotVariables();
    }

    Real dPlotFileTime0 = ParallelDescriptor::second();


    const std::string pltfile = BoxLib::Concatenate(plot_file_root,
                                                    region_steps.getRoot(),file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "PLOTFILE: file = " << pltfile << '\n';

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "PLOTFILE: file = " << pltfile << '\n';

  BoxLib::StreamRetry sretry(pltfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  while(sretry.TryFileOutput()) {

    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(pltfile, 0755))
            BoxLib::CreateDirectoryFailed(pltfile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = pltfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec = 0;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        old_prec = HeaderFile.precision(15);
    }

    //
    // Get our aggregate regions for plotting.
    // This should eventually be upgraded 
    // when we have more sophisticated plt/chk logic.
    //
    PArray<AmrRegion> plot_levels(PArrayManage);
    aggregate_descendants(ROOT_ID, plot_levels);
    
    //
    // Write plotfiles for our plot regions.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        plot_levels[k].set_ancestors(plot_levels);
        plot_levels[k].writePlotFile(pltfile, HeaderFile);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);
        if (!HeaderFile.good())
            BoxLib::Error("Amr::writePlotFile() failed");
    }

    last_plotfile = region_steps.getRoot();

    if (verbose > 0)
    {
        const int IOProc        = ParallelDescriptor::IOProcessorNumber();
        Real      dPlotFileTime = ParallelDescriptor::second() - dPlotFileTime0;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(dPlotFileTime,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Write plotfile time = " << dPlotFileTime << "  seconds" << '\n';
#endif
    }
    ParallelDescriptor::Barrier();

  }  // end while

}

void
Amr::checkInput ()
{
    if (max_level < 0)
        BoxLib::Error("checkInput: max_level not set");
    //
    // Check that blocking_factor is a power of 2.
    //
    for (int i = 0; i < max_level; i++)
    {
        int k = blocking_factor[i];
        while ( k > 0 && (k%2 == 0) )
            k /= 2;
        if (k != 1)
            BoxLib::Error("Amr::checkInputs: blocking_factor not power of 2");
    }
    //
    // Check level dependent values.
    //
    int i;
    for (i = 0; i < max_level; i++)
    {
        if (MaxRefRatio(i) < 2 || MaxRefRatio(i) > 12)
            BoxLib::Error("checkInput bad ref_ratios");
    }
    const Box& domain = geom[0].Domain();
    if (!domain.ok())
        BoxLib::Error("level 0 domain bad or not set");
    //
    // Check that domain size is a multiple of blocking_factor[0].
    //
    for (i = 0; i < BL_SPACEDIM; i++)
    {
        int len = domain.length(i);
        if (len%blocking_factor[0] != 0)
            BoxLib::Error("domain size not divisible by blocking_factor");
    }
    //
    // Check that max_grid_size is even.
    //
    for (i = 0; i < max_level; i++)
    {
        if (max_grid_size[i]%2 != 0)
            BoxLib::Error("max_grid_size is not even");
    }

    //
    // Check that max_grid_size is a multiple of blocking_factor at every level.
    //
    for (i = 0; i < max_level; i++)
    {
        if (max_grid_size[i]%blocking_factor[i] != 0)
            BoxLib::Error("max_grid_size not divisible by blocking_factor");
    }

    if (!Geometry::ProbDomain().ok())
        BoxLib::Error("checkInput: bad physical problem size");

    if (max_level > 0) 
       if (regrid_int[0] <= 0)
          BoxLib::Error("checkinput: regrid_int not defined and max_level > 0");

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
       std::cout << "Successfully read inputs file ... " << '\n';
}

void
Amr::init (Real strt_time,
           Real stop_time)
{
    if (!restart_file.empty() && restart_file != "init")
    {
        restart(restart_file);
    }
    else
    {
        initialInit(strt_time,stop_time);
        checkPoint();
        if (plot_int > 0 || plot_per > 0)
            writePlotFile();
    }
#ifdef HAS_XGRAPH
    if (first_plotfile)
    {
        first_plotfile = false;
        coarseRegion().setPlotVariables();
    }
#endif
}

void
Amr::readProbinFile (int& init)
{
    BL_PROFILE("Amr::readProbinFile()");
    //
    // Populate integer array with name of probin file.
    //
    int probin_file_length = probin_file.length();

    Array<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
        probin_file_name[i] = probin_file[i];

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
       std::cout << "Starting to read probin ... " << std::endl;

    const int nAtOnce = probinit_natonce;
    const int MyProc  = ParallelDescriptor::MyProc();
    const int NProcs  = ParallelDescriptor::NProcs();
    const int NSets   = (NProcs + (nAtOnce - 1)) / nAtOnce;
    const int MySet   = MyProc/nAtOnce;

    Real piStart = 0, piEnd = 0, piStartAll = ParallelDescriptor::second();

    for (int iSet = 0; iSet < NSets; ++iSet)
    {
        if (MySet == iSet)
        {
            //
            // Call the pesky probin reader.
            //
            piStart = ParallelDescriptor::second();
            FORT_PROBINIT(&init,
                          probin_file_name.dataPtr(),
                          &probin_file_length,
                          Geometry::ProbLo(),
                          Geometry::ProbHi());
            piEnd = ParallelDescriptor::second();
            const int iBuff     = 0;
            const int wakeUpPID = (MyProc + nAtOnce);
            const int tag       = (MyProc % nAtOnce);
            if (wakeUpPID < NProcs)
                ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
        if (MySet == (iSet + 1))
        {
            //
            // Next set waits.
            //
            int iBuff;
            int waitForPID = (MyProc - nAtOnce);
            int tag        = (MyProc % nAtOnce);
            ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
        }
    }

    if (verbose > 1)
    {
        const int IOProc     = ParallelDescriptor::IOProcessorNumber();
        Real      piTotal    = piEnd - piStart;
        Real      piTotalAll = ParallelDescriptor::second() - piStartAll;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(piTotal,    IOProc);
        ParallelDescriptor::ReduceRealMax(piTotalAll, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "MFRead::: PROBINIT max time   = " << piTotal    << '\n';
            std::cout << "MFRead::: PROBINIT total time = " << piTotalAll << '\n';
        }
#endif
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "Successfully read probin file: \"" << probin_file << "\"\n";
}

void
Amr::initialInit (Real strt_time,
                  Real stop_time)
{
    checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    int init = true;

    if (!probin_file.empty()) {
        readProbinFile(init);
    }

#ifdef BL_SYNC_RANTABLES
    int iGet(0), iSet(1);
    const int iTableSize(64);
    Real *RanAmpl = new Real[iTableSize];
    Real *RanPhase = new Real[iTableSize];
    FORT_SYNC_RANTABLES(RanPhase, RanAmpl, &iGet);
    ParallelDescriptor::Bcast(RanPhase, iTableSize);
    ParallelDescriptor::Bcast(RanAmpl, iTableSize);
    FORT_SYNC_RANTABLES(RanPhase, RanAmpl, &iSet);
    delete [] RanAmpl;
    delete [] RanPhase;
#endif

    cumtime = strt_time;
    //
    // Define base level grids.
    //
    defBaseLevel(strt_time);
    //
    // Compute dt and set time levels of all grid data.
    //
    coarseRegion().computeInitialDt(finest_level,
                                    sub_cycle,
                                    n_cycle,
                                    ref_ratio,
                                    dt_region,
                                    stop_time);
    //
    // The following was added for multifluid.
    //
    /// I Don't feel like fixing this atm. Multifluid can use it if they really
    /// don't want to use the standard computeInitialDt model.
    //Real dt0   = dt_region.getRoot();
    //dt_min[0]  = dt_region.getRoot();
    //n_cycle[0] = 1;

    //std::list<int> structure = amr_regions.getStructure(ROOT_ID);
    
    //for (int lev = 1; lev <= max_level; lev++)
    //{
        //dt0           /= n_cycle[lev];
        //dt_level[lev]  = dt0;
        //dt_min[lev]    = dt_level[lev];
    //}

    if (max_level > 0)
        bldFineLevels(strt_time);
    
    //
    // Compute dt and set time levels of all grid data.
    // Only need to re-do this if we actually created new fine levels.
    //
    if (finest_level > 0)
        coarseRegion().computeInitialDt(finest_level,
                                        sub_cycle,
                                        n_cycle,
                                        ref_ratio,
                                        dt_region,
                                        stop_time);
    ID id;
    RegionIterator prit = amr_regions.getIterator(Prefix);
    for (; !prit.isFinished(); ++prit)
    {
        id = prit.getID();
        (*prit)->setTimeLevel(strt_time,dt_region.getData(id),dt_region.getData(id));
    }
    
    prit = amr_regions.getIterator(Prefix);
    for (; !prit.isFinished(); ++prit)
        (*prit)->post_regrid(ROOT_ID,finest_level);
    //
    // Perform any special post_initialization operations.
    //
    prit = amr_regions.getIterator( Prefix);
    for (; !prit.isFinished(); ++prit)
        (*prit)->post_init(stop_time);

    TreeIterator<int> iit = region_count.getIterator();
    for (; !iit.isFinished(); ++iit)
    {
        id = iit.getID();
        region_count.setData(id,0);
    }

    if (ParallelDescriptor::IOProcessor())
    {
       if (verbose > 1)
       {
           std::cout << "INITIAL GRIDS \n";
           printGridInfo(std::cout,0,finest_level);
       }
       else if (verbose > 0)
       { 
           std::cout << "INITIAL GRIDS \n";
           printGridSummary(std::cout,0,finest_level);
       }
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        gridlog << "INITIAL GRIDS \n";
        printGridInfo(gridlog,0,finest_level);
    }

#ifdef USE_STATIONDATA
    station.init(amr_regions, geom[finestLevel()].CellSize());
    station.findGrid(amr_regions,geom);
#endif
}

void
Amr::restart (const std::string& filename)
{
    BL_PROFILE("Amr::restart()");

    Real dRestartTime0 = ParallelDescriptor::second();

    VisMF::SetMFFileInStreams(mffile_nstreams);

    int i;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "restarting calculation from file: " << filename << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Init problem dependent data.
    //
    int init = false;

    readProbinFile(init);
    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Open the checkpoint header file for reading.
    //
    std::string File = filename;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    Array<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);
    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    
    int         spdim;
    bool        new_checkpoint_format = false;
    std::string first_line;

    std::getline(is,first_line);

    if (first_line == CheckPointVersion)
    {
        new_checkpoint_format = true;
        is >> spdim;
    }
    else
    {
        BoxLib::Error("Currently unable to read this checkpoint version");
        spdim = atoi(first_line.c_str());
    }

    if (spdim != BL_SPACEDIM)
    {
        std::cerr << "Amr::restart(): bad spacedim = " << spdim << '\n';
        BoxLib::Abort();
    }

    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    int num_regions;
    is >> num_regions;
    // Read the structure
    std::list<int> structure;
    int temp;
    for (int i = 0; i < num_regions; i++)
    {
        is >> temp;
        structure.push_back(temp);
    }
    
    Array<Box> inputs_domain(max_level+1);
    for (int lev = 0; lev <= max_level; lev++)
    {
       Box bx(geom[lev].Domain().smallEnd(),geom[lev].Domain().bigEnd());
       inputs_domain[lev] = bx;
    }

    if (max_level >= mx_lev) {

        for (i = 0; i <= mx_lev; i++) is >> geom[i];
        for (i = 0; i <  mx_lev; i++) is >> ref_ratio[i];
        
    }
    else
    {
        if (ParallelDescriptor::IOProcessor())
            BoxLib::Warning("Amr::restart(): max_level is lower than before");

        // These are just used to hold the extra stuff we have to read in.
        Geometry   geom_dummy;
        //Real       real_dummy;
        //int         int_dummy;
        IntVect intvect_dummy;

        for (i = 0          ; i <= max_level; i++) is >> geom[i];
        for (i = max_level+1; i <= mx_lev   ; i++) is >> geom_dummy;

        for (i = 0        ; i <  max_level; i++) is >> ref_ratio[i];
        for (i = max_level; i <  mx_lev   ; i++) is >> intvect_dummy;
    }
    
    
    dt_region.buildFromStructure(structure);
    for(TreeIterator<Real> it = dt_region.getIterator(Prefix); !it.isFinished(); ++it)
        is >> *it;
        
    finest_level = dt_region.getFinestLevel();
    BL_ASSERT(finest_level <= max_level);

    dt_reg_min.buildFromStructure(structure);
    for(TreeIterator<Real> it = dt_reg_min.getIterator(Prefix); !it.isFinished(); ++it)
        is >> *it;

    Tree<int>  n_cycle_in;
    n_cycle_in.buildFromStructure(structure);
    for(TreeIterator<int> it = n_cycle_in.getIterator(Prefix); !it.isFinished(); ++it)
        is >> *it;
        
    // n_cycle has not been properly set yet, so we do that now.
    n_cycle.buildFromStructure(structure);
    if (subcycling_mode == "Optimal" || subcycling_mode == "Auto")
    {
        // if subcycling mode is Optimal, n_cycle is set dynamically.
        // We'll initialize it to be Auto subcycling.
        n_cycle.setRoot(1);
        TreeIterator<int> it = n_cycle.getIterator(Prefix);
        for (++it; !it.isFinished(); ++it)
        {
            (*it) = MaxRefRatio(it.getLevel()-1);
        } 
    }
    else if (subcycling_mode == "Manual")
    {
        n_cycle.setRoot(1);
        TreeIterator<int> it = n_cycle.getIterator(Prefix);
        for (++it; !it.isFinished(); ++it)
        {
            (*it) = manual_n_cycle[it.getLevel()];
        } 
    }
    else if (subcycling_mode == "None")
    {
        n_cycle.setRoot(1);
        TreeIterator<int> it = n_cycle.getIterator(Prefix);
        for (++it; !it.isFinished(); ++it)
        {
            (*it) = 1;
        } 
    }
    else
    {
        BoxLib::Abort("Unknown Subcycling Mode\n");
    }
    bool any_changed = false;

    if (subcycling_mode == "Manual" || subcycling_mode == "Auto" || subcycling_mode == "None")
    {
        for(TreeIterator<int> it = n_cycle_in.getIterator(Prefix); !it.isFinished(); ++it)
        {
            ID id =it.getID();
            if (n_cycle.getData(id) != n_cycle_in.getData(id))
            {
                any_changed = true;
                if (verbose > 0 && ParallelDescriptor::IOProcessor())
                    std::cout << "Warning: n_cycle has changed at region " << id << 
                                " from " << n_cycle_in.getData(id) << " to " << n_cycle.getData(id) << std::endl;;
            }
        }
    }
    //If we change n_cycle then force a full regrid from level 0 up
    if (max_level > 0 && any_changed)
    {
       region_count.setRoot(regrid_int[0]);
       if ((verbose > 0) && ParallelDescriptor::IOProcessor())
           std::cout << "Warning: This forces a full regrid " << std::endl;
    }

    region_steps.buildFromStructure(structure);
    for(TreeIterator<int> it = region_steps.getIterator(Prefix); !it.isFinished(); ++it)
        is >> *it;
    region_count.buildFromStructure(structure);
    for(TreeIterator<int> it = region_count.getIterator(Prefix); !it.isFinished(); ++it)
        is >> *it;

    if (regrid_on_restart && max_level > 0)
        region_count.setRoot(regrid_int[0]);

    checkInput();
    //
    // Read levels.
    //
    amr_regions.buildFromStructure(structure);
    
    for (RegionIterator it = getRegionIterator(Prefix); !it.isFinished(); ++it)
    { 
        amr_regions.setData(it.getID(), (*levelbld)());
        amr_regions.getData(it.getID()).restart(*this, is);
    }
    //
    // Build any additional data structures.
    // 
    for (RegionIterator it = getRegionIterator(Prefix); !it.isFinished(); ++it)
    { 
        (*it)->post_restart();
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
       Box restart_domain(geom[lev].Domain());
       if (! (inputs_domain[lev] == restart_domain) )
       {
          if (ParallelDescriptor::IOProcessor())
          {
             std::cout << "Problem at level " << lev << '\n';
             std::cout << "Domain according to     inputs file is " <<  inputs_domain[lev] << '\n';
             std::cout << "Domain according to checkpoint file is " << restart_domain      << '\n';
             std::cout << "Amr::restart() failed -- box from inputs file does not equal box from restart file" << std::endl;
          }
          BoxLib::Abort();
       }
    }

#ifdef USE_STATIONDATA
    station.init(amr_regions, geom[finestLevel()].CellSize());
    station.findGrid(amr_regions,geom);
#endif

    if (verbose > 0)
    {
        Real dRestartTime = ParallelDescriptor::second() - dRestartTime0;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Restart time = " << dRestartTime << " seconds." << '\n';
#endif
    }
}

void
Amr::checkPoint ()
{
    if (!checkpoint_files_output) return;

    BL_PROFILE("Amr::checkPoint()");

    VisMF::SetNOutFiles(checkpoint_nfiles);
    //
    // In checkpoint files always write out FABs in NATIVE format.
    //
    FABio::Format thePrevFormat = FArrayBox::getFormat();

    FArrayBox::setFormat(FABio::FAB_NATIVE);

    Real dCheckPointTime0 = ParallelDescriptor::second();

    const std::string ckfile = BoxLib::Concatenate(check_file_root,region_steps.getRoot(),file_name_digits);

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "CHECKPOINT: file = " << ckfile << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "CHECKPOINT: file = " << ckfile << '\n';


  BoxLib::StreamRetry sretry(ckfile, abort_on_stream_retry_failure,
                             stream_max_tries);

  while(sretry.TryFileOutput()) {

    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(ckfile, 0755))
            BoxLib::CreateDirectoryFailed(ckfile);
    //Create the Level directories for the Region data
    for (int i = 0; i < finest_level; i++)
    {
        std::string levfile = BoxLib::Concatenate(ckfile+"/Level_", i, 1);
        if (!BoxLib::UtilCreateDirectory(levfile, 0755))
            BoxLib::CreateDirectoryFailed(levfile);
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = ckfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec = 0, i;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);

        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);

        old_prec = HeaderFile.precision(17);

        HeaderFile << CheckPointVersion << '\n'
                   << BL_SPACEDIM       << '\n'
                   << cumtime           << '\n'
                   << max_level         << '\n';
        // We output the region tree structure in list form. 
        // (See the tree classes for details)
        // First the number of regions
        std::list<int> structure = amr_regions.getStructure();
        HeaderFile << structure.size() << '\n';
        for (std::list<int>::iterator it = structure.begin(); it!= structure.end(); ++it)
        {
            HeaderFile << *it << ' ';
        }
        HeaderFile << '\n';
        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) HeaderFile << geom[i]        << ' ';
        HeaderFile << '\n';
        for (i = 0; i < max_level; i++)  HeaderFile << ref_ratio[i]   << ' ';
        HeaderFile << '\n';
        TreeIterator<Real> rit = dt_region.getIterator(Prefix);
        for (; !rit.isFinished(); ++rit)
            HeaderFile << *rit  << ' ';
        HeaderFile << '\n';
        for (rit = dt_reg_min.getIterator(Prefix); !rit.isFinished(); ++rit)
            HeaderFile << *rit  << ' ';
        HeaderFile << '\n';
        TreeIterator<int> iit = n_cycle.getIterator(Prefix);
        for ( ; !iit.isFinished(); ++iit)
            HeaderFile << *iit  << ' ';
        HeaderFile << '\n';
        for (iit = region_steps.getIterator(Prefix); !iit.isFinished(); ++iit)
            HeaderFile << *iit  << ' ';
        HeaderFile << '\n';
        for (iit = region_count.getIterator(Prefix); !iit.isFinished(); ++iit)
            HeaderFile << *iit  << ' ';
        HeaderFile << '\n';
    }

    RegionIterator it = amr_regions.getIterator(Prefix);
    for ( ; !it.isFinished(); ++it)
    {
        (*it)->checkPoint(ckfile, HeaderFile);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile.precision(old_prec);

        if (!HeaderFile.good())
            BoxLib::Error("Amr::checkpoint() failed");
    }

    last_checkpoint = region_steps.getRoot();

#ifdef USE_SLABSTAT
    //
    // Dump out any SlabStats MultiFabs.
    //
    AmrRegion::get_slabstat_lst().checkPoint(getRegions(), region_steps.getRoot());
#endif
    //
    // Don't forget to reset FAB format.
    //
    FArrayBox::setFormat(thePrevFormat);

    if (verbose > 0)
    {
        Real dCheckPointTime = ParallelDescriptor::second() - dCheckPointTime0;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(dCheckPointTime,
	                            ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "checkPoint() time = " << dCheckPointTime << " secs." << '\n';
#endif
    }
    ParallelDescriptor::Barrier();

  }  // end while

}

void
Amr::RegridOnly (Real time)
{
    BL_ASSERT(regrid_on_restart == 1);

    int lev_top = std::min(finest_level, max_level-1);

    for (int i = 0; i <= lev_top; i++)
       regrid(i,time);

    if (plotfile_on_restart)
	writePlotFile();

    if (checkpoint_on_restart)
       checkPoint();

}

void
Amr::timeStep (AmrRegion& base_region,
               Real time,
               int  iteration,
               int  niter,
               Real stop_time)
{
    BL_PROFILE("Amr::timeStep()");

    int level = base_region.Level();
    ID base_id = base_region.getID();

    //
    // Allow regridding of level 0 calculation on restart.
    //
    if (max_level == 0 && regrid_on_restart)
    {
        regrid_on_restart = 0;
        //
        // Coarsening before we split the grids ensures that each resulting
        // grid will have an even number of cells in each direction.
        //
        BoxArray lev0(BoxLib::coarsen(geom[0].Domain(),2));
        //
        // Now split up into list of grids within max_grid_size[0] limit.
        //
        lev0.maxSize(max_grid_size[0]/2);
        //
        // Now refine these boxes back to level 0.
        //
        lev0.refine(2);

        //
        // If use_efficient_regrid flag is set, then test to see whether we in fact 
        //    have just changed the level 0 grids. If not, then don't do anything more here.
        //
        if ( !( (use_efficient_regrid == 1) && (lev0 == coarseRegion().boxArray()) ) ) 
        {
            //
            // Construct skeleton of new level.
            //
            AmrRegion* a = (*levelbld)(*this,ROOT_ID,geom[0],lev0,cumtime);

            a->init(coarseRegion());
            amr_regions.clearData(ROOT_ID);
            amr_regions.setRoot( a);
            // Updated everyone's ancestors
            for(RegionIterator it = getRegionIterator(); !it.isFinished(); ++it)
            {
                (*it)->get_ancestors().set(0,a);
            }

            coarseRegion().post_regrid(ROOT_ID,0);

            if (ParallelDescriptor::IOProcessor())
            {
               if (verbose > 1)
               {
                  printGridInfo(std::cout,0,finest_level);
               }
               else if (verbose > 0)
               {
                  printGridSummary(std::cout,0,finest_level);
               }
            }

            if (record_grid_info && ParallelDescriptor::IOProcessor())
                printGridInfo(gridlog,0,finest_level);
        }
        else
        {
            if (verbose > 0 && ParallelDescriptor::IOProcessor())
                std::cout << "Regridding at level 0 but grids unchanged " << std::endl;
        }
    }
    else
    {
        //See if any subregions want to regrid in prefix fashion
        RegionIterator it = amr_regions.getIterator(base_id,Prefix);
        for ( ; !it.isFinished(); ++it)
        {
            //const int old_finest = finest_level;
            int subiteration = it.getID() == base_id ? iteration : 1;
            if (okToRegrid(it.getID(), subiteration))
            {
                regrid(it.getID(),time);

                //
                // Compute new dt after regrid if at level 0 and compute_new_dt_on_regrid.
                //
                if ( compute_new_dt_on_regrid && ((*it)->Level() == 0) )
                {
                    int post_regrid_flag = 1;
                    coarseRegion().computeNewDt(finest_level,
                                              sub_cycle,
                                              n_cycle,
                                              ref_ratio,
                                              dt_reg_min,
                                              dt_region,
                                              stop_time, 
                                              post_regrid_flag);
                }
            }
        }
    }
    //
    // Check to see if should write plotfile.
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && !(restart_file.empty()) )
    {
	plotfile_on_restart = 0;
	writePlotFile();
    }
    //
    // Advance grids at this level.
    //
    
    Real my_dt = dt_region.getData(base_id);
    
    // We now put the print statements in the advance routines.
    Real dt_new = base_region.advance(time,my_dt,iteration,niter);
    
    dt_reg_min.setData(base_id, iteration == 1 ? dt_new : std::min(dt_reg_min.getData(base_id),dt_new));

    region_steps.setData(base_id, region_steps.getData(base_id) + 1);
    region_count.setData(base_id, region_count.getData(base_id) + 1);

#ifdef USE_STATIONDATA
    station.report(time+my_dt,base_region);
#endif

#ifdef USE_SLABSTAT
    AmrRegion::get_slabstat_lst().update(base_region,time,my_dt);
#endif
    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
        RegionList children;
        amr_regions.getChildrenOfNode(base_region.getID(), children);
        for (RegionList::iterator child = children.begin(); child != children.end(); ++child)
        {
            ID c_id = (*child)->getID();
            if (sub_cycle)
            {
                const int ncycle = n_cycle.getData(c_id);
                for (int i = 1; i <= ncycle; i++)
                {
                    timeStep(**child,time+(i-1)*dt_region.getData(c_id),i,ncycle,stop_time);
                }
            }
            else
            {
                timeStep(**child,time,1,1,stop_time);
            }
        }
    }

    base_region.post_timestep(iteration);
}

void
Amr::coarseTimeStep (Real stop_time)
{
    BL_PROFILE("Amr::coarseTimeStep()");

    const Real run_strt = ParallelDescriptor::second() ;
    //
    // Compute new dt.
    //
    if (regionSteps(ROOT_ID) > 0)
    {
        int post_regrid_flag = 0;
        coarseRegion().computeNewDt(finest_level,
                                  sub_cycle,
                                  n_cycle,
                                  ref_ratio,
                                  dt_reg_min,
                                  dt_region,
                                  stop_time,
                                  post_regrid_flag);
    }
    else
    {
        coarseRegion().computeInitialDt(finest_level,
                                      sub_cycle,
                                      n_cycle,
                                      ref_ratio,
                                      dt_region,
                                      stop_time);
    }
    timeStep(coarseRegion(),cumtime,1,1,stop_time);

    cumtime += dt_region.getRoot();

    coarseRegion().postCoarseTimeStep(cumtime);

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

        long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
        long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;

        ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
        //
        // Reset to zero to calculate high-water-mark for next timestep.
        //
        BoxLib::total_bytes_allocated_in_fabs_hwm = 0;

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "\nFAB byte spread across MPI nodes for timestep: ["
                      << min_fab_bytes
                      << " ... "
                      << max_fab_bytes
                      << "]\n";
        }
    }

    BL_PROFILE_ADD_STEP(level_steps[0]);
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nSTEP = "
                 << region_steps.getRoot()
                  << " TIME = "
                  << cumtime
                  << " DT = "
                  << dt_region.getRoot()
                  << '\n'
                  << std::endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "
               << region_steps.getRoot()
               << " TIME = "
               << cumtime
               << " DT = "
               << dt_region.getRoot()
               << '\n';
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << region_steps.getRoot() << " " << cumtime << " " << dt_region.getRoot() << '\n';

    int check_test = 0;
    if (check_per > 0.0)
    {
      const int num_per_old = cumtime / check_per;
      const int num_per_new = (cumtime+dt_region.getRoot()) / check_per;

      if (num_per_old != num_per_new)
	{
	 check_test = 1;
	}
    }

    int to_stop       = 0;    
    int to_checkpoint = 0;
    if (ParallelDescriptor::IOProcessor())
    {
        FILE *fp;
        if ((fp=fopen("dump_and_continue","r")) != 0)
        {
            remove("dump_and_continue");
            to_checkpoint = 1;
            fclose(fp);
        }
        else if ((fp=fopen("stop_run","r")) != 0)
        {
            remove("stop_run");
            to_stop = 1;
            fclose(fp);
        }
        else if ((fp=fopen("dump_and_stop","r")) != 0)
        {
            remove("dump_and_stop");
            to_checkpoint = 1;
            to_stop = 1;
            fclose(fp);
        }
    }
    ParallelDescriptor::Bcast(&to_checkpoint, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&to_stop,       1, ParallelDescriptor::IOProcessorNumber());

    if(to_stop == 1 && to_checkpoint == 0) {  // prevent main from writing files
      last_checkpoint = region_steps.getRoot();
      last_plotfile   = region_steps.getRoot();
    }
    if ((check_int > 0 && region_steps.getRoot() % check_int == 0) || check_test == 1
	|| to_checkpoint)
        checkPoint();

    if (writePlotNow() || to_checkpoint)
    {
        writePlotFile();
    }

    bUserStopRequest = to_stop;
    if (to_stop)
    {
        ParallelDescriptor::Barrier();
        if(ParallelDescriptor::IOProcessor()) {
          if (to_checkpoint)
          {
            std::cerr << "Stopped by user w/ checkpoint" << std::endl;
          }
          else
          {
            std::cerr << "Stopped by user w/o checkpoint" << std::endl;
          }
	}
    }
}

bool
Amr::writePlotNow()
{
    int plot_test = 0;
    if (plot_per > 0.0)
    {
#ifdef BL_USE_NEWPLOTPER
      Real rN(0.0);
      Real rR = modf(cumtime/plot_per, &rN);
      if (rR < (dt_region.getRoot()*0.001))
#else
      const int num_per_old = (cumtime-dt_region.getRoot()) / plot_per;
      const int num_per_new = (cumtime            ) / plot_per;

      if (num_per_old != num_per_new)
#endif
	{
	  plot_test = 1;
	}
    }

    return ( (plot_int > 0 && region_steps.getRoot() % plot_int == 0) || 
              plot_test == 1 || coarseRegion().writePlotNow() );
} 

void
Amr::defBaseLevel (Real strt_time)
{
    //
    // Check that base domain has even number of zones in all directions.
    //
    const Box& domain = geom[0].Domain();
    IntVect d_length  = domain.size();
    const int* d_len  = d_length.getVect();

    for (int idir = 0; idir < BL_SPACEDIM; idir++)
        if (d_len[idir]%2 != 0)
            BoxLib::Error("defBaseLevel: must have even number of cells");
    //
    // Coarsening before we split the grids ensures that each resulting
    // grid will have an even number of cells in each direction.
    //
    BoxArray lev0(1);

    lev0.set(0,BoxLib::coarsen(domain,2));
    //
    // Now split up into list of grids within max_grid_size[0] limit.
    //
    lev0.maxSize(max_grid_size[0]/2);
    //
    // Now refine these boxes back to level 0.
    //
    lev0.refine(2);

    //
    // If (refine_grid_layout == 1) and (Nprocs > ngrids) then break up the 
    //    grids into smaller chunks
    //
    Array<BoxArray> new_grids(1);
    new_grids[0] = lev0;
    impose_refine_grid_layout(0,0,new_grids);
    lev0 = new_grids[0];

    //
    // Now build level 0 grids.
    //
    amr_regions.setRoot((*levelbld)(*this,ROOT_ID,geom[0],lev0,strt_time));

    lev0.clear();
    //
    // Now init level 0 grids with data.
    //
    coarseRegion().initData();
}


void 
Amr::restructure (ID base_region, std::list<int> structure, bool do_regions)
{
    // Note: This clear might have been problematic save for the fact
    // that regrid has already evicted the children
    if (do_regions)
    {
        amr_regions.clearChildrenOfNode(base_region); 
        amr_regions.buildFromStructure(base_region, structure);
    }
    
    n_cycle.clearChildrenOfNode(base_region);
    n_cycle.buildFromStructure(base_region, structure);
    
    TreeIterator<int> iit = n_cycle.getIterator(base_region, Prefix);
    for (++iit; !iit.isFinished(); ++iit)
    {
        if (subcycling_mode == "None")
            (*iit) = 1;
        else if (subcycling_mode == "Auto" || subcycling_mode == "Optimal")
        {
            (*iit) = MaxRefRatio(iit.getLevel()-1);
        }
        else if (subcycling_mode == "Manual")
        {
            (*iit) = manual_n_cycle[iit.getLevel()];
        }
        else
        {
            BoxLib::Abort("Invalid Subcycling mode");
        }
    }
    
    region_count.clearChildrenOfNode(base_region);
    region_count.buildFromStructure(base_region, structure);
    iit = region_count.getIterator(base_region, Prefix);
    for (; !iit.isFinished(); ++iit)
    {
        (*iit) = 0;
    }
    
    region_steps.clearChildrenOfNode(base_region);
    region_steps.buildFromStructure(base_region, structure);
    iit = region_steps.getIterator(base_region, Prefix);
    for (++iit; !iit.isFinished(); ++iit)
    {
        (*iit) = 0;
    }
    
    dt_region.clearChildrenOfNode(base_region);
    dt_region.buildFromStructure(base_region, structure);
    Real dt_base = dt_region.getData(base_region);
    TreeIterator<Real> rit = dt_region.getIterator(base_region,Prefix);
    for ( ++rit; !rit.isFinished(); ++rit)
    {
        *rit = dt_region.getData(rit.getParentID())/n_cycle.getData(rit.getID());
    }
    
    dt_reg_min.clearChildrenOfNode(base_region);
    dt_reg_min.buildFromStructure(base_region, structure);
    dt_base = dt_reg_min.getData(base_region);
    rit = dt_reg_min.getIterator(base_region, Prefix);
    for (  ++rit; !rit.isFinished(); ++rit)
    {
        *rit = dt_reg_min.getData(rit.getParentID())/n_cycle.getData(rit.getID());
    }
    
    getRegion(base_region).restructure(structure);
}

void
Amr::regrid (int  lbase,
             Real time,
             bool initial)
{

    if (multi_region)
        BoxLib::Abort("Calling lbase regrid in a multiregion simulation\n");
    ID region_id(lbase+1);
    regrid(region_id, time, initial);
}

void
Amr::regrid (ID base_id,
             Real time,
             bool initial)
{
    BL_PROFILE("Amr::regrid()");

    AmrRegion& base_region = getRegion(base_id);
    int lbase = base_region.Level();

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "REGRID: at base region " << base_id << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "REGRID: at base region " << base_id << '\n';

    // This list tracks all regions that need to call post_regrid.
    RegionList touched_regions;
    
    //
    // Create the level hierarchy with the data we can see.
    //
    PArray<AmrRegion> active_levels;
    active_levels.resize(finest_level+1);
    Array<DistributionMapping> dms(finest_level+1);
    for (int i = 0; i <= lbase; i++)
    {
        ID ancestor_id(base_id);
        ancestor_id.resize(i+1);
        active_levels.set(i,&getRegion(ancestor_id));
        touched_regions.push_back(&getRegion(ancestor_id));
    }
    
    PArray<AmrRegion> descendants(PArrayManage);
    aggregate_descendants(base_id, descendants);
    for (int i = lbase+1; i <= finest_level; i++)
    {
        active_levels.set(i,&descendants[i-lbase]);
        active_levels[i].set_ancestors(active_levels);
    }
    for(int i = 0; i < finest_level; i++)
        dms[i] = active_levels[i].get_distribution_map();

    //
    // Compute positions of new grids.
    //
    int             new_finest;
    Array<BoxArray> new_grid_places(max_level+1);
    
    if (lbase <= std::min(finest_level,max_level-1))
        grid_places(lbase, active_levels, dms, time, new_finest, new_grid_places);
    bool regrid_level_zero =
        (lbase == 0 && new_grid_places[0] != coarseRegion().boxArray()) && (!initial);
    const int start = regrid_level_zero ? 0 : lbase+1;

    //
    // If use_efficient_regrid flag is set, then test to see whether we in fact 
    // have changed the grids at any of the levels through the regridding process.
    // If not, then don't do anything more here.
    //
    if (use_efficient_regrid == 1 && !regrid_level_zero && (finest_level == new_finest) && !(initial && multi_region))
    {
        bool grids_unchanged = true; 
        for (int lev = start; lev <= finest_level && grids_unchanged; lev++)
        {
            if (new_grid_places[lev] != boxArray(lev)) grids_unchanged = false;
        }
        if (grids_unchanged) 
        {
            if (verbose > 0 && ParallelDescriptor::IOProcessor())
                std::cout << "Regridding at base region " << base_id << " but grids unchanged " << std::endl;
            return;
        }
    }
    
    std::list<int> old_structure = amr_regions.getStructure(base_id);
    //
    // We now remove the descendants from the hierarchy without deleting them
    // and put them in evictees.
    // Note that evictees is doubly managed so that it will delete
    // all the necessary Regions.
    //
    RegionList evictees(PListManage);
    amr_regions.extractChildrenOfNode(base_region.getID(), evictees);
    if (regrid_level_zero)
    {
        evictees.push_back(amr_regions.removeData(ROOT_ID));
    }

    //
    // Reclaim old-time grid space for all remain levels > lbase.
    //
    for (RegionList::iterator it = evictees.begin(); it != evictees.end(); it++)
            (*it)->removeOldData();

    finest_level = new_finest;

    if (lbase == 0)
    {
        MultiFab::FlushSICache();
        Geometry::FlushPIRMCache();
        FabArrayBase::CPC::FlushCache();
        DistributionMapping::FlushCache();
    }
    
    //
    // Determine the new region hierarchy using FOF clustering
    //
    Tree<BoxArray> grid_tree;
    if (start == 0)
        grid_tree.setRoot(new_grid_places[0]);
    else 
        grid_tree.setRoot(base_region.boxArray());
    for (int lev = lbase + 1; lev <= new_finest; lev++) 
    {
        // Cluster the new grids into regions
        std::list<BoxArray> clusters;
        if (region_creation == "FOF")
        {
            FOFCluster(0, lev, new_grid_places[lev], clusters);
        }
        else if (region_creation == "Single")
        {
            clusters.push_back(new_grid_places[lev]);
        }
        else if (region_creation == "Application")
        {
            coarseRegion().cluster(base_id, lev, new_grid_places[lev], clusters);
        }
        else
            BoxLib::Abort("Unrecognized Region creation mode\n");
            
        // Build the new clusters into a tree
        for (std::list<BoxArray>::iterator clust_it = clusters.begin(); clust_it != clusters.end(); ++clust_it)
        {
            BoxArray cba = *clust_it;
            cba.coarsen(refRatio(lev-1));
            TreeIterator<BoxArray> it = grid_tree.getIterator(lev-lbase-1);
            bool added = false;
            for ( ; !it.isFinished(); ++it)
            {
                if ((*it).contains(cba))
                {
                    grid_tree.addChildToNode(it.getID(), *clust_it);
                    added = true;
                    break;
                }
            }
            if (!added)
                BoxLib::Abort("Failed to find parent for new region");
        }
    }
    
    //
    // Restructure to reflect the new hierarchy
    // 
    restructure(base_id, grid_tree.getStructure(), false);
    
    std::list<int> new_structure = n_cycle.getStructure(base_id);
    
    
    //
    // Check to see if the structure has changed. Note, this check only
    // recognizes identical trees--if the leaves have been re-ordered
    // it may think it's a new structure.
    // This could be updated to check for fully isometric structures.
    //
    if(old_structure != new_structure)
    {
        if (verbose > 0 && ParallelDescriptor::IOProcessor())
        {
            Tree<std::string> st_tree;
            st_tree.buildFromStructure(n_cycle.getStructure());
            for(TreeIterator<std::string> sit = st_tree.getIterator(); !sit.isFinished(); ++sit)
                *sit = sit.getID().toString();
            std::cout << "Region Structure Changed--New Structure:\n" << st_tree.toString();
        }
    }
    
    //
    // Define the new grids
    // 
    TreeIterator<BoxArray> it = grid_tree.getIterator(Prefix);
    // Only touch lbase if we're regridding level 0.
    if (start != 0)
        ++it;
    for(; !it.isFinished(); ++it)
    {
        // Configure id's.
        ID id = it.getID();
        ID new_id = base_id;
        new_id.resize(id.size() + lbase);
        for (int i = 1; i < id.size(); i++)
            new_id[lbase+i] = id[i];
        int lev = new_id.level();
        
        //
        // Construct skeleton of new level.
        //
        AmrRegion* a = (*levelbld)(*this,new_id,geom[lev],*it,cumtime);
        if (initial)
        {
            //
            // We're being called on startup from bldFineLevels().
            // NOTE: The initData function may use a fillPatch, and so needs to
            //       be officially inserted into the hierarchy prior to the call.
            //
            ID parent_id = new_id.parent();
            amr_regions.addChildToNode(parent_id,a);
            a->initData();
        }
        else if (active_levels.defined(lev))
        {
            //
            // Init with data from old structure then remove old structure.
            // NOTE: The init function may use a filPatch from the old level,
            //       which therefore needs remain in the hierarchy during the call.
            //
            a->init(active_levels[lev]);
        }
        else
        {
            //
            // This is a new level of refinement within base_region,
            //
            a->init();
        }
        //
        // Add the new region to the hierarchy
        //
        if (!initial)
        {
        ID parent_id = new_id.parent();
        amr_regions.addChildToNode(parent_id,a);
        }
        if (lev > 0)
            touched_regions.push_back(a);
    }
    
    //
    // Build any additional data structures at the base region or higher after grid generation.
    //
    for (RegionList::iterator it = touched_regions.begin(); it != touched_regions.end(); it++)
        (*it)->post_regrid(base_id,new_finest);
        
    /// Is there any problem with changing dt here?
    //
    // Optimal subcycling codes must update the subcycling pattern for the subtree
    //
    if (subcycling_mode == "Optimal")
        base_region.computeRestrictedDt(base_id, n_cycle, dt_region);

#ifdef USE_STATIONDATA
    station.findGrid(amr_regions,geom);
#endif
    //
    // Report creation of new grids.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        printGridInfo(runlog,start,finest_level);
    }
    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        if (lbase == 0)
            gridlog << "STEP = " << region_steps.getRoot() << ' ';

        gridlog << "TIME = "
                << time
                << " : REGRID  with lbase = "
                << lbase
                << '\n';

        printGridInfo(gridlog,start,finest_level);
    }
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        if (lbase == 0)
            std::cout << "STEP = " << region_steps.getRoot() << ' ';

        std::cout << "TIME = "
                  << time
                  << " : REGRID  with base region = "
                  << base_id
                  << std::endl;
                  
        if (verbose > 1)
        {
           printGridInfo(std::cout,start,finest_level);
        }
        else
        {
           printGridSummary(std::cout,start,finest_level);
        }
    }
}

void
Amr::printGridInfo (std::ostream& os,
                    int           min_lev,
                    int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = boxArray(lev);
        int                       numgrid = bs.size();
        long                      ncells  = cellCount(lev);
        double                    ntot    = geom[lev].Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';


        RegionIterator it = amr_regions.getIterator(lev);
        for ( ; !it.isFinished(); ++it)
        {
            const BoxArray& ba = boxArray(it.getID());
            const DistributionMapping& map = getRegion(it.getID()).get_new_data(0).DistributionMap();
            for (int k = 0; k < ba.size(); k++)
            {
                const Box& b = ba[k];

                os << ' ' << lev << ": " << b << "   ";
                    
                for (int i = 0; i < BL_SPACEDIM; i++)
                    os << b.length(i) << ' ';

                os << ":: " << map[k] << '\n';
            }
        }
    }

    os << std::endl; // Make sure we flush!
}

void
Amr::printGridSummary (std::ostream& os,
                       int           min_lev,
                       int           max_lev)
{
    for (int lev = min_lev; lev <= max_lev; lev++)
    {
        const BoxArray&           bs      = boxArray(lev);
        int                       numgrid = bs.size();
        long                      ncells  = cellCount(lev);
        double                    ntot    = geom[lev].Domain().d_numPts();
        Real                      frac    = 100.0*(Real(ncells) / ntot);

        os << "  Level "
           << lev
           << "   "
           << numgrid
           << " grids  "
           << ncells
           << " cells  "
           << frac
           << " % of domain"
           << '\n';
    }

    os << std::endl; // Make sure we flush!
}

void
Amr::ProjPeriodic (BoxList&        blout,
                   const Geometry& geom)
{
    //
    // Add periodic translates to blout.
    //
    Box domain = geom.Domain();

    BoxList blorig(blout);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    D_TERM( nist , =njst , =nkst ) = -1;
    D_TERM( niend , =njend , =nkend ) = +1;

    int ri,rj,rk;
    for (ri = nist; ri <= niend; ri++)
    {
        if (ri != 0 && !geom.isPeriodic(0))
            continue;
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,ri*domain.length(0));
        for (rj = njst; rj <= njend; rj++)
        {
            if (rj != 0 && !geom.isPeriodic(1))
                continue;
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,rj*domain.length(1));
            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk != 0 && !geom.isPeriodic(2))
                    continue;
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,rk*domain.length(2));

                BoxList tmp(blorig);
                tmp.intersect(domain);
                blout.catenate(tmp);
 
                if (rk != 0 && geom.isPeriodic(2))
                    blorig.shift(2,-rk*domain.length(2));
            }
            if (rj != 0 && geom.isPeriodic(1))
                blorig.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && geom.isPeriodic(0))
            blorig.shift(0,-ri*domain.length(0));
    }
}

void
Amr::grid_places (int               lbase,
                  PArray<AmrRegion>& active_levels,
                  Array<DistributionMapping>& dms,
                  Real              time,
                  int&              new_finest,
                  Array<BoxArray>&  new_grids)
{
    BL_PROFILE("Amr::grid_places()");

    int i, max_crse = std::min(finest_level,max_level-1);

    const Real strttime = ParallelDescriptor::second();

    if (lbase == 0)
    {
        //
        // Recalculate level 0 BoxArray in case max_grid_size[0] has changed.
        // This is done exactly as in defBaseLev().
        //
        BoxArray lev0(1);

        lev0.set(0,BoxLib::coarsen(geom[0].Domain(),2));
        //
        // Now split up into list of grids within max_grid_size[0] limit.
        //
        lev0.maxSize(max_grid_size[0]/2);
        //
        // Now refine these boxes back to level 0.
        //
        lev0.refine(2);

        new_grids[0] = lev0;
    }

    if ( time == 0. && !initial_grids_file.empty() )
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min(new_finest,initial_ba.size());

        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = initial_ba[lev-1].size();
            for (i = 0; i < ngrid; i++)
            {
                Box bx(initial_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    else if ( !regrid_grids_file.empty() )
    {
        new_finest = std::min(max_level,(finest_level+1));
        new_finest = std::min(new_finest,regrid_ba.size());
        for (int lev = 1; lev <= new_finest; lev++)
        {
            BoxList bl;
            int ngrid = regrid_ba[lev-1].size();
            for (i = 0; i < ngrid; i++)
            {
                Box bx(regrid_ba[lev-1][i]);
                if (lev > lbase)
                    bl.push_back(bx);
            }
            if (lev > lbase)
                new_grids[lev].define(bl);
        }
        return;
    }

    //
    // Construct problem domain at each level.
    //
    Array<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Array<IntVect> rr_lev(max_level);
    Array<Box>     pc_domain(max_level);  // Coarsened problem domain.
    for (i = 0; i <= max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            bf_lev[i][n] = std::max(1,blocking_factor[i]/ref_ratio[i][n]);
    }
    for (i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for (i = lbase; i <= max_crse; i++)
        pc_domain[i] = BoxLib::coarsen(geom[i].Domain(),bf_lev[i]);
    //
    // Construct proper nesting domains.
    //
    Array<BoxList> p_n(max_level);      // Proper nesting domain.
    Array<BoxList> p_n_comp(max_level); // Complement proper nesting domain.

    BoxList bl(active_levels[lbase].boxArray());
    bl.simplify();
    bl.coarsen(bf_lev[lbase]);
    p_n_comp[lbase].complementIn(pc_domain[lbase],bl);
    p_n_comp[lbase].simplify();
    p_n_comp[lbase].accrete(n_proper);
    Amr::ProjPeriodic(p_n_comp[lbase], Geometry(pc_domain[lbase]));
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].simplify();
    bl.clear();
    for (i = lbase+1; i <= max_crse; i++)
    {
        p_n_comp[i] = p_n_comp[i-1];

        // Need to simplify p_n_comp or the number of grids can too large for many levels.
        p_n_comp[i].simplify();

        p_n_comp[i].refine(rr_lev[i-1]);
        p_n_comp[i].accrete(n_proper);

        Amr::ProjPeriodic(p_n_comp[i], Geometry(pc_domain[i]));

        p_n[i].complementIn(pc_domain[i],p_n_comp[i]);
        p_n[i].simplify();
    }
    //
    // Now generate grids from finest level down.
    //
    new_finest = lbase;

    for (int levc = max_crse; levc >= lbase; levc--)
    {
        int levf = levc+1;
        //
        // Construct TagBoxArray with sufficient grow factor to contain
        // new levels projected down to this level.
        //
        int ngrow = 0;

        if (levf < new_finest)
        {
            BoxArray ba_proj(new_grids[levf+1]);

            ba_proj.coarsen(ref_ratio[levf]);
            ba_proj.grow(n_proper);
            ba_proj.coarsen(ref_ratio[levc]);

            BoxArray levcBA = active_levels[levc].boxArray();

            while (!levcBA.contains(ba_proj))
            {
                BoxArray tmp = levcBA;
                tmp.grow(1);
                levcBA = tmp;
                ngrow++;
            }
        }
        TagBoxArray tags(active_levels[levc].boxArray(),dms[levc], n_error_buf[levc]+ngrow);

        active_levels[levc].errorEst(tags,
                                 TagBox::CLEAR,TagBox::SET,time,
                                 n_error_buf[levc],ngrow);
        //
        // If new grids have been constructed above this level, project
        // those grids down and tag cells on intersections to ensure
        // proper nesting.
        //
        // NOTE: this loop replaces the previous code:
        //      if (levf < new_finest) 
        //          tags.setVal(ba_proj,TagBox::SET);
        // The problem with this code is that it effectively 
        // "buffered the buffer cells",  i.e., the grids at level
        // levf+1 which were created by buffering with n_error_buf[levf]
        // are then coarsened down twice to define tagging at
        // level levc, which will then also be buffered.  This can
        // create grids which are larger than necessary.
        //
        if (levf < new_finest)
        {
            int nerr = n_error_buf[levf];

            BoxList bl_tagged(new_grids[levf+1]);
            bl_tagged.simplify();
            bl_tagged.coarsen(ref_ratio[levf]);
            //
            // This grows the boxes by nerr if they touch the edge of the
            // domain in preparation for them being shrunk by nerr later.
            // We want the net effect to be that grids are NOT shrunk away
            // from the edges of the domain.
            //
            for (BoxList::iterator blt = bl_tagged.begin(), End = bl_tagged.end();
                 blt != End;
                 ++blt)
            {
                for (int idir = 0; idir < BL_SPACEDIM; idir++)
                {
                    if (blt->smallEnd(idir) == geom[levf].Domain().smallEnd(idir))
                        blt->growLo(idir,nerr);
                    if (blt->bigEnd(idir) == geom[levf].Domain().bigEnd(idir))
                        blt->growHi(idir,nerr);
                }
            }
            Box mboxF = BoxLib::grow(bl_tagged.minimalBox(),1);
            BoxList blFcomp;
            blFcomp.complementIn(mboxF,bl_tagged);
            blFcomp.simplify();
            bl_tagged.clear();

            const IntVect iv = IntVect(D_DECL(nerr/ref_ratio[levf][0],
                                              nerr/ref_ratio[levf][1],
                                              nerr/ref_ratio[levf][2]));
            blFcomp.accrete(iv);
            BoxList blF;
            blF.complementIn(mboxF,blFcomp);
            BoxArray baF(blF);
            blF.clear();
            baF.grow(n_proper);
            //
            // We need to do this in case the error buffering at
            // levc will not be enough to cover the error buffering
            // at levf which was just subtracted off.
            //
            for (int idir = 0; idir < BL_SPACEDIM; idir++) 
            {
                if (nerr > n_error_buf[levc]*ref_ratio[levc][idir]) 
                    baF.grow(idir,nerr-n_error_buf[levc]*ref_ratio[levc][idir]);
            }

            baF.coarsen(ref_ratio[levc]);

            tags.setVal(baF,TagBox::SET);
        }
        //
        // Buffer error cells.
        //
        tags.buffer(n_error_buf[levc]+ngrow);
        //
        // Coarsen the taglist by blocking_factor/ref_ratio.
        //
        int bl_max = 0;
        for (int n=0; n<BL_SPACEDIM; n++)
            bl_max = std::max(bl_max,bf_lev[levc][n]);
        if (bl_max > 1) 
            tags.coarsen(bf_lev[levc]);
        //
        // Remove or add tagged points which violate/satisfy additional 
        // user-specified criteria.
        //
        active_levels[levc].manual_tags_placement(tags, bf_lev);
        //
        // Map tagged points through periodic boundaries, if any.
        //
        tags.mapPeriodic(Geometry(pc_domain[levc]));
        //
        // Remove cells outside proper nesting domain for this level.
        //
        tags.setVal(p_n_comp[levc],TagBox::CLEAR);
        //
        // Create initial cluster containing all tagged points.
        //
        long     len = 0;
        IntVect* pts = tags.collate(len);

        tags.clear();

        if (len > 0)
        {
            //
            // Created new level, now generate efficient grids.
            //
            new_finest = std::max(new_finest,levf);
            //
            // Construct initial cluster.
            //
            ClusterList clist(pts,len);
            clist.chop(grid_eff);
            BoxDomain bd;
            bd.add(p_n[levc]);
            clist.intersect(bd);
            bd.clear();
            //
            // Efficient properly nested Clusters have been constructed
            // now generate list of grids at level levf.
            //
            BoxList new_bx;
            clist.boxList(new_bx);
            new_bx.refine(bf_lev[levc]);
            new_bx.simplify();
            BL_ASSERT(new_bx.isDisjoint());
            IntVect largest_grid_size;
            for (int n = 0; n < BL_SPACEDIM; n++)
                largest_grid_size[n] = max_grid_size[levf] / ref_ratio[levc][n];
            //
            // Ensure new grid boxes are at most max_grid_size in index dirs.
            //
            new_bx.maxSize(largest_grid_size);

#ifdef BL_FIX_GATHERV_ERROR
	      int wcount = 0, iLGS = largest_grid_size[0];

              while (new_bx.size() < 64 && wcount++ < 4)
              {
                  iLGS /= 2;
                  if (ParallelDescriptor::IOProcessor())
                  {
                      std::cout << "BL_FIX_GATHERV_ERROR:  using iLGS = " << iLGS
                                << "   largest_grid_size was:  " << largest_grid_size[0]
                                << '\n';
                      std::cout << "BL_FIX_GATHERV_ERROR:  new_bx.size() was:   "
                                << new_bx.size() << '\n';
                  }

                  new_bx.maxSize(iLGS);

                  if (ParallelDescriptor::IOProcessor())
                  {
                      std::cout << "BL_FIX_GATHERV_ERROR:  new_bx.size() now:   "
                                << new_bx.size() << '\n';
                  }
	      }
#endif
            //
            // Refine up to levf.
            //
            new_bx.refine(ref_ratio[levc]);
            BL_ASSERT(new_bx.isDisjoint());
            new_grids[levf].define(new_bx);
        }
        //
        // Don't forget to get rid of space used for collate()ing.
        //
        delete [] pts;
    }

    // If Nprocs > Ngrids and refine_grid_layout == 1 then break up the grids
    //    into smaller chunks for better load balancing
    impose_refine_grid_layout(lbase,new_finest,new_grids);

    if (verbose > 0)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifndef BL_PROFILING
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "grid_places() time: " << stoptime << '\n';
#endif
    }
}

void
Amr::bldFineLevels (Real strt_time)
{
    finest_level = 0;

    Array<BoxArray> grids(max_level+1);
    ////
    //// Get initial grid placement.
    ////
    
    //
    // Setup data that will help us create the "Tree"
    // 
    ID parent_id = ROOT_ID;
    ID new_id = parent_id;
    std::list<int> structure;
    structure.push_back(1);
    structure.push_back(0);
    do
    {
        int new_finest;
        PArray<AmrRegion> active_levels;
        ID id;
        Array<DistributionMapping> dms(finest_level + 1);
        active_levels.resize(finest_level+1);
        for (int i = 0; i <= finest_level; i++)
        {
            id.resize(i+1);
            id[i] = 0;
            active_levels.set(i,&getRegion(id));
            dms[i] = active_levels[i].get_distribution_map();
        }
        grid_places(finest_level,active_levels,dms, strt_time,new_finest,grids);
        if (new_finest <= finest_level) break;
        finest_level = new_finest;
        
        new_id.resize(finest_level + 1);
        new_id[finest_level] = 0;
        //
        // Create a new level and link with others.
        // Construct skeleton of new level.
        //
        restructure(parent_id, structure);
        AmrRegion* level = (*levelbld)(*this,
                                      new_id,
                                      geom[new_finest],
                                      grids[new_finest],
                                      strt_time);
        amr_regions.setData(new_id, level);
        level->initData();
        
        parent_id = new_id;
    }
    while (finest_level < max_level);
    //
    // Iterate grids to ensure fine grids encompass all interesting gunk.
    //     but only iterate if we did not provide a grids file.
    //
    if ( regrid_grids_file.empty() || (strt_time == 0.0 && !initial_grids_file.empty()) )  
      {
	bool grids_the_same;

	const int MaxCnt = 4;

	int count = 0;

	do
	  {
	    for (int i = 0; i <= finest_level; i++)
	      grids[i] = boxArray(i);

	    regrid(ROOT_ID,strt_time,true);
        
	    grids_the_same = true;

        //Note: This only tracks whether the overall grid 
        // structure is the same. At the moment this <=> the region
        // grid structure being the same. However, if region creation
        // becomes more complicated this will need to update.
	    for (int i = 0; i <= finest_level && grids_the_same; i++)
            if (!(grids[i] == boxArray(i)))
                grids_the_same = false;

	    count++;
	  }
	while (!grids_the_same && count < MaxCnt);
      }
}

void
Amr::initSubcycle ()
{
    ParmParse pp("amr");
    sub_cycle = true;
    if (pp.contains("nosub"))
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "Warning: The nosub flag has been deprecated.\n ";
            std::cout << "... please use subcycling_mode to control subcycling.\n";
        }
        int nosub;
        pp.query("nosub",nosub);
        if (nosub > 0)
            sub_cycle = false;
        else
            BoxLib::Error("nosub <= 0 not allowed.\n");
        subcycling_mode = "None";
    }
    else 
    {
        subcycling_mode = "Auto";
        pp.query("subcycling_mode",subcycling_mode);
    }
    
    if (subcycling_mode == "None")
    {
        sub_cycle = false;
        TreeIterator<int> it = n_cycle.getIterator();
        for (; !it.isFinished(); ++it)
        {
            (*it) = 1;
        }
    }
    else if (subcycling_mode == "Manual")
    {
        int cnt = pp.countval("subcycling_iterations");

        if (cnt == 1)
        {
            //
            // Set all values to the single available value.
            //
            int cycles = 0;

            pp.get("subcycling_iterations",cycles);
            manual_n_cycle.resize(max_level);
            manual_n_cycle[0] = 1; // coarse level is always 1 cycle
            for (int i = 1; i <=max_level; i++)
            {
                manual_n_cycle[i] = cycles;
            }
        }
        else if (cnt > 1)
        {
            //
            // Otherwise we expect a vector of max_level values.
            //
            
            pp.getarr("subcycling_iterations",manual_n_cycle,0,max_level+1);
            if (manual_n_cycle[0] != 1)
            {
                BoxLib::Error("First entry of subcycling_iterations must be 1");
            }
        }
        else
        {
            BoxLib::Error("Must provide a valid subcycling_iterations if mode is Manual");
        }
        // Check validity
        for (int i = 1; i <=max_level; i++)
        {
            if (manual_n_cycle[i] > MaxRefRatio(i-1))
                BoxLib::Error("subcycling iterations must always be <= ref_ratio");
            if (manual_n_cycle[i] <= 0)
                BoxLib::Error("subcycling iterations must always be > 0");
        }
        // Set n_cycle from manual n_cycle
        TreeIterator<int> it = n_cycle.getIterator(Prefix);
        for (++it; !it.isFinished(); ++it)
        {
            *it = manual_n_cycle[it.getLevel()];
        }
    }
    else if (subcycling_mode == "Auto")
    {
        n_cycle.setRoot(1);
        TreeIterator<int> it = n_cycle.getIterator(Prefix);
        ++it;
        for (; !it.isFinished(); ++it)
        {
            (*it) = MaxRefRatio(it.getLevel()-1);
        } 
    }
    else if (subcycling_mode == "Optimal")
    {
        // if subcycling mode is Optimal, n_cycle is set dynamically.
        // We'll initialize it to be Auto subcycling.
        n_cycle.setRoot(1);
        TreeIterator<int> it = n_cycle.getIterator(Prefix);
        for (++it; !it.isFinished(); ++it)
        {
            (*it) = MaxRefRatio(it.getLevel()-1);
        }  
    }
    else
    {
        std::string err_message = "Unrecognzied subcycling mode: " + subcycling_mode + "\n";
        BoxLib::Error(err_message.c_str());
    }
}

void
Amr::initPltAndChk ()
{
    ParmParse pp("amr");

    pp.query("checkpoint_files_output", checkpoint_files_output);
    pp.query("plot_files_output", plot_files_output);

    pp.query("plot_nfiles", plot_nfiles);
    pp.query("checkpoint_nfiles", checkpoint_nfiles);
    //
    // -1 ==> use ParallelDescriptor::NProcs().
    //
    if (plot_nfiles       == -1) plot_nfiles       = ParallelDescriptor::NProcs();
    if (checkpoint_nfiles == -1) checkpoint_nfiles = ParallelDescriptor::NProcs();
    
    check_file_root = "chk";
    pp.query("check_file",check_file_root);

    check_int = -1;
    int got_check_int = pp.query("check_int",check_int);

    check_per = -1.0;
    int got_check_per = pp.query("check_per",check_per);

    if (got_check_int == 1 && got_check_per == 1)
    {
        BoxLib::Error("Must only specify amr.check_int OR amr.check_per");
    }

    plot_file_root = "plt";
    pp.query("plot_file",plot_file_root);

    plot_int = -1;
    int got_plot_int = pp.query("plot_int",plot_int);

    plot_per = -1.0;
    int got_plot_per = pp.query("plot_per",plot_per);

    if (got_plot_int == 1 && got_plot_per == 1)
    {
        BoxLib::Error("Must only specify amr.plot_int OR amr.plot_per");
    }

    stream_max_tries = 4;
    pp.query("stream_max_tries",stream_max_tries);
    stream_max_tries = std::max(stream_max_tries, 1);

    abort_on_stream_retry_failure = false;
    pp.query("abort_on_stream_retry_failure",abort_on_stream_retry_failure);

}

void
Amr::initRegions()
{
    ParmParse pp("amr");
    multi_region = false;
    pp.query("multi_region", multi_region);
    region_creation = "FOF";
    pp.query("region_creation", region_creation);
    if (region_creation != "FOF" && region_creation != "Single" && region_creation != "Application" )
        BoxLib::Error("Unrecognized region creation mode\n");
    if (!multi_region)
        region_creation = "Single";
}

bool
Amr::okToRegrid (ID region_id, int iteration)
{
    bool ok = true;
    int level = region_id.level();
    if (level == max_level)
        return false;
    ok = ok && amr_regions.getData(region_id).okToRegrid(iteration);
    return region_count.getData(region_id) >= regrid_int[level] && ok;
}

Real
Amr::computeOptimalSubcycling (Tree<int>& best, Tree<Real> dt_max, Tree<Real>& est_work, Tree<int>& cycle_max)
{
    BL_ASSERT(best.getStructure() == cycle_max.getStructure());
    // internally these represent the total number of steps at a level, 
    // not the number of cycles
    Tree<int> cycles;
    cycles.buildFromStructure(best.getStructure());
    cycles.setRoot(1);
    Real best_ratio = 1e200;
    int nprocs = ParallelDescriptor::NProcs();
    Real best_dt = 0;
    Real ratio;
    Real dt;
    Real work;
    Real cur_dt;
    int cur_cyc;
    int cur_max;
    int par_cyc;
    int limit = 1;
    ID id;
    ID parent_id;
    
    const Real strttime = ParallelDescriptor::second();
    
    // Limit dt_max by the fully iterated child dts.
    TreeIterator<Real> dt_it = dt_max.getIterator();
    for ( ; !dt_it.isFinished(); ++dt_it)
    {
        id = dt_it.getID();
        if (id.level() == 0)
            break;
        parent_id = dt_it.getParentID();
        dt_max.setData(id, std::min(dt_max.getData(id), dt_max.getData(parent_id)));
        Real cmax = dt_max.getData(id)*cycle_max.getData(id);
        dt_max.setData(parent_id, std::min(dt_max.getData(parent_id), cmax));
    }
    
    // This provides a memory efficient way to test all candidates
    TreeIterator<int> pti = cycle_max.getIterator(Prefix);
    long temp_cand;
    TreeIterator<int> ptic = cycle_max.getIterator();
    for (; !ptic.isFinished(); ++ptic)
        limit *= *ptic;
    // We do the next part in parallel.
    int m_start = (limit * ParallelDescriptor::MyProc() ) / nprocs;
    int m_end = (limit * (ParallelDescriptor::MyProc() + 1)) / nprocs;
    for (int candidate = m_start; candidate < m_end; candidate++)
    {
        temp_cand = candidate;
        dt = dt_max.getRoot();
        work = est_work.getRoot();
        pti = cycle_max.getIterator(Prefix);
        for (++pti; !pti.isFinished(); ++pti)
        {
            //get the id for this node
            id = pti.getID();
            parent_id = pti.getParentID();
            
            // All this gettin is probably inefficient. It can be streamlined if needed.
            par_cyc = cycles.getData(parent_id);
            cur_dt = dt_max.getData(id);
            cur_max = cycle_max.getData(id);
            
            // grab the relevant "digit" and shift over.
            cur_cyc = (1 + temp_cand%cur_max) * par_cyc;
            temp_cand /= cur_max;
            // check if we're subcycling too much or too little
            // This may require more nuance with nosub specific checks.
            if ( (cur_dt * (cur_cyc-1) >= dt) || ((cur_cyc < cur_max) && ((cur_cyc+1) * cur_dt <= dt)) )
            {
                work = 1e100;
                break;
            }
            cycles.setData(id,cur_cyc);
            dt = std::min(dt, cur_cyc * cur_dt);
            work += cur_cyc * est_work.getData(id);
        }
        ratio = work/dt;
        if (ratio < best_ratio) 
        {
            for (pti = cycle_max.getIterator(); !pti.isFinished(); ++pti)
            {
                id = pti.getID();
                best.setData(id,cycles.getData(id));
            }
            best_ratio = ratio;
            best_dt = dt;
        }
    }
    // Find the best ratio among all processors.
    Array<Real> ratios(nprocs);
    int best_proc = 0;
    ParallelDescriptor::Gather(&best_ratio, 1, ratios.dataPtr(), 1,  ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor())
    {
        Real cur_best = 1e100;
        for (int i = 0; i < nprocs; i++)
        {
            if (ratios[i] < cur_best)
            {
                cur_best = ratios[i];
                best_proc = i;
            }
        }
    }
    ParallelDescriptor::Bcast(&best_proc, 1, ParallelDescriptor::IOProcessorNumber());
    // Now have that processor send out its pattern
    int N = best.numElements();
    Array<int> best_arr(N);
    if (best_proc == ParallelDescriptor::MyProc())
    {
        best_arr = best.writeToArray();
        ParallelDescriptor::Bcast(best_arr.dataPtr(), N, best_proc);
    }
    else
    {
        ParallelDescriptor::Bcast(best_arr.dataPtr(), N, best_proc);
        best.readFromArray(best_arr);
    }
    ParallelDescriptor::Bcast(&best_dt, 1, best_proc);
    //
    // Now we convert best back to n_cycles format
    //
    pti = cycle_max.getIterator();
    for (; !pti.isFinished(); ++pti)
    {
        //get the id for this node
        id = pti.getID();
        if (id.level() == 0) break;
        //get the parent id by cutting the last digit.
        parent_id = id.parent();
        best.setData(id,best.getData(id)/best.getData(parent_id));
    }
    
    if (verbose > 0)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "computeOptimalSubcycling() time: " << stoptime << '\n';
        }
    }
    return best_dt;
}

void
Amr::setRestrictedSubcycling(ID  base_region,
                             Tree<int>&   ncycle, 
                             Tree<Real>& dtregion, 
                             Tree<Real>& dt_max, 
                             Tree<int>&  cycle_max)
{
    
    TreeIterator<int> it = ncycle.getIterator(base_region, Prefix);
    ID id;
    ID parent_id;
    // Limit dt_max by the fully iterated child dts.
    TreeIterator<Real> dt_it = dt_max.getIterator(base_region);

    for ( ; !dt_it.isFinished(); ++dt_it)
    {
        id = dt_it.getID();
        if (id.size() == base_region.size())
            break;
        parent_id = dt_it.getParentID();
        Real cmax = dt_max.getData(id)*cycle_max.getData(id);
        dt_max.setData(parent_id, std::min(dt_max.getData(parent_id), cmax));
    }
   
    // BL_ASSERT(dtregion.getData(base_region) <= dt_max.getData(base_region));

    // We change the dt of the parent instead of failing here ...
    if (dtregion.getData(base_region) > dt_max.getData(base_region))
    {
       std::cout << "... DT REGION OF BASE REGION IS " << dtregion.getData(base_region) << std::endl;
       std::cout << "... DT    MAX OF BASE REGION IS " << dt_max.getData(base_region) << std::endl;
       std::cout << "... DT REGION OF BASE REGION > DT MAX OF BASE REGION" << std::endl;
       std::cout << "... so RESETTING DTREGION OF BASE REGION TO " << dt_max.getData(base_region) << std::endl;

       dtregion.setData(base_region,dt_max.getData(base_region));
       std::cout << "... DT REGION OF BASE REGION IS NOW         " << dtregion.getData(base_region) << std::endl;
    }


    for (++it ; !it.isFinished(); ++it)
    {
        id = it.getID();
        parent_id = it.getParentID();
        TreeIterator<Real> dt_it = dt_max.getIterator(id);

        // Necessary dt / allowed dt
        int new_cycle = int(ceil(dtregion.getData(parent_id) / dt_max.getData(id)));
        if (new_cycle > cycle_max.getData(id))
        {
            BoxLib::Abort("Invalid dt after regrid: unable to subcycle enough\n");
        }
        ncycle.setData(id, new_cycle);
        dtregion.setData(id, dtregion.getData(parent_id) / new_cycle);
    }
}

void
Amr::FindMaxDt(Real& dt_0, Tree<int> n_cycle, Tree<Real> dt_level)
{
    Tree<Real> dt_max(dt_level);
    ID id;
    ID parent_id;
    TreeIterator<Real> dt_it = dt_max.getIterator();
    for (;!dt_it.isFinished(); ++dt_it)
    {
        id = dt_it.getID();
        if (id.size() > 1) //not at root
        {
            parent_id = dt_it.getParentID();
            dt_max.setData(parent_id, std::min(dt_max.getData(parent_id), n_cycle.getData(id)*dt_max.getData(id)));
        }
    }
    dt_0 = std::min(dt_0, dt_max.getRoot());
}

void
Amr::FOFCluster (int d, int level, BoxArray boxes, std::list<BoxArray>& cluster_list )
{
    std::list<BoxList>::iterator it;
    int N = boxes.size();
    std::list<BoxList> clusters;
    Array<IntVect> pisect;
    Geometry& gm = geom[level];
    for (int i = 0;i < N; i++)
    {
        BoxList* cluster_id = 0;
        const Box box = boxes[i];
        Box bgrown = boxes[i];
        // We grow the box by d+1; this way it will intersect boxes 
        // that are within distance d of the real box.
        bgrown.grow(d+1);
        for (it = clusters.begin(); it != clusters.end(); it++)
        {
            pisect.resize(0);
            bool friends = false;
            if (!BoxLib::intersect(*it, bgrown).isEmpty())
            {
                friends = true;
            }
            else if (gm.isAnyPeriodic())
            {
                for (BoxList::iterator bit = (*it).begin(); bit!=(*it).end(); ++bit)
                {
                    gm.periodicShift(bgrown, *bit, pisect);
                    if (!pisect.empty())
                    {
                        friends = true;
                        break;
                    }
                }
            }
            
            if (friends)
            {
                if (cluster_id == 0)
                {
                    // This is the first cluster this box belongs to
                    cluster_id = &(*it);
                    it->push_back(box);
                }
                else
                {
                    // This box belongs to multiple clusters; combine them.
                    cluster_id->join(*it);
                    clusters.erase(it);
                    it--;
                }
            }
                
                
        }
        if (cluster_id == 0)
        {
            // We didn't find a cluster. Make a new one.
            BoxList* new_cluster = new BoxList(box);
            clusters.push_back(*new_cluster);
        }
    }
;
    for (std::list<BoxList>::iterator it = clusters.begin(); it != clusters.end(); ++it)
    {
        BoxArray ba(*it);
        cluster_list.push_back(ba);
    }
    return; 
}

AmrRegion*
Amr::build_blank_region()
{
    return (*levelbld)();
}

void
Amr::aggregate_descendants(const ID id, PArray<AmrRegion>& aggregates)
{
    int base_level = id.size()-1;
    int num_levels = finest_level - base_level + 1;
    // Initialize the lists of descendants
    PArray<RegionList> descendant_list(PArrayManage);
    descendant_list.resize(num_levels);
    for (int i = 0; i < num_levels; i++)
    {
        RegionList* rl = new RegionList(PListNoManage);
        descendant_list.set(i,rl);
    }
    // Fill the lists
    int max_level = 0;
    RegionIterator it = amr_regions.getIterator(id);
    for ( ; !it.isFinished(); ++it)
    {
        int cur_level = (*it)->Level();
        max_level = std::max(max_level, cur_level);
        descendant_list[cur_level- base_level].push_back(*it);
    }
    // Create Aggregates
    BL_ASSERT(max_level <= finest_level);
    aggregates.resize(max_level - base_level + 1);
    for (int i = 0; i < max_level - base_level + 1; i++)
    {
        // The user is responsible for deleting this pointer
        // Ideally aggregates should be a managed array.
        AmrRegion* a = build_blank_region();
        a->define(descendant_list[i], this);
        aggregates.set(i,a);
    }
}

ID
Amr::whichRegion(int level, IntVect cell)
{
    RegionIterator it = amr_regions.getIterator(level);
    for ( ; !it.isFinished(); ++it)
    {
        if ((*it)->boxArray().contains(cell))
            return (*it)->getID();
    }
    BoxLib::Abort("Unable to find region containing specified IntVect");
    return ROOT_ID; // to remove compiler warnings.
}

MultiFab* 
Amr::createMultiFab(ID     base_region, 
                            int             lev, 
                            ExecutionTree*  exec_tree, 
                            int             ncomp,
                            int             ngrow)
{
    MultiFab* result = new MultiFab();
    PList<MultiFab> mf_list;
    ExecutionTreeIterator it = exec_tree->getIterator(base_region, lev);
    for ( ; !it.isFinished(); ++it)
    {
        mf_list.push_back(new MultiFab(amr_regions.getData(it.getID()).boxArray(), ncomp, ngrow));
    }
    result->define(mf_list);
    return result;
}

DistributionMapping
Amr::createDM(ID     base_region, 
              int             lev,
              ExecutionTree* exec_tree)
{
    Array<DistributionMapping> da;
    int nprocs = ParallelDescriptor::NProcs();
    int box_count = 0;
    if (exec_tree)
    {
        ExecutionTreeIterator it = exec_tree->getIterator(base_region, lev);
        for ( ; !it.isFinished(); ++it)
        {
            // inefficient
            da.resize(da.size() + 1);
            BoxArray ba = boxArray(it.getID());
            box_count += ba.size();
            da[da.size()-1] = DistributionMapping(ba, nprocs);
        }
    }
    else
    {
        RegionIterator it = getRegionIterator(base_region, lev);
        for ( ; !it.isFinished(); ++it)
        {
            // inefficient
            da.resize(da.size() + 1);
            BoxArray ba = boxArray(it.getID());
            box_count += ba.size();
            da[da.size()-1] = DistributionMapping(ba, nprocs);
        }
    }
    DistributionMapping dm;
    dm.define(box_count, da);
    return dm;
}

void 
Amr::impose_refine_grid_layout (int lbase, int new_finest, Array<BoxArray>& new_grids)
{
    const int NProcs = ParallelDescriptor::NProcs();
    if (NProcs > 1 && refine_grid_layout)
    {
        //
        // Chop up grids if fewer grids at level than CPUs.
        // The idea here is to make more grids on a given level
        // to spread the work around.
        //
        for (int cnt = 1; cnt <= 4; cnt *= 2)
        {
            for (int i = lbase; i <= new_finest; i++)
            {
                const int ChunkSize = max_grid_size[i]/cnt;

                IntVect chunk(D_DECL(ChunkSize,ChunkSize,ChunkSize));
                //
                // We go from Z -> Y -> X to promote cache-efficiency.
                //
                for (int j = BL_SPACEDIM-1; j >= 0 ; j--)
                {
                    chunk[j] /= 2;

                    if ( (new_grids[i].size() < NProcs) && (chunk[j]%blocking_factor[i] == 0) )
                    {
                        new_grids[i].maxSize(chunk);
                    }
                }
            }
        }
    }
}
