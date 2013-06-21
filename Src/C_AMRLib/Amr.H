
#ifndef _Amr_H_
#define _Amr_H_ 

#include <fstream>

#include <Box.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>
#include <PArray.H>
#include <PList.H>
#include <PTree.H>
#include <Tree.H>
#include <ParmParse.H>

#ifdef USE_STATIONDATA
#include <StationData.H>
#endif

class AmrRegion;
class LevelBld;
class BoxDomain;
typedef PTreeIterator<AmrRegion> RegionIterator;

//
//
// Manages hierarchy of levels for time-dependent AMR computations.
//
// The Amr class is designed to manage parts of the computation  which do
// not belong on a single level, like establishing and updating the hierarchy
// of levels, global timestepping, and managing the different AmrRegions
//

class Amr
{
public:
    // Typedef the list of regions for amr_level stuff.
    typedef PList<AmrRegion> RegionList;
    
    //
    // The constructor.
    //
    Amr ();
    //
    // The destructor.
    //
    virtual ~Amr ();
    //
    // set up t_nodes (for SDC) after construction. Must be called before init if doing SDC.
    //
    void set_t_nodes (const Array<Real>& tnodes);
    //
    // Init data after construction. Must be called before timestepping.
    //
    virtual void init (Real strt_time, Real stop_time);
    //
    // return blocking_factor
    //
    int blockingFactor (int lev) const { return blocking_factor[lev]; }
    //
    // Largest allowable grid.
    //
    int maxGridSize (int lev) const { return max_grid_size[lev]; }
    //
    // Subcycle in time?
    //
    int subCycle () const { return sub_cycle; }
    //
    // How are we subcycling?
    //
    const std::string& subcyclingMode() const { return subcycling_mode; }
    //
    // Maximum allowable level.
    //
    int maxLevel () const { return max_level; }
    //
    // Finest currently existing level.
    //
    int finestLevel () const { return finest_level; }
    //
    // Physical time.
    //
    Real cumTime () const { return cumtime; }
    void setCumTime (Real t) {cumtime = t;}
    //
    // Physical time this simulation started
    //
    Real startTime () const { return start_time; }
    void setStartTime (Real t) {start_time = t;}
    //
    // Time step at specified level.
    // Only allowed when multi_region = 0
    //
    Real dtLevel (int level) const;
    //
    // Tree of timesteps at all regions.
    //
    const Tree<Real>& dtRegion () { return dt_region; }
    //
    // Timestep at a specific region.
    //
    Real dtRegion(ID id) const { return dt_region.getData(id); }
    //
    // Array of inner node times on specified level 
    //
    const Array<Real>& tNodes (int level) const {
	return (multi_level_sdc) ? t_nodes[level] : t_nodes[0]; }
    //
    // Multi-level SDC?
    //
    bool isMultiLevSDC () const { return multi_level_sdc; }
    //
    // Geometry object.
    //
    const Geometry& Geom (int level) const { return geom[level]; }
    //
    // Number of subcycled time steps.
    // Only allowed when multi_region = 0
    //
    int nCycle (int lev) const;
    //
    // Number of subcycled time steps.
    //
    int nCycle (ID region) const { return n_cycle.getData(region); }
    //
    // Refinement ratio.
    //
    IntVect refRatio (int _level) const { return ref_ratio[_level]; }
    //
    // Maximum refinement ratio in any direction.
    //
    int MaxRefRatio (int level) const;
    //
    // Refinement ratios between all levels.
    //
    const Array<IntVect>& refRatio () const { return ref_ratio; }
    //
    // Number of time steps at specified level.
    // Only allowed when multi_region = 0
    //
    int levelSteps (int lev) const;
    //
    // Number of time steps at specified region.
    //
    int regionSteps (ID region) const;
    //
    // Which step are we at for the specified level?
    // Only allowed when multi_region = 0
    //
    int levelCount (int lev) const;
    //
    // Whether to regrid right after restart
    //
    bool RegridOnRestart () const;
    //
    // Interval between regridding.
    //
    int regridInt (int lev) const { return regrid_int[lev]; }
    //
    // Number of time steps between checkpoint files.
    //
    int checkInt () const { return check_int; }
    //
    // Time between checkpoint files.
    //
    Real checkPer() const { return check_per; }
    //
    // Number of time steps between plot files.
    //
    int plotInt () const { return plot_int; }
    //
    // Time between plot files.
    //
    Real plotPer () const { return plot_per; }
    //
    // The names of state variables to output in the
    // plotfile.  They can be set using the amr.plot_vars variable
    // in a ParmParse inputs file. 
    //
    static const std::list<std::string>& statePlotVars () { return state_plot_vars; }
    //
    // Is the string the name of a variable in state_plot_vars?
    //
    static bool isStatePlotVar (const std::string& name);
    //
    // If the string is not the name of a variable in state_plot_vars,
    // add it to state_plot_vars.  
    //
    static void addStatePlotVar (const std::string& name);
    //
    // Remove the string from state_plot_vars.
    //
    static void deleteStatePlotVar (const std::string& name);
    //
    // Clear the list of state_plot_vars.
    //
    static void clearStatePlotVarList ();
    //
    //  Fill the list of state_plot_vars with all of the state quantities.
    //
    static void fillStatePlotVarList ();
    //
    //  Write out plotfiles (True/False)?
    //
    static bool Plot_Files_Output ();
    //
    // The names of derived variables to output in the
    // plotfile.  They can be set using the amr.derive_plot_vars 
    // variable in a ParmParse inputs file.
    //
    static const std::list<std::string>& derivePlotVars () { return derive_plot_vars; }
    //
    // Is the string the name of a variable in derive_plot_vars?
    //
    static bool isDerivePlotVar (const std::string& name);
    //
    // If the string is not the name of a variable in 
    // derive_plot_vars, add it to derive_plot_vars.  
    //
    static void addDerivePlotVar (const std::string& name);
    //
    // Remove the string from derive_plot_vars.
    //
    static void deleteDerivePlotVar (const std::string& name);
    //
    // Clear the list of derive_plot_vars.
    //
    static void clearDerivePlotVarList ();
    //
    //  Fill the list of derive_plot_vars with all derived quantities.
    //
    static void fillDerivePlotVarList ();

    static void Initialize ();
    static void Finalize ();
    //
    // Number of buffer cells in error estimator.
    //
    int nErrorBuf (int lev) const { return n_error_buf[lev]; }
    //
    // Minimum allowable grid efficiency.
    //
    Real gridEff () const { return grid_eff; }
    //
    // Number of cells to define proper nesting 
    //
    int nProper () const { return n_proper; }
    //
    // List of grids in various portions of the domain.
    //
    BoxArray boxArray (int lev) const;
    BoxArray boxArray (ID base_region, int lev) const;
    BoxArray boxArray (ID base_region, int lev, ExecutionTree* exec_tree) const;
    BoxArray boxArray (ID region_id) const;
    //
    // Creates an empty multifab built on the boxArrays of the specified domain
    // portion. The importance of this method is that the resulting multifab's
    // distribution map is the union of the maps for the component boxArrays,
    // rather than a separate map for the new, larger mf.
    //
    MultiFab* createMultiFab(ID     base_region, 
                            int             lev, 
                            ExecutionTree*  exec_tree, 
                            int             ncomp,
                            int             ngrow);
    
    //
    // Creates a Distribution Mapping to represent the aggregated
    // distribution maps at the specified level.
    //
    DistributionMapping createDM(ID     base_region, 
                                 int             lev,
                                 ExecutionTree* exec_tree = 0);
    //
    // gets the AmrRegion at a given level. 
    // Only allowed when multi_region = 0
    //
    AmrRegion&  getLevel (int lev);
    //
    // Total number of cells.
    //
    long cellCount ();
    //
    // Number of cells at given level.
    //
    long cellCount (int lev);
    //
    // Number of cells at given region.
    //
    long cellCount (ID region_id);
    //
    // Total number of grids.
    //
    int numGrids ();
    //
    // Number of grids at given level.
    //
    int numGrids (int lev);
    //
    // More work to be done?
    //
    int okToContinue ();
    //
    // Regrid only!
    //
    void RegridOnly (Real time);
    //
    // Should we regrid this region?
    //
    bool okToRegrid (ID region_id, int iteration = 0);
    //
    // Do a complete integration cycle.
    //
    virtual void coarseTimeStep (Real stop_time);
    //
    // Retrieve derived data FOR THE ENTIRE LEVEL. User is responsible for deleting pointer.
    //
    MultiFab* derive (const std::string& name,
                      Real           time,
                      int            lev,
                      int            ngrow);
    //
    // Name of the restart file.
    //
    const std::string& theRestartFile () const { return restart_file; }
    //
    // The ith datalog file.  Do with it what you want.
    //
    std::ostream& DataLog (int i);
    //
    // How many datalogs have been opened
    //
    int NumDataLogs ();
    //
    // Compute the optimal subcycling pattern.
    // This assumes that anything less than cycle_max[i] is a valid
    // number of subcycles at level[i]. For example 
    // if ref_ratio[i] = cycle_max[i] = 4, then 1,2,3,4 are all valid 
    // values for n_cycles[i]
    //
    // This method is essentially static apart from the verbosity check.
    //
    Real computeOptimalSubcycling(Tree<int>& best,
                                  Tree<Real> dt_max, 
                                  Tree<Real>& est_work, 
                                  Tree<int>& cycle_max);
             
    //
    // Compute and set a new subcycling pattern for the subtree of base_region.
    // Descendents subcycle just enough to respect their max dt's while
    // allowing the base_region to keep its dt the same.
    //
    // This method is typically called after a regrid on a subtree.
    //
    static void setRestrictedSubcycling(ID  base_region,
                             Tree<int>&   n_cycle, 
                             Tree<Real>& dt_region, 
                             Tree<Real>& dt_max, 
                             Tree<int>&  cycle_max);
    
    //
    // Uses a friends-of-friends clustering algorithm with min distance d
    // to identify clusters in boxes. Distance is computed as the minimum '
    // distance between cells in the 2 boxes.
    // 
    void FOFCluster(int d, int level, BoxArray boxes, std::list<BoxArray>& cluster_list);
    
    AmrRegion* build_blank_region();
    
    void aggregate_descendants(const ID id, PArray<AmrRegion>& aggregates);
    
    //
    // Returns a reference to the coarse region of the simulation.
    //
    AmrRegion& coarseRegion() { return *(amr_regions.getRoot()); }
    
    //
    // Return a reference to the specified region.
    //
    AmrRegion& getRegion(const ID id) { return amr_regions.getData(id); }
    
    //
    // Locate the region that contains the specified cell. Aborts on failure.
    //
    ID whichRegion(int level, IntVect cell);
    
    //
    // PTree of all the AmrRegions.
    //
    PTree<AmrRegion>& getRegions();
    const PTree<AmrRegion>& getRegions() const;
    
    RegionIterator getRegionIterator(IteratorMode mode);
    RegionIterator getRegionIterator(int lev = -1, IteratorMode mode = Postfix);
    RegionIterator getRegionIterator(const ID region_id, int lev = -1, IteratorMode mode = Postfix);
    
    void FindMaxDt(Real& dt_0, Tree<int> n_cycle, Tree<Real> dt_level);
                                         
    //
    // Write the plot file to be used for visualization.
    //
    virtual void writePlotFile ();
    int stepOfLastPlotFile () const {return last_plotfile;}
    //
    // Write current state into a chk* file.
    //
    virtual void checkPoint ();
    int stepOfLastCheckPoint () const {return last_checkpoint;}

protected:
    //
    // Used by grid_places().
    //
    static void ProjPeriodic (BoxList&        bd,
                              const Geometry& geom);

    //
    // Array of BoxArrays read in to initially define grid hierarchy
    //
    static Array<BoxArray> initial_ba;

    //
    // Array of BoxArrays read in to externally define grid hierarchy at each regrid
    //
    static Array<BoxArray> regrid_ba;

    //
    // Initialize grid hierarchy.
    //
    void initialInit (Real strt_time, Real stop_time);
    //
    // Read the probin file.
    //
    void readProbinFile (int& init);
    //
    // Check for valid input.
    //
    void checkInput ();
    //
    // Restart from a checkpoint file.
    //
    void restart (const std::string& filename);
    //
    // Define and initialize coarsest level.
    //
    void defBaseLevel (Real start_time);
    //
    // Define and initialize refined levels.
    //
    void bldFineLevels (Real start_time);
    //
    // Rebuild grid hierarchy finer than lbase.
    //
    void regrid (int  lbase,
                 Real time,
                 bool initial = false);
    //
    // Rebuild grid hierarchy for descendants of base_region.
    //
    void regrid (ID base_id,
                 Real time,
                 bool initial = false);
                 
    //
    // Restructures our data trees in preparation for recieving new data
    // Normally called in regrid.
    //
    void restructure(ID base_region, std::list<int> structure, bool do_regions = true);
    
    //
    // Define new grid locations (called from regrid) and put into new_grids.
    //
    void grid_places (int               lbase,
                      PArray<AmrRegion>& active_levels,
                      Array<DistributionMapping>& dms,
                      Real              time,
                      int&              new_finest,
                      Array<BoxArray>&  new_grids);

    //
    // If ngrids < nprocs, break up the grids to better load balance
    //
    void impose_refine_grid_layout (int              lbase,
                                    int              new_finest,
                                    Array<BoxArray>& new_grids);

    //
    // Do a single timestep on the base region.
    //
    void timeStep (AmrRegion& base_region,
                   Real dt,
                   int  iter,
                   int  ncycle,
                   Real stop_time);

    //
    // Whether to write a plotfile now
    //
    bool writePlotNow ();

    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

    void printGridSummary (std::ostream& os,
                           int           min_lev,
                           int           max_lev);

    void setRecordGridInfo (const std::string&);

    void setRecordRunInfo (const std::string&);

    void setRecordRunInfoTerse (const std::string&);

    void setRecordDataInfo (int i, const std::string&);
    
    //
    // The following methods contain the parmParse code for their section.
    //
    void initSubcycle();
    void initPltAndChk();
    void initRegions();

    //
    // The data ...
    //
    std::string      regrid_grids_file;   // Grids file that will bypass regridding.
    std::string      initial_grids_file;  // Grids file that will bypass regridding only at initialization.
    int              max_level;    // Maximum allowed level.
    int              finest_level; // Current finest level.
    int              multi_level_sdc; // multi-level SDC?
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]
    PTree<AmrRegion> amr_regions;
    Array<Geometry>  geom;         // Index extent of domain at all levels.
    Real             cumtime;      // Physical time variable.
    Real             start_time;   // Physical time this simulation started.

    Tree<Real>       dt_region;     // Timestep for this region.
    Tree<int>        region_steps;  // This is unreliable right now--I'm not sure what use it is
    Tree<int>        region_count;
    Tree<int>        n_cycle;
    Array<int>       manual_n_cycle;

    Array<Array<Real> > t_nodes;   // Time of (SDC) inner nodes in the range of (0,1)

    std::string      subcycling_mode; //Type of subcycling to use.
    Tree<Real>       dt_reg_min; 

    bool             isPeriodic[BL_SPACEDIM];  // Domain periodic?
    Array<int>       regrid_int;      // Interval between regridding.
    Array<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Array<int>       blocking_factor; // Blocking factor in grid generation (by level).
    Array<int>       max_grid_size;   // Maximum allowable grid size (by level).
    Real             grid_eff;        // Grid efficiency.
    int              n_proper;        // # cells required for proper nesting.
    int              last_checkpoint; // Step number of previous checkpoint.
    int              check_int;       // How often checkpoint (# time steps).
    Real             check_per;       // How often checkpoint (units of time).
    std::string      check_file_root; // Root name of checkpoint file.
    int              last_plotfile;   // Step number of previous plotfile.
    int              plot_int;        // How often plotfile (# of time steps)
    Real             plot_per;        // How often plotfile (in units of time)
    int              file_name_digits; // How many digits to use in the plotfile and checkpoint names
    std::string      plot_file_root;  // Root name of plotfile.
#ifdef USE_STATIONDATA
    StationData      station;
#endif
    int              verbose;
    int              record_grid_info;
    int              record_run_info;
    int              record_run_info_terse;
    std::ofstream    gridlog;
    std::ofstream    runlog;
    std::ofstream    runlog_terse;
    PArray<std::ofstream> datalog;
    int              sub_cycle;
    std::string      restart_file;
    std::string      probin_file;
    LevelBld*        levelbld;

    // Region Controls
    // Allow multiple regions per level?
    bool             multi_region;
    // How do we create regions?
    std::string      region_creation;

    bool             abort_on_stream_retry_failure;
    int              stream_max_tries;

    bool             bUserStopRequest;
    //
    // The static data ...
    //
    static std::list<std::string> state_plot_vars;  // State Vars to dump to plotfile 
    static std::list<std::string> derive_plot_vars; // Derived Vars to dump to plotfile 
    static bool                   first_plotfile;
    

private:
    //
    // Disallowed
    //
    Amr (const Amr& rhs);
    Amr& operator= (const Amr& rhs);
};

inline
PTree<AmrRegion>&
Amr::getRegions()
{
    return amr_regions;
}

inline
const PTree<AmrRegion>&
Amr::getRegions() const
{
    return amr_regions;
}

inline
RegionIterator
Amr::getRegionIterator(int lev, IteratorMode mode)
{
    return amr_regions.getIterator(mode, lev);
}

inline
RegionIterator
Amr::getRegionIterator(IteratorMode mode)
{
    return amr_regions.getIterator(mode);
}

inline
RegionIterator
Amr::getRegionIterator(const ID region_id, int lev, IteratorMode mode)
{
    return amr_regions.getIterator(region_id, mode, lev);
}

#endif /*_Amr_H_*/
