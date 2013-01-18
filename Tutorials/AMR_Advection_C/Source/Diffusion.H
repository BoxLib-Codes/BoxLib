#ifndef _Diffusion_H_
#define _Diffusion_H_

#include <AmrRegion.H>
#include <FluxRegister.H>
#include <Tree.H>

class Diffusion {

public:

  Diffusion(Amr* Parent, int _finest_level, BCRec* _phys_bc);
  virtual ~Diffusion ();

  void read_params ();

  void install_region (ID          region_id,
                       AmrRegion*  region_data_to_install,
                       MultiFab&   volume,
                       MultiFab*   area);

  // Level 0 version 
  void applyop(ID region_id,MultiFab& species,
               MultiFab& DiffTerm, PArray<MultiFab>& diff_coef);

  // Level > 0 version 
  void applyop(ID region_id,MultiFab& species, MultiFab& Crse,
               MultiFab& DiffTerm, PArray<MultiFab>& diff_coef);

  void make_mg_bc();

  void restructure(ID base_region, std::list<int> structure);

  void GetCrsePhi(ID region_id,
                  MultiFab& phi_crse,
                  Real time);

  void zeroPhiFluxReg (ID region_id);


protected:
  //
  // Pointers to amr,amrlevel.
  //
  Amr*             master;
  PTree<AmrRegion> region_data;

  PTree<FluxRegister> phi_flux_reg;
  //
  // BoxArray at each level
  //
  Tree<BoxArray> grids;
  //
  // Volume and area fractions.
  //
  PTree<MultiFab> volume;
  Tree<MultiFab*> area;

  int mg_bc[2*BL_SPACEDIM];

  BCRec*       phys_bc;

  static int   verbose;
  static int   stencil_type;

#if (BL_SPACEDIM < 3)
  void applyMetricTerms(ID region,MultiFab& Rhs, PArray<MultiFab>& coeffs);
  void unweight_cc(ID region,MultiFab& cc);
#endif

  ExecutionTree* init_solve_tree();
};
#endif

