
#include <SDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>

#include <cassert>

#include <stdio.h>

using namespace std;


extern "C" {
  void mlsdc_amr_interpolate(void *F, void *G, void *ctxF, void *ctxG)
  {

  }

  void mlsdc_amr_restrict(void *F, void *G, void *ctxF, void *ctxG)
  {

  }
}


void SDCAmr::timeStep (int  level,
                         Real time,
                         int  iteration,
                         int  niter,
                         Real stop_time)
{
  // assert(level == 0);

  // XXX

  // for (int k=0; k<max_iters; k++) {
  //   sdc_mg_sweep(&mg, time, dt_level[0], 0);
  // }

  cout << "!!! USING SDC TIME STEP !!!" << endl;

  Amr::timeStep(level, time, iteration, niter, stop_time);
}


SDCAmr::SDCAmr (sdc_level_bld_f bld)
{
    Initialize();
    InitAmr();
    InitSDCAmr(bld);
}

void
SDCAmr::InitSDCAmr(sdc_level_bld_f bld)
{
  // get parameters
  //
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 22;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 3;

  sdc_level_bld = bld;
  
  //
  // build MultiFab encapsulation for SDCLib
  //
  encap_ctx.ncomp = 1;
  encap_ctx.ngrow = 1;
  // build_encap();

  //
  // build multigrid sdc sweeper, add coarsest level
  //
  sweepers.push_back(sdc_level_bld(0));

  sdc_mg_build(&mg, max_level);
  sdc_mg_add_level(&mg, sweepers[0], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  
  sdc_mg_setup(&mg);
  sdc_mg_allocate(&mg);
}


SDCAmr::~SDCAmr()
{
  sdc_mg_destroy(&mg);
}
