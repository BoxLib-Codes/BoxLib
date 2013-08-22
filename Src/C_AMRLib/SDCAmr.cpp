
#include <SDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>

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
  // if (sweepers[0] == NULL) {
  //   build_encap(0);
  //   sweepers[0] = sdc_level_bld(*this, 0);
  //   sdc_mg_add_level(&mg, sweepers[0], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  //   sdc_mg_setup(&mg);
  //   sdc_mg_allocate(&mg);
  // }

  // assert(level == 0);

  // XXX

  // for (int k=0; k<max_iters; k++) {
  //   sdc_mg_sweep(&mg, time, dt_level[0], 0);
  // }

  cout << "!!! USING SDC TIME STEP !!!" << endl;

  Amr::timeStep(level, time, iteration, niter, stop_time);
}

void SDCAmr::regrid (int  lbase,
		     Real time,
		     bool initial)
{
  sdc_mg_reset(&mg);
  // for (unsigned int i=0; i<=finest_level; i++) {
  //   if (sweepers[i] != NULL) {
  //     printf("%d %p %p\n", i, sweepers[i], sweepers[i]->destroy);
  //     sweepers[i]->destroy(sweepers[i]);
  //     delete (mf_encap_t*) encaps[i]->ctx;
  //     delete encaps[i];
  //     sweepers[i] = NULL;
  //   }
  // }

  Amr::regrid(lbase, time, initial);

  for (int lev=0; lev<=finest_level; lev++) {
    encaps[lev]   = build_encap(lev);
    sweepers[lev] = sdc_level_bld(*this, lev);
    sdc_mg_add_level(&mg, sweepers[lev], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  }
  sdc_mg_setup(&mg);
  sdc_mg_allocate(&mg);
  cout << "done" << endl;
}

SDCAmr::SDCAmr (sdc_level_bld_f bld)
{
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 22;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 3;

  sdc_level_bld = bld;
  sdc_mg_build(&mg, max_level);

  sweepers.resize(max_level);
  encaps.resize(max_level);

  for (unsigned int i=0; i<max_level; i++)
    sweepers[i] = NULL;
}


SDCAmr::~SDCAmr()
{
  sdc_mg_destroy(&mg);
}
