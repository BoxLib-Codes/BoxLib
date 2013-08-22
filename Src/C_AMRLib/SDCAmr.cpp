
#include <SDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>

#include <stdio.h>

using namespace std;


BEGIN_EXTERN_C

void mlsdc_amr_interpolate(void *F, void *G, void *ctxF, void *ctxG)
{
  MultiFab& UF = *((MultiFab*) F);
  MultiFab& UG = *((MultiFab*) G);
  AmrLevel& levelF = *((AmrLevel*) ctxF);
  AmrLevel& levelG = *((AmrLevel*) ctxG);

}

void mlsdc_amr_restrict(void *F, void *G, void *ctxF, void *ctxG)
{
  MultiFab& UF = *((MultiFab*) F);
  MultiFab& UG = *((MultiFab*) G);
  AmrLevel& levelF = *((AmrLevel*) ctxF);
  AmrLevel& levelG = *((AmrLevel*) ctxG);

}

END_EXTERN_C


void SDCAmr::timeStep (int  level,
                         Real time,
                         int  iteration,
                         int  niter,
                         Real stop_time)
{
  if (level != 0) {
    cout << "NOT ON LEVEL 0" << endl;
  }

  // set initial conditions
  for (int lev=0; lev<=finest_level; lev++) {
    AmrLevel& level = getLevel(lev);
    const DescriptorList& dl = level.get_desc_lst();
    for (int st=0; st<dl.size(); st++) {
      MultiFab& Unew = level.get_new_data(st);
      // XXX: this assumes that the encapsulation is a multifab
      MultiFab& U0   = *((MultiFab*)mg.sweepers[lev]->nset->Q[0]);
      U0.copy(Unew);
    }
  }

  // spread and iterate
  sdc_mg_spread(&mg, time, dtLevel(0));
  for (int k=0; k<max_iters; k++) {
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
      std::cout << "MLSDC iteration: " << k << std::endl;
    sdc_mg_sweep(&mg, time, dt_level[0], 0);
  }

  // grab final solution
  for (int lev=0; lev<=finest_level; lev++) {
    AmrLevel& level = getLevel(lev);
    const DescriptorList& dl = level.get_desc_lst();
    for (int st=0; st<dl.size(); st++) {
      MultiFab& Unew = level.get_new_data(st);
      // XXX: this assumes that the encapsulation is a multifab
      int nnodes = mg.sweepers[lev]->nset->nnodes;
      MultiFab& Uend = *((MultiFab*)mg.sweepers[lev]->nset->Q[nnodes-1]);
      Unew.copy(Uend);
    }
  }
}

void SDCAmr::rebuild_mlsdc()
{
  // reset previous and clear sweepers etc
  sdc_mg_reset(&mg);
  for (unsigned int i=0; i<=max_level; i++) {
    if (sweepers[i] != NULL) {
      sweepers[i]->destroy(sweepers[i]);
      delete (mf_encap_t*) encaps[i]->ctx;
      delete encaps[i];
      sweepers[i] = NULL;
    }
  }

  // rebuild
  for (int lev=0; lev<=finest_level; lev++) {
    encaps[lev]   = build_encap(lev);
    sweepers[lev] = sdc_level_bld(*this, lev);
    sdc_mg_add_level(&mg, sweepers[lev], mlsdc_amr_interpolate, mlsdc_amr_restrict);
  }
  sdc_mg_setup(&mg);
  sdc_mg_allocate(&mg);

  if (verbose > 0 && ParallelDescriptor::IOProcessor())
    std::cout << "Rebuilt MLSDC: " << mg.nlevels << std::endl;
}

// void SDCAmr::regrid (int  lbase,
// 		     Real time,
// 		     bool initial)
// {
//   Amr::regrid(lbase, time, initial);
//   rebuild_mlsdc();
// }

SDCAmr::SDCAmr (sdc_level_bld_f bld)
{
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 22;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 3;

  sdc_level_bld = bld;
  sdc_mg_build(&mg, max_level+1);

  sweepers.resize(max_level+1);
  encaps.resize(max_level+1);

  for (unsigned int i=0; i<=max_level; i++)
    sweepers[i] = NULL;
}


SDCAmr::~SDCAmr()
{
  sdc_mg_destroy(&mg);
}
