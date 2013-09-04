
#include <SDCAmr.H>
#include <MultiFab.H>
#include <ParmParse.H>
#include <StateDescriptor.H>
#include <AmrLevel.H>
#include <Interpolater.H>
#include <FabArray.H>
#include <stdio.h>

using namespace std;

BEGIN_EXTERN_C

/*
 * Spatial interpolation is done by...
 */
void mlsdc_amr_interpolate(void *F, void *G, void *vctxF, void *vctxG)
{
  sdc_level_ctx* ctxF = (sdc_level_ctx*) vctxF;
  sdc_level_ctx* ctxG = (sdc_level_ctx*) vctxG;
  MultiFab& UF        = *((MultiFab*) F);
  MultiFab& UG        = *((MultiFab*) G);
  Amr& amr            = *ctxF->amr;
  AmrLevel& levelF    = amr.getLevel(ctxF->level);
  AmrLevel& levelG    = amr.getLevel(ctxG->level);

  const DescriptorList& dl = levelF.get_desc_lst();
  for (int st=0; st<dl.size(); st++) {
    int ncomp = dl[st].nComp();

    const Array<BCRec>& bcs = dl[st].getBCs();

    Array<BCRec> bcr(ncomp);

    const IntVect& fine_ratio = amr.refRatio(ctxG->level);
    Interpolater& mapper      = *dl[st].interp();

    const Box& fine_domain = levelF.get_state_data(st).getDomain();

    for (MFIter mfi(UF); mfi.isValid(); ++mfi) { 

      const Box& fine_box = mfi.validbox();

      FArrayBox fine(fine_box);
      FArrayBox crse(mapper.CoarseBox(fine_box, fine_ratio));
      for (int j=0; j<UG.size(); j++) crse.copy(UG[j]);

      cout << crse.box() << endl;
      cout << fine.box() << endl;
      cout << fine_box << endl;
      cout << ncomp << endl;
      // cout << levelG.Geom() << endl;
      // cout << levelF.Geom() << endl;

      BoxLib::setBC(fine_box, fine_domain, st, 0, ncomp, bcs, bcr);

      const Geometry& crse_geom = levelG.Geom();
      const Geometry& fine_geom = levelF.Geom();


      cout << "interpolating..." << endl;
      mapper.interp(crse, 0, fine, 0, ncomp,
                    fine_box, fine_ratio, crse_geom, fine_geom, bcr, st, st);
      cout << "interpolating... done." << endl;
    }
  }

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

  // set intial conditions...
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
  sdc_mg_spread(&mg, time, dtLevel(0), 0);
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
  for (unsigned int lev=0; lev<=max_level; lev++) {
    if (sweepers[lev] != NULL) {
      delete sweepers[lev]->nset->ctx;
      sweepers[lev]->destroy(sweepers[lev]);
      delete (mf_encap*) encaps[lev]->ctx;
      delete encaps[lev];
      sweepers[lev] = NULL;
    }
  }

  // rebuild
  for (int lev=0; lev<=finest_level; lev++) {
    sdc_level_ctx* ctx = new sdc_level_ctx;
    ctx->amr   = this;
    ctx->level = lev;
    encaps[lev]   = build_encap(lev);
    sweepers[lev] = sdc_sweeper_bld(lev);
    sweepers[lev]->nset->ctx   = ctx; // XXX: need to free this...
    sweepers[lev]->nset->encap = encaps[lev];
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

SDCAmr::SDCAmr (sdc_sweeper_bld_f bld)
{
  ParmParse ppsdc("mlsdc");
  if (!ppsdc.query("max_iters", max_iters)) max_iters = 22;
  if (!ppsdc.query("max_trefs", max_trefs)) max_trefs = 3;

  sdc_sweeper_bld = bld;
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
