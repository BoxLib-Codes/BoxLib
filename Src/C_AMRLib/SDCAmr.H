
#ifndef _SDCAMR_H_
#define _SDCAMR_H_

#define BEGIN_EXTERN_C extern "C" {
#define END_EXTERN_C   }

#include <MultiFab.H>
#include <Amr.H>
#include <sdc.h>

using namespace std;

/*
 * MultiFab encapsulation context for SDCLib.
 */

typedef struct {
  const BoxArray* ba;
  int       ncomp;
  int       ngrow;
} mf_encap;

typedef struct {
  Amr *amr;
  int level;
} sdc_level_ctx;

/*
 * Multigrid SDC+AMR class.
 */

class SDCAmr;

typedef sdc_sweeper* (*sdc_sweeper_bld_f)(int level);

class SDCAmr : public Amr
{
public:

  int max_iters;                // maximum number of sdc iterations
  int max_trefs;                // maximum number of time refinements

  sdc_mg mg;
  sdc_sweeper_bld_f sdc_sweeper_bld;
  vector<sdc_sweeper*> sweepers;
  vector<sdc_encap*>   encaps;

  SDCAmr(sdc_sweeper_bld_f);
  ~SDCAmr();

  sdc_encap* build_encap(int lev);
  void rebuild_mlsdc();

  virtual void timeStep(int  level,
                        Real time,
                        int  iteration,
                        int  niter,
                        Real stop_time);

  virtual void regrid (int  lbase,
		       Real time,
		       bool initial = false) {
    Amr::regrid(lbase, time, initial);
    rebuild_mlsdc();
  }

  virtual void initialInit (Real strt_time, Real stop_time) {
    Amr::initialInit(strt_time, stop_time);
    //    rebuild_mlsdc();
  }


};

#endif
