
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
} mf_encap_t;

/*
 * Multigrid SDC+AMR class.
 */

class SDCAmr;

typedef sdc_sweeper_t* (*sdc_level_bld_f)(SDCAmr& amr, int level);

class SDCAmr : public Amr
{
public:

  int max_iters;                // maximum number of sdc iterations
  int max_trefs;                // maximum number of time refinements

  sdc_mg_t mg;
  sdc_level_bld_f sdc_level_bld;
  vector<sdc_sweeper_t*> sweepers;
  vector<sdc_encap_t*>   encaps;

  SDCAmr(sdc_level_bld_f);
  ~SDCAmr();

  sdc_encap_t* build_encap(int lev);
  virtual void timeStep(int  level,
                        Real time,
                        int  iteration,
                        int  niter,
                        Real stop_time);
  virtual void regrid (int  lbase,
		       Real time,
		       bool initial = false);

};

#endif
