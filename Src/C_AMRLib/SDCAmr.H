
#ifndef _SDCAMR_H_
#define _SDCAMR_H_

#define BEGIN_EXTERN_C extern "C" {
#define END_EXTERN_C   }

#include <MultiFab.H>
#include <Amr.H>
#include <sdc.h>


/*
 * MultiFab encapsulation context for SDCLib.
 */

typedef struct {
  BoxArray& ba;
  int       ncomp;
  int       ngrow;
} mf_encap_t;

/*
 * Multigrid SDC+AMR class.
 */

using namespace std;

class SDCLevelBld;

class SDCAmr : public Amr
{
public:

  int max_iters;                // maximum number of sdc iterations
  int max_trefs;                // maximum number of time refinements

  sdc_mg_t    mg;
  sdc_encap_t encap;
  mf_encap_t  encap_ctx;

  SDCAmr();
  ~SDCAmr();
  void InitSDCAmr();

  SDCLevelBld* sdc_level_bld;
  vector<sdc_sweeper_t*> sweepers;

  void build_encap();

  virtual void timeStep(int  level,
                        Real time,
                        int  iteration,
                        int  niter,
                        Real stop_time);


};

#endif
