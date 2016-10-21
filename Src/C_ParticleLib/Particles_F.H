#ifndef _Particles_F_H_

#define _Particles_F_H_

#include <BLFort.H>

#ifdef USE_PARTICLES

extern "C"
{

BL_FORT_PROC_DECL(PART_SUMMASSDOWN,part_summassdown)
    (BL_FORT_FAB_ARG(crse_fab),
     const BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMMASSUP,part_summassup)
    (const BL_FORT_FAB_ARG(crse_fab),
     BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMDENSUP,part_sumdensup)
    (const BL_FORT_FAB_ARG(crse_fab),
     BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

#if (BL_SPACEDIM==2)
    BL_FORT_PROC_DECL(BICUBICINTERP,bicubicinterp)
    (const Real* x, const Real* y, const Real* fab,
     const int* lo, const int* hi,
     const Real* plo, const Real* phi, const Real* dx, Real* z);
#endif

#if (BL_SPRAY_PARTICLES)
    BL_FORT_PROC_DECL(SPRAY_PART_SRC,spray_part_src)
        (Real* pstate, Real* fld_at_part, Real* pstate_src,
         int* n_part, int* n_state, int* n_fld_comp);

    BL_FORT_PROC_DECL(SPRAY_PART_SRC_TO_FLD_SRC,spray_part_src_to_fld_src)
        (Real* pstate, Real* fld_at_part, Real* pstate_src, Real* fld_src,
         int* n_part, int* n_state, int* n_fld_comp);
#endif
};

#endif

#endif /*_Particles_F_H_*/
