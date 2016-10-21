! ::
! :: ----------------------------------------------------------
! :: Compute source term for particle state update
! :: ----------------------------------------------------------
! ::
      subroutine spray_part_src(pstate, fld_at_part, pstate_src,
     &                          n_part, n_state, n_fld_comp)
      implicit none

      integer          n_part, n_state, n_fld_comp
      double precision pstate(n_part,n_state)
      double precision pstate_src(n_part,n_state)
      double precision fld_at_part(n_part, n_fld_comp)

      return

      end subroutine spray_part_src

! ::
! :: ----------------------------------------------------------
! :: Compute source term for field date from particle sources
! :: ----------------------------------------------------------
! ::
      subroutine spray_part_src(pstate, fld_at_part, pstate_src,
     &                          fld_src,
     &                          n_part, n_state, n_fld_comp)
      implicit none

      integer          n_part, n_state, n_fld_comp
      double precision pstate(n_part,n_state)
      double precision pstate_src(n_part,n_state)
      double precision fld_at_part(n_part, n_fld_comp)
      double precision fld_src(n_part, n_fld_comp)


      return


      end subroutine spray_part_src
