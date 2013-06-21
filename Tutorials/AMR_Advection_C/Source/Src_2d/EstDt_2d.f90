     subroutine estdt(u,u_l1,u_l2,u_h1,u_h2,lo,hi,dx,dt)

     use meth_params_module, only : NVAR, UX, UY, UFA

     implicit none

     integer          :: u_l1,u_l2,u_h1,u_h2
     integer          :: lo(2), hi(2)
     double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
     double precision :: dx(2),dt

     double precision :: xvel,yvel,dt1,dt2,dt_loc
     double precision :: s_max
     integer          :: i,j

     s_max = -1.d20
     do j = lo(2),hi(2)
         do i = lo(1),hi(1)
             s_max = max(s_max,u(i,j,UFA))
         end do
     end do

     dt_loc = 1.d200
     do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            xvel = u(i,j,UX)
            yvel = u(i,j,UY)

            dt1 = dx(1)/abs(xvel)
            dt2 = dx(2)/abs(yvel)

            dt_loc = min(dt_loc,dt1,dt2)

         enddo
      enddo

      dt = min(dt,dt_loc)

      end subroutine estdt
