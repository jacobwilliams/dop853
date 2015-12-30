! * * * * * * * * * * * * * * * * * * * * * * * * *
! --- driver for dop853 on van der pol's equation
! * * * * * * * * * * * * * * * * * * * * * * * * *

program dr_dop853
    use dop853_module
    use dop853_constants

      implicit none

      real(wp) :: atol , rpar , rtol , tol , work , x ,     &
           & xend , y
      integer i , idid , iout , ipar , itol , iwork , j , liwork ,      &
            & lwork , n , ndgl , nrd
      parameter (ndgl=2,nrd=2)
      parameter (lwork=11*ndgl+8*nrd+21,liwork=nrd+21)
      dimension y(ndgl) , work(lwork) , iwork(liwork)
      dimension rtol(1) , atol(1) , rpar(1) , ipar(1)

! --- dimension of the system
      n = 2
      rpar = 1.0d-3
! --- output routine (and dense output) is used during integration
      iout = 3
! --- initial values
      x = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
! --- endpoint of integration
      xend = 2.0d0
! --- required (relative) tolerance
      tol = 1.0d-9
      itol = 0
      rtol = tol
      atol = tol
! --- default values for parameters
      do i = 1 , 10
         iwork(i) = 0
         work(i) = 0.d0
      enddo
      iwork(5) = n
      iwork(4) = 1
! --- call of the subroutine dopri8
      call dop853(n,fvpol,x,y,xend,rtol,atol,itol,solout,iout,work,     &
                & lwork,iwork,liwork,idid)
! --- print final solution
      write (6,99001) x , y(1) , y(2)
99001 format (1x,'x =',f5.2,'    y =',2e18.10)
! --- print statistics
      write (6,99002) tol
99002 format ('       tol=',d8.2)
      write (6,99003) (iwork(j),j=17,20)
99003 format (' fcn=',i5,' step=',i4,' accpt=',i4,' rejct=',i3)

      contains

      subroutine solout(nr,xold,x,y,n,con,icomp,nd,irtrn,xout)
          ! --- prints solution at equidistant output-points
          ! --- by using "contd8", the continuous collocation solution

      implicit none

      real(wp) :: con , x , xold , xout , y
      integer icomp , irtrn , n , nd , nr
      dimension y(n) , con(8*nd) , icomp(nd)

      if ( nr==1 ) then
         write (6,99001) x , y(1) , y(2) , nr - 1
         xout = 0.1d0
      else
 50      if ( x>=xout ) then
            write (6,99001) xout , contd8(1,xout,con,icomp,nd) ,        &
                          & contd8(2,xout,con,icomp,nd) , nr - 1
            xout = xout + 0.1d0
            goto 50
         endif
      endif

99001 format (1x,'x =',f5.2,'    y =',2e18.10,'    nstep =',i4)

      end subroutine solout

      subroutine fvpol(n,x,y,f)
          ! --- right-hand side of van der pol's equation

      implicit none

      integer,intent(in) :: n
      double precision,intent(in) :: x
      double precision,dimension(n),intent(in) :: y
      double precision,dimension(n),intent(out) :: f

      real(wp) :: eps

      eps = rpar(1)
      f(1) = y(2)
      f(2) = ((1-y(1)**2)*y(2)-y(1))/eps

      end subroutine fvpol

end program dr_dop853
