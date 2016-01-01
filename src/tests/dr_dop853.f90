! * * * * * * * * * * * * * * * * * * * * * * * * *
! --- driver for dop853 on van der pol's equation
! * * * * * * * * * * * * * * * * * * * * * * * * *

program dr_dop853
    use dop853_module
    use dop853_constants
    use iso_fortran_env, only: output_unit

      implicit none

      type(dop853_class) :: prop

      integer :: i
      integer,parameter :: n = 2  !dimension of the system
      integer,dimension(n),parameter :: icomp = [(i,i=1,n)]

      real(wp) :: atol , rpar , rtol , tol , work , x ,     &
           & xend , y
      integer idid , iout , ipar , itol , iwork , j , liwork ,      &
            & lwork , ndgl , nrd
      parameter (ndgl=2,nrd=2)
      parameter (lwork=11*ndgl+8*nrd+21,liwork=nrd+21)
      dimension y(ndgl) , work(lwork) , iwork(liwork)
      dimension rtol(1) , atol(1) , rpar(1) , ipar(1)

      integer :: nfcn,nstep,naccpt,nrejct
      logical :: status_ok

! --- dimension of the system
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

      call prop%initialize(ICOMP=icomp,NSTIFF=1,STATUS_OK=status_ok)
      call prop%integrate(n,fvpol,x,y,xend,rtol,atol,solout,iout,idid)
      call prop%info(nfcn,nstep,naccpt,nrejct)

! --- print final solution
      write (6,99001) x , y(1) , y(2)
99001 format (1x,'x =',f5.2,'    y =',2e18.10)
! --- print statistics
      write (6,99002) tol
99002 format ('       tol=',d8.2)
      write (output_unit,'(A,I5,A,I4,A,I4,A,I3)') &
              ' fcn=',nfcn,' step=',nstep,' accpt=',naccpt,' rejct=',nrejct

      contains

      subroutine solout(nr,xold,x,y,n,con,icomp,nd,irtrn,xout)
          ! --- prints solution at equidistant output-points
          ! --- by using "contd8", the continuous collocation solution

      implicit none

      integer,intent(in)                  :: nr
      real(wp),intent(in)                 :: xold
      real(wp),intent(in)                 :: x
      integer,intent(in)                  :: n
      integer,intent(in)                  :: nd
      real(wp),dimension(n),intent(in)    :: y
      real(wp),dimension(8*nd),intent(in) :: con
      integer,dimension(nd),intent(in)    :: icomp
      integer,intent(inout)               :: irtrn
      real(wp),intent(out)                :: xout

     ! real(wp) :: con , x , xold , xout , y
     ! integer icomp , irtrn , n , nd , nr
     ! dimension y(n) , con(8*nd) , icomp(nd)

      if ( nr==1 ) then
         write (6,99001) x , y(1) , y(2) , nr - 1
         xout = 0.1d0
      else
 50      if ( x>=xout ) then
            write (6,99001) xout , prop%contd8(1,xout,con,icomp,nd) ,        &
                          & prop%contd8(2,xout,con,icomp,nd) , nr - 1
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
