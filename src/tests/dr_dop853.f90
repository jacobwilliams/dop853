!*****************************************************************************************
!>
!  Driver for [[dop853]] on van der Pol's equation.
!
    program dr_dop853

    use dop853_module
    use dop853_constants
    use iso_fortran_env, only: output_unit

    implicit none

    integer :: i  !! counter

    integer,parameter              :: n     = 2            !! dimension of the system
    integer,dimension(n),parameter :: icomp = [(i,i=1,n)]  !! indices where we need dense output
    integer,parameter              :: iout  = 3            !! output routine (and dense output) is used during integration
    real(wp),parameter             :: tol   = 1.0e-9_wp    !! required (relative) tolerance

    type(dop853_class) :: prop
    real(wp),dimension(n) :: y
    real(wp),dimension(1) :: rtol,atol
    real(wp) :: x,xend
    integer :: idid,j,nfcn,nstep,naccpt,nrejct
    logical :: status_ok

    x    = 0.0_wp           ! initial values
    y    = [2.0_wp,0.0_wp]  !
    xend = 2.0_wp           ! endpoint of integration
    rtol = tol
    atol = tol
    call prop%initialize(   fcn       = fvpol,    &
                            solout    = solout,   &
                            icomp     = icomp,    &
                            nstiff    = 1,        &
                            status_ok = status_ok )
    !all other parameters use defaults

    if (status_ok) then

      call prop%integrate(n,x,y,xend,rtol,atol,iout,idid)
      call prop%info(nfcn,nstep,naccpt,nrejct)

      ! print final solution
      write (output_unit,'(1X,A,F5.2,A,2E18.10)') &
                'x =',x ,'    y =',y(1),y(2)

      ! print statistics
      write (output_unit,'(A,D8.2)') '       tol=',tol
      write (output_unit,'(A,I5,A,I4,A,I4,A,I3)') &
                ' fcn=',nfcn,' step=',nstep,' accpt=',naccpt,' rejct=',nrejct

    else
      write(output_unit,'(A)') 'error calling INITIALIZE.'
    end if

    contains
!*****************************************************************************************

!*******************************************************************************
    subroutine solout(me,nr,xold,x,y,n,irtrn,xout)

    !! prints solution at equidistant output-points
    !! by using [[contd8]], the continuous collocation solution.

    implicit none

    class(dop853_class),intent(inout) :: me
    integer,intent(in)                :: nr
    real(wp),intent(in)               :: xold
    real(wp),intent(in)               :: x
    integer,intent(in)                :: n
    real(wp),dimension(n),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(wp),intent(out)              :: xout   !! the point where we want the next output reported

    if ( nr==1 ) then
        write (output_unit,'(1X,A,F5.2,A,2E18.10,A,I4)') &
                'x =',x,&
                '    y =',y(1),y(2),&
                '    nstep =',nr - 1
        xout = 0.1_wp
    else
        do
            if ( x<xout ) exit
            write (output_unit,'(1X,A,F5.2,A,2E18.10,A,I4)') &
                     'x =',xout,&
                     '    y =',&
                     prop%contd8(1,xout),&
                     prop%contd8(2,xout),&
                     '    nstep =',nr - 1
            xout = xout + 0.1_wp
        end do
    end if

    end subroutine solout
!*******************************************************************************

!*******************************************************************************
    subroutine fvpol(me,n,x,y,f)

    !! right-hand side of van der Pol's equation

    implicit none

    class(dop853_class),intent(inout)         :: me
    integer,intent(in)                        :: n
    double precision,intent(in)               :: x
    double precision,dimension(n),intent(in)  :: y
    double precision,dimension(n),intent(out) :: f

    real(wp),parameter :: eps = 1.0e-3_wp

    f(1) = y(2)
    f(2) = ((1.0_wp-y(1)**2)*y(2)-y(1))/eps

    end subroutine fvpol
!*******************************************************************************

    end program dr_dop853
!*****************************************************************************************
