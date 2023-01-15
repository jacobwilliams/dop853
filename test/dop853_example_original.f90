!*****************************************************************************************
!>
!  Original Hairer test case.

    program dop853_example_original

    use dop853_module, wp => dop853_wp
    use iso_fortran_env, only: output_unit

    implicit none

    integer,parameter               :: n     = 2                !! dimension of the system
    real(wp),parameter              :: tol   = 1.0e-9_wp        !! integration tolerance
    real(wp),parameter              :: x0    = 0.0_wp           !! initial `x` value
    real(wp),parameter              :: xf    = 2.0_wp           !! endpoint of integration
    real(wp),dimension(n),parameter :: y0    = [2.0_wp,0.0_wp]  !! initial `y` value

    type(dop853_class) :: prop
    real(wp),dimension(n) :: y
    real(wp),dimension(1) :: rtol,atol
    real(wp) :: x
    integer :: idid
    logical :: status_ok

    x = x0
    y = y0
    rtol = tol ! set tolerances
    atol = tol !
    call prop%initialize(   fcn       = fvpol,    &
                            nstiff    = 1, &
                            n         = n,        &
                            status_ok = status_ok )
    if (.not. status_ok) error stop 'initialization error'

    call prop%integrate(x,y,xf,rtol,atol,iout=0,idid=idid)

    write (output_unit,'(1X,A,F6.2,A,2E18.10)') &
                       'x =',x ,'    y =',y(1),y(2)

contains

    subroutine fvpol(me,x,y,f)

    !! Right-hand side of van der Pol's equation

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: x
    real(wp),dimension(:),intent(in)  :: y
    real(wp),dimension(:),intent(out) :: f

    real(wp),parameter :: eps = 1.0e-3_wp

    f(1) = y(2)
    f(2) = ((1-y(1)**2)*y(2)-y(1))/eps

    end subroutine fvpol

end program dop853_example_original
