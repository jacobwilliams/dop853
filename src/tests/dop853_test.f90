!*****************************************************************************************
!>
!  Driver for [[dop853]] on van der Pol's equation.
!
!### See also
!  * Based on [dr_dop853.f](http://www.unige.ch/~hairer/prog/nonstiff/dr_dop853.f)
!
!@note This requires [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran).

    program dop853_test

    use dop853_module
    use dop853_constants
    use iso_fortran_env, only: output_unit
    use pyplot_module

    implicit none

    integer,parameter               :: n     = 2                !! dimension of the system
    integer,dimension(n),parameter  :: icomp = [1,2]            !! indices where we need dense output
    integer,parameter               :: iout  = 3                !! output routine (and dense output) is used during integration
    real(wp),parameter              :: tol   = 1.0e-12_wp       !! required (relative) tolerance
    real(wp),parameter              :: x0    = 0.0_wp           !! initial values
    real(wp),parameter              :: xf    = 100.0_wp         !! endpoint of integration
    real(wp),dimension(n),parameter :: y0    = [0.0_wp,0.1_wp]  !! initial values
    real(wp),parameter              :: dx    = 0.01_wp          !! time step for dense output

    type(dop853_class) :: prop
    real(wp),dimension(n) :: y
    real(wp),dimension(1) :: rtol,atol
    real(wp) :: x,xend
    integer :: i,idid,j,nfcn,nstep,naccpt,nrejct
    logical :: status_ok
    type(pyplot) :: plt
    real(wp),dimension(:),allocatable :: t_vec, y_vec, yp_vec

    x = x0
    y = y0
    xend = xf
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

        t_vec  = [t_vec,x]        !last point
        y_vec  = [y_vec,y(1)]
        yp_vec = [yp_vec,y(2)]

        ! print final solution
        write (output_unit,'(1X,A,F6.2,A,2E18.10)') &
                'x =',x ,'    y =',y(1),y(2)

        ! print statistics
        write (output_unit,'(A,D8.2)') '       tol=',tol
        write (output_unit,'(A,I5,A,I4,A,I4,A,I3)') &
                ' fcn=',nfcn,' step=',nstep,' accpt=',naccpt,' rejct=',nrejct

        ! plot:
        call plt%initialize(grid=.true.,xlabel='y(x)',&
                          ylabel='y''(x)',axis_equal=.true.,&
                          title='van der Pol''s Equation ($\mu = 0.2$)',legend=.true.)
        call plt%add_plot(y_vec,yp_vec,label='Forward',&
                          linestyle='r-',linewidth=2)
        call plt%savefig('dop853_forward.png')
        call plt%destroy()

        deallocate(t_vec )
        deallocate(y_vec )
        deallocate(yp_vec)

        !--------------------------------------------------

        write(*,*) ''
        write(*,*) 'backwards test'
        write(*,*) ''
        call prop%destroy()
        call prop%initialize(fcn=fvpol,solout=solout2,status_ok = status_ok)
        call prop%integrate(n,x,y,x0,rtol,atol,iout=1,idid=idid)

        write(*,*) ''
        write(*,*) 'error:', norm2(y-y0)
        write(*,*) ''

        ! plot:
        call plt%initialize(grid=.true.,xlabel='y(x)',&
                        ylabel='y''(x)',axis_equal=.true.,&
                        title='van der Pol''s Equation ($\mu = 0.2$)',legend=.true.)

        call plt%add_plot(y_vec,yp_vec,label='Backward',&
                        linestyle='r-',linewidth=2)

        call plt%savefig('dop853_backward.png')
        call plt%destroy()

    else
      write(output_unit,'(A)') 'error calling INITIALIZE.'
    end if

    contains
!*****************************************************************************************

!*******************************************************************************
!>
!  Prints solution at equidistant output-points
!  by using [[contd8]], the continuous collocation solution.
!
!@note This routine uses Fortran 2008 LHS automatic allocations.

    subroutine solout(me,nr,xold,x,y,n,irtrn,xout)

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
        write (output_unit,'(1X,A,F6.2,A,2E18.10,A,I4)') &
                'x =',x,&
                '    y =',y(1),y(2),&
                '    nstep =',nr - 1
        xout = dx
    else
        do
            if ( x<xout ) exit
            write (output_unit,'(1X,A,F6.2,A,2E18.10,A,I4)') &
                     'x =',xout,&
                     '    y =',&
                     prop%contd8(1,xout),&
                     prop%contd8(2,xout),&
                     '    nstep =',nr - 1
            xout = xout + dx
        end do
    end if

    if (allocated(t_vec)) then !fortran 2008 lhs allocations
        t_vec  = [t_vec,x]
        y_vec  = [y_vec,y(1)]
        yp_vec = [yp_vec,y(2)]
    else
        t_vec  = [x]
        y_vec  = [y(1)]
        yp_vec = [y(2)]
    end if

    end subroutine solout
!*******************************************************************************

!*******************************************************************************
!>
!  Prints a normal step from [[dop853]] (for `iout=1`).
!
!@note This routine uses Fortran 2008 LHS automatic allocations.

    subroutine solout2(me,nr,xold,x,y,n,irtrn,xout)

    implicit none

    class(dop853_class),intent(inout) :: me
    integer,intent(in)                :: nr
    real(wp),intent(in)               :: xold
    real(wp),intent(in)               :: x
    integer,intent(in)                :: n
    real(wp),dimension(n),intent(in)  :: y
    integer,intent(inout)             :: irtrn
    real(wp),intent(out)              :: xout   !! not used for `iout=1`.

    write (output_unit,'(1X,A,F6.2,A,2E18.10,A,I4)') &
                    'x =',x,&
                    '    y =',y(1),y(2),&
                    '    nstep =',nr - 1

    if (allocated(t_vec)) then
        t_vec  = [t_vec,x]
        y_vec  = [y_vec,y(1)]
        yp_vec = [yp_vec,y(2)]
    else
        t_vec  = [x]
        y_vec  = [y(1)]
        yp_vec = [y(2)]
    end if

    end subroutine solout2
!*******************************************************************************

!*******************************************************************************
!>
! Right-hand side of van der Pol's equation:
! \( y'' = \mu(1-y^2)y' - y \).
!
!### Reference
! * Weisstein, Eric W. "[van der Pol Equation](http://mathworld.wolfram.com/vanderPolEquation.html)."
!   From MathWorld--A Wolfram Web Resource.
!
!@note The [original code](http://www.unige.ch/~hairer/prog/nonstiff/dr_dop853.f)
!      had a slightly different equation.

    subroutine fvpol(me,n,x,y,f)

    implicit none

    class(dop853_class),intent(inout) :: me
    integer,intent(in)                :: n
    real(wp),intent(in)               :: x
    real(wp),dimension(n),intent(in)  :: y  !! \( [y, y' ] \)
    real(wp),dimension(n),intent(out) :: f  !! \( [y', y'' ] \)

    real(wp),parameter :: mu  = 0.2_wp  !! \( \mu \) in van der Pol's equation.

    f(1) = y(2)
    f(2) = mu*(1.0_wp-y(1)**2)*y(2) - y(1)

    end subroutine fvpol
!*******************************************************************************

    end program dop853_test
!*****************************************************************************************
