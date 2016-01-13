# dop853

This is a modern Fortran (2003/2008) implementation of Hairer's DOP853 ODE solver. The original FORTRAN 77 code has been extensively refactored, and is now object-oriented and thread-safe, with an easy-to-use class interface.  DOP853 is an explicit Runge-Kutta method of order 8(5,3) due to Dormand & Prince (with stepsize control and dense output).

This project is hosted on [GitHub](https://github.com/jacobwilliams/dop853).

## Example

Basic use of the solver is shown here. The main methods in the `dop853_class` are `initialize()` and `integrate()`.

```fortran
    program dop853_example

    use dop853_module
    use dop853_constants

    implicit none

    integer,parameter               :: n   = 2               !! dimension of the system
    real(wp),parameter              :: tol = 1.0e-12_wp      !! integration tolerance
    real(wp),parameter              :: x0  = 0.0_wp          !! initial x value
    real(wp),parameter              :: xf  = 100.0_wp        !! endpoint of integration
    real(wp),dimension(n),parameter :: y0  = [0.0_wp,0.1_wp] !! initial y value

    type(dop853_class)    :: prop
    real(wp),dimension(n) :: y
    real(wp),dimension(1) :: rtol,atol
    real(wp)              :: x
    integer               :: idid
    logical               :: status_ok

    x    = x0   ! initial conditions
    y    = y0   !
    rtol = tol  ! set tolerances
    atol = tol  !

    !initialize the integrator:
    call prop%initialize(fcn=fvpol,n=n,status_ok=status_ok)
    if (.not. status_ok) error stop 'initialization error'

    !now, perform the integration:
    call prop%integrate(x,y,xf,rtol,atol,iout=0,idid=idid)

    !print solution:
    write (output_unit,'(1X,A,F6.2,A,2E18.10)') &
                'x =',x ,' y =',y(1),y(2)

contains

    subroutine fvpol(me,x,y,f)
    !! Right-hand side of van der Pol's equation

    implicit none

    class(dop853_class),intent(inout) :: me
    real(wp),intent(in)               :: x
    real(wp),dimension(:),intent(in)  :: y
    real(wp),dimension(:),intent(out) :: f

    real(wp),parameter :: mu  = 0.2_wp

    f(1) = y(2)
    f(2) = mu*(1.0_wp-y(1)**2)*y(2) - y(1)

    end subroutine fvpol

end program dop853_example
```

The result is:

```
x =100.00 y = -0.1360372426E+01  0.1325538438E+01
```

For dense output, see the example in the `src/tests` directory.

# Documentation

The latest API documentation for the `master` branch can be found [here](http://jacobwilliams.github.io/dop853/). This was generated from the source code using [FORD](https://github.com/cmacmackin/ford) (note that the `build.sh` script will also generate these files if FORD is installed).

## References

1. E. Hairer, S.P. Norsett and G. Wanner, "[Solving ordinary
   Differential Equations I. Nonstiff Problems](http://www.unige.ch/~hairer/books.html)", 2nd edition.
   Springer Series in Computational Mathematics,
   Springer-Verlag (1993).
2. Ernst Hairer's website: [Fortran and Matlab Codes](http://www.unige.ch/~hairer/software.html)

## License

* [Original license for Hairer's codes](http://www.unige.ch/~hairer/prog/licence.txt).
* The updates are released under a [similar BSD-style license](https://raw.githubusercontent.com/jacobwilliams/dop853/master/LICENSE).
