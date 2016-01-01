program dr_dop853

    !! driver for [[dop853]] on van der pol's equation

    use dop853_module
    use dop853_constants
    use iso_fortran_env, only: output_unit

    implicit none

    integer :: i

    integer,parameter :: n = 2  !dimension of the system
    integer,dimension(n),parameter :: icomp = [(i,i=1,n)]

    type(dop853_class) :: prop
    real(wp) :: tol , x , xend, eps
    real(wp),dimension(n) :: y
    real(wp),dimension(1) :: rtol,atol
    integer :: idid,iout,itol,j,nfcn,nstep,naccpt,nrejct
    logical :: status_ok

    eps  = 1.0e-3_wp
    iout = 3          ! output routine (and dense output) is used during integration
    x    = 0.0_wp     ! initial values
    y(1) = 2.0_wp
    y(2) = 0.0_wp
    xend = 2.0_wp     ! endpoint of integration
    tol  = 1.0e-9_wp  ! required (relative) tolerance
    rtol = tol
    atol = tol
    call prop%initialize(ICOMP=icomp,NSTIFF=1,STATUS_OK=status_ok,hinitial=1.0_wp)
    !all other parameters use defaults
    
    if (status_ok) then

      call prop%integrate(n,fvpol,x,y,xend,rtol,atol,solout,iout,idid)
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

    subroutine solout(nr,xold,x,y,n,con,icomp,nd,irtrn,xout)

    !! prints solution at equidistant output-points
    !! by using `contd8`, the continuous collocation solution

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
                     'x =',xout,'    y =',&
                     prop%contd8(1,xout,con,icomp,nd),&
                     prop%contd8(2,xout,con,icomp,nd),&
                     '    nstep =',nr - 1
            xout = xout + 0.1_wp
        end do
    end if

    end subroutine solout

    subroutine fvpol(n,x,y,f)

    !! right-hand side of van der pol's equation

    implicit none

    integer,intent(in)                        :: n
    double precision,intent(in)               :: x
    double precision,dimension(n),intent(in)  :: y
    double precision,dimension(n),intent(out) :: f

    f(1) = y(2)
    f(2) = ((1.0_wp-y(1)**2)*y(2)-y(1))/eps

    end subroutine fvpol

end program dr_dop853
