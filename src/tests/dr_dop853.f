C * * * * * * * * * * * * * * * * * * * * * * * * *
C --- DRIVER FOR DOPRI5 ON VAN DER POL'S EQUATION
C * * * * * * * * * * * * * * * * * * * * * * * * *
cfeh dr_dop853 dop853
c        include 'dop853.f'
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NDGL=2,NRD=2)
        PARAMETER (LWORK=11*NDGL+8*NRD+21,LIWORK=NRD+21)
        DIMENSION Y(NDGL),WORK(LWORK),IWORK(LIWORK)
        EXTERNAL FVPOL,SOLOUT
	  DIMENSION RTOL(1),ATOL(1),RPAR(1),IPAR(1)
C --- DIMENSION OF THE SYSTEM
        N=2
        RPAR=1.0D-3
C --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
        IOUT=3
C --- INITIAL VALUES
        X=0.0D0
        Y(1)=2.0D0
        Y(2)=0.0D0
C --- ENDPOINT OF INTEGRATION
        XEND=2.0D0
C --- REQUIRED (RELATIVE) TOLERANCE
        TOL=1.0D-9
        ITOL=0
        RTOL=TOL
        ATOL=TOL
C --- DEFAULT VALUES FOR PARAMETERS
        DO 10 I=1,10
        IWORK(I)=0
  10    WORK(I)=0.D0
        IWORK(5)=N
        IWORK(4)=1
C --- CALL OF THE SUBROUTINE DOPRI8
        CALL DOP853(N,FVPOL,X,Y,XEND,
     &                  RTOL,ATOL,ITOL,
     &                  SOLOUT,IOUT,
     &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
C --- PRINT FINAL SOLUTION
        WRITE (6,99) X,Y(1),Y(2)
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10)
C --- PRINT STATISTICS
        WRITE (6,90) TOL
 90     FORMAT('       tol=',D8.2)
        WRITE (6,91) (IWORK(J),J=17,20)
 91     FORMAT(' fcn=',I5,' step=',I4,' accpt=',I4,' rejct=',I3)
        STOP
        END
C
        SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
     &                                          RPAR,IPAR,IRTRN,XOUT)
C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
C --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),CON(8*ND),ICOMP(ND)
        IF (NR.EQ.1) THEN
           WRITE (6,99) X,Y(1),Y(2),NR-1
           XOUT=0.1D0
        ELSE
 10        CONTINUE
           IF (X.GE.XOUT) THEN
              WRITE (6,99) XOUT,CONTD8(1,XOUT,CON,ICOMP,ND),
     &                     CONTD8(2,XOUT,CON,ICOMP,ND),NR-1
              XOUT=XOUT+0.1D0
              GOTO 10
           END IF
        END IF
 99     FORMAT(1X,'X =',F5.2,'    Y =',2E18.10,'    NSTEP =',I4)
        RETURN
        END
C
        SUBROUTINE FVPOL(N,X,Y,F,RPAR,IPAR)
C --- RIGHT-HAND SIDE OF VAN DER POL'S EQUATION
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION Y(N),F(N)
        EPS=RPAR
        F(1)=Y(2)
        F(2)=((1-Y(1)**2)*Y(2)-Y(1))/EPS
        RETURN
        END
