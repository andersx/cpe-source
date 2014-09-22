************************************************************************
C     Eigenvalue-solver for IBM RS6000/SP2 machines with esslp2 library
C
C     Parameters: 
C 
C       NA      (I) :  Dimension of A 
C       NB      (I) :  Dimension of B 
C       N       (I) :  Dimension of Problem  
C       A       (I) :  Matrix A 
C               (O) :  Eigenvector matrix  
C       B       (I) :  Matrix B 
C       EW      (O) :  Eigenvalues 
C       H       (-) :  Dummy vector 
C       AUX     (-) :  Auxiliary vector 
C       IEV     (I) :  0: No eigenvectors  
C       IER     (O) :  Dummy variable  
C 
C ********************************************************************** 
      SUBROUTINE EWEVGE (NA,NB,N,A,B,EW,H,IEV,IORD,IER) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION  A(NA,N),B(NB,N),EW(N),H(N),AUX(3*NA)

      CALL DSYGV(IEV,'V','L',N,A,NA,B,NB,EW,AUX,3*NA,IER) 

      RETURN
      END
