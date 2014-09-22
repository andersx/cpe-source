C
C ************************************************************
C
C LOGCGR VERSION DIRK POREZAG AUGUST 1994
C
       SUBROUTINE LOGCGR(IUNIT,IWARN,ILGRET,IMODE,ILINE,ISTEP)
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*23 LINE
        LOGICAL ILGRET
        SAVE
C
        IF (IWARN.EQ.1) THEN
          WRITE(IUNIT,*) 'Desired accuracy is small, single ',
     &                   'precision numerics may cause problems'
        END IF
        IF (IWARN.EQ.2) THEN
          WRITE(IUNIT,*) 'Accuracy too low to work with ',
     &                   'energies, checking derivatives only'
        END IF
        IF (IWARN.EQ.3) THEN
          WRITE(IUNIT,*) 'Maximum in (a,c) detected in line ',
     &                   'minimization'
        END IF
        IF (IWARN.EQ.4) THEN
          WRITE(IUNIT,*) 'Maximum in (a,b) and minimum in (b,c) ',
     &                   'detected in line minimization'
        END IF
        IF (IWARN.EQ.5) THEN
          WRITE(IUNIT,*) 'Maximum in (b,c) and minimum in (a,b) ',
     &                   'detected in line minimization'
        END IF
        IF (ILGRET) RETURN
        IF (IMODE.EQ.0) THEN
         LINE='Error: Imode=0         '
        ELSE IF (IMODE.EQ.1) THEN
         LINE='Expanding interval     '
        ELSE IF (IMODE.EQ.2) THEN
         LINE='Quadratic interpolation'
        ELSE IF (IMODE.EQ.3) THEN
         LINE='Linear interpolation   '
        ELSE IF (IMODE.EQ.4) THEN
         LINE='Bisection              '
        END IF
        WRITE(IUNIT,100) ISTEP,ILINE,LINE
  100   FORMAT('Istep= ',I2,', Iline= ',I2,',  ',A23)
        RETURN
       END
C
C ************************************************************
C
C CGRAD VERSION DIRK POREZAG AUGUST 1994
C
C INPUT: MODE:    0 -> CGRAD STEP
C                 1 -> CLOSE FILES AND RETURN
C                 2 -> CLOSE FILES, DELETE CGRAD AND RETURN
C        NPAR:    NUMBER OF RELEVANT DEGREES OF FREEDOM
C        GTOL:    MAXIMAL FORCE FOR CONVERGENCE
C        ENERGY:  ENERGY FOR CURRENT STRUCTURE
C        Y:       COORDINATES FOR CURRENT STRUCTURE
C        G:       GRADIENT FOR CURRENT STRUCTURE
C
       SUBROUTINE CGRAD(MODE,NPAR,GTOL,ENERGY,Y,G)
       IMPLICIT REAL*8 (A-H,O-Z)
       INCLUDE 'maxima.inc'
       PARAMETER (MAXLINE=15)
       PARAMETER (MAX2=MAXLINE+2)
       DIMENSION Y(*),G(*)
       LOGICAL   FIRST,LOPF
       LOGICAL   IRESET,ICONVERGE,LINPOS,LINNEG,ILGRET,IEMIN
       DIMENSION GOLD(MAXOPT),GBEST(MAXOPT),H(MAXOPT)
       DIMENSION YOLD(MAXOPT),YBEST(MAXOPT),U(MAXOPT)
       DIMENSION GAMMA(MAX2),FUNCT(MAX2),DERIV(MAX2)
       CHARACTER*10 STRDUM
       SAVE
       DATA FIRST /.TRUE./
       DATA EPS   /1.0D-6/
       DATA XTOL  /1.0D-5/
       DATA FACCU /1.0D-4/
       DATA EACCU /5.0D-6/
C
C CHECK IF NPAR IS IN BOUNDS
C 
       IF ((NPAR.GT.MAXOPT).OR.(NPAR.LT.1)) THEN
        PRINT *,'cgrad: (NPAR > MAXOPT) or (NPAR < 1)'
        STOP
       END IF
C
C FILE OPERATIONS ACCORDING TO MODE
C
       IF (MODE.GT.0) THEN
        INQUIRE(75,OPENED=LOPF)
        IF (LOPF) CLOSE(75)
        INQUIRE(77,OPENED=LOPF)
        IF (LOPF) CLOSE(77)
        IF (MODE.EQ.2) THEN
         OPEN(75,FILE='cgrad',FORM='UNFORMATTED',STATUS='UNKNOWN')
         CLOSE(75,STATUS='DELETE')
        END IF
        RETURN
       END IF
C
C IF (FIRST), OPEN FILES AND CHECK IF THEY ARE VALID
C
       IRESET= .FALSE.
       ILGRET= .TRUE.
       IF (FIRST) THEN
        FIRST= .FALSE.
        OPEN(75,FILE='cgrad',FORM='UNFORMATTED',STATUS='UNKNOWN')
        OPEN(77,FILE='cgrlog',FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(75)
        REWIND(77)
C
C CHECK IF FILE CGRAD IS OKAY. IF NOT, START NEW RELAXATION
C
        ATOL=0.5*GTOL
        ISTEP=0
        ILINE=0
        NPARRD=0
        DELTA=0.0d0
        READ(75,END=10) ISTEP,ILINE,NPARRD,DELTA
        IF (NPARRD.NE.NPAR) GOTO 10
        GOTO 20 
   10   IRESET= .TRUE.
        ISTEP=0
        ILINE=0
   20   CONTINUE
C
C READ IN CGRLOG DATA. IF (IRESET), REMOVE THE OLD FILE   
C
        DMAG= 1.61D0 
        SMALL=0.10D0
        READ(77,*,END=30) SMALL,DMAG
   30   CONTINUE 
        IF (IRESET) THEN
         CLOSE(77,STATUS='DELETE')
         OPEN(77,FILE='cgrlog',FORM='FORMATTED',STATUS='UNKNOWN') 
         REWIND(77)
         WRITE(77,*) SMALL,DMAG
         DELTA=SMALL
        ELSE
   34     READ(77,'(A10)',END=36) STRDUM
          GOTO 34
   36    BACKSPACE(77)
        END IF
C
C MAKE ENTRY TO CGRLOG, IF NECESSARY
C
        IF (GTOL.LT.FACCU) CALL LOGCGR(77,1,ILGRET,0,0,0)
       END IF
C
C SET UP VALUES FOR CONJUGATE GRADIENT MINIMIZATION
C
       REWIND(75)
       IF (.NOT.IRESET) THEN
        READ(75) ISTEP,ILINE,NPARRD,DELTA
        IF (ILINE.GE.1) GO TO 100
        READ(75)(YOLD(I),I=1,NPAR)
        READ(75)(YBEST(I),I=1,NPAR)
        READ(75)(GOLD(I),I=1,NPAR)
        READ(75)(GBEST(I),I=1,NPAR)
        READ(75)(H(I),I=1,NPAR)
       END IF 
   40  CONTINUE
C
C FIND CONJUGATE DIRECTION OF TRAVEL (VECTOR H)
C USING POLAK-RIBIERE FORMULA
C
       ISTEP=ISTEP+1
       ILINE=0
       IF(ISTEP.GT.NPAR) THEN
        ISTEP=1
        DELTA=0.2*DELTA
       END IF
       IF(ISTEP.EQ.1)THEN
C
C START WITH NEGATIVE GRADIENT
C
        DO IPAR=1,NPAR
         H(IPAR)= -G(IPAR)
        END DO  
       ELSE
C
C CONSTRUCTION OF NEW H USING POLAK-RIBIERE
C
        GAM=0.0d0
        GDIV=0.0d0
        DO IPAR=1,NPAR
         GAM=GAM+(G(IPAR)-GOLD(IPAR))*G(IPAR)
         GDIV=GDIV+GOLD(IPAR)**2
        END DO
        GAM=GAM/GDIV 
        DO IPAR=1,NPAR
         H(IPAR)= -G(IPAR)+GAM*H(IPAR)
        END DO
       END IF
       HNRM=0.0d0
       DO IPAR=1,NPAR
        HNRM=HNRM+H(IPAR)**2
       END DO
       HNRM=1.0/SQRT(HNRM)
       DO IPAR=1,NPAR
        U(IPAR)=H(IPAR)*HNRM
       END DO  
C
C SETUP FOR INITIAL BEST VALUES
C
       EBEST=1.0D30
       FBEST=1.0D30
       DO IPAR=1,NPAR
        YBEST(IPAR)=Y(IPAR)
        GBEST(IPAR)=G(IPAR)
       END DO  
       GAMMA(1)= 0.0d0
       FUNCT(1)= 0.0d0
       DERIV(1)= 0.0d0
       REWIND(75)
       WRITE(75)ISTEP,ILINE,NPAR,DELTA
       WRITE(75)(Y(I),I=1,NPAR)
       WRITE(75)(YBEST(I),I=1,NPAR)
       WRITE(75)(G(I),I=1,NPAR)
       WRITE(75)(GBEST(I),I=1,NPAR)
       WRITE(75)(H(I),I=1,NPAR)
       WRITE(75)(U(I),I=1,NPAR)
       WRITE(75)(GAMMA(I),I=1,ILINE+1)
       WRITE(75)(FUNCT(I),I=1,ILINE)
       WRITE(75)(DERIV(I),I=1,ILINE)
       WRITE(75)EBEST,FBEST
C
C BEGIN/CONTINUE LINE MINIMIZTION IN DIRECTION OF U:
C
  100  CONTINUE
       REWIND(75)
       READ(75)ISTEP,ILINE,NPARRD,DELTA
       READ(75)(YOLD(I),I=1,NPAR)
       READ(75)(YBEST(I),I=1,NPAR)
       READ(75)(GOLD(I),I=1,NPAR)
       READ(75)(GBEST(I),I=1,NPAR)
       READ(75)(H(I),I=1,NPAR)
       READ(75)(U(I),I=1,NPAR)
       READ(75)(GAMMA(I),I=1,ILINE+1)
       READ(75)(FUNCT(I),I=1,ILINE)
       READ(75)(DERIV(I),I=1,ILINE)
       READ(75)EBEST,FBEST
C
C SETUP FOR LINE MINIMIZATION
C
       IEMIN=.TRUE.
       EDIFF=EACCU*ABS(EBEST)
       ILINE=ILINE+1
       FUNCT(ILINE)=ENERGY
       DERIV(ILINE)=0.0d0
       DO IPAR=1,NPAR
        DERIV(ILINE)=DERIV(ILINE)+U(IPAR)*G(IPAR)
       END DO  
C
C CHANGE BEST VALUES IF NECESSARY
C
       IF (ABS(EBEST-ENERGY).LT.EDIFF) THEN
        IF (ABS(DERIV(ILINE)).LT.FBEST) THEN
         EBEST=ENERGY
         FBEST=ABS(DERIV(ILINE))
         DO IPAR=1,NPAR
          YBEST(IPAR)=Y(IPAR)
          GBEST(IPAR)=G(IPAR)
         END DO
        END IF
       ELSE IF(ENERGY.LE.EBEST)THEN
        EBEST=ENERGY
        FBEST=ABS(DERIV(ILINE))
        DO IPAR=1,NPAR
         YBEST(IPAR)=Y(IPAR)
         GBEST(IPAR)=G(IPAR)
        END DO 
       END IF
C
C CALCULATE NEXT VALUE OF GAMMA:
C
       IMODE=0
       IF(ILINE.EQ.1)THEN
        IMODE=1
        GAMMA(1)=0.0d0
        IF(DERIV(1).GT.0.0)THEN
         GAMMA(2)= -DELTA
        ELSE
         GAMMA(2)=  DELTA
        END IF
       ELSE 
C
C CHECK WHETHER BRACKETS ALREADY AVAILABLE
C 
        LINPOS=.FALSE.
        LINNEG=.FALSE.
        DO JLINE=1,ILINE
         IF (DERIV(JLINE).GT.0.0) THEN
          LINPOS=.TRUE.
         ELSE
          LINNEG=.TRUE.
         END IF
        END DO
C
C IF NO BRACKETING, MAGNIFY INTERVAL
C
        IF (.NOT.(LINPOS.AND.LINNEG)) THEN
         IMODE=1
         GAMMA(ILINE+1)=DMAG*GAMMA(ILINE)+GAMMA(2)
        ELSE
C
C DO BETTER UPDATE, FOR (ILINE.EQ.2) TRY LINEAR INTERPOLATION 
C
         IF(ILINE.EQ.2)THEN
          GAMMA(3)=0.5*GAMMA(2)
          IF (ABS(DERIV(1)-DERIV(2)).GT.EPS) THEN
           IMODE=3
           GAMMA(3)=GAMMA(1)*DERIV(2)-GAMMA(2)*DERIV(1)
           GAMMA(3)=GAMMA(3)/(DERIV(2)-DERIV(1))
          ELSE
           IMODE=4
           GAMMA(3)=0.5*GAMMA(2)
          END IF
         ELSE 
C
C NOW ILINE IS GREATER THAN OR EQUAL TO 3
C CHECK IF LAST ITERATION LEAD TO GOOD RESULT
C
          ICONVERGE=.FALSE.
          IF (ABS(DERIV(ILINE)).LT.ATOL) ICONVERGE=.TRUE.
C
C CHECK FOR APPARENT PILE UP OF POINTS NEAR THE MINIMUM:
C
          IF (ILINE.GT.3) THEN
           DISTANCE=0.0d0
           DO KLINE=1,2
            DO JLINE=KLINE+1,3
             DISTANCE=DISTANCE+ABS(GAMMA(KLINE)-GAMMA(JLINE))
            END DO
           END DO
           IF(DISTANCE.LT.2.0*XTOL)ICONVERGE=.TRUE.
          END IF
C
C IF LINEMIN CONVERGED, CALCULATE NEW DIRECTION 
C
          IF(ICONVERGE) THEN
           GOTO 40 
          ELSE
C
C FIRST TRY QUADRATIC INTERPOLATION OF DERIVATIVES, 
C THEN LINEAR AND IF NOTHING HELPS BISECTION        
C
C SORT GAMMA, FUNCT AND DERIV WITH RESPECT TO FUNCT
C (BEST GAMMA -> GAMMA(0))
C          
           DO JLINE=1,ILINE
            DO KLINE=JLINE+1,ILINE
             IF (FUNCT(KLINE).LT.FUNCT(JLINE)) THEN
              CALL SWAP(GAMMA(KLINE),GAMMA(JLINE))
              CALL SWAP(FUNCT(KLINE),FUNCT(JLINE))
              CALL SWAP(DERIV(KLINE),DERIV(JLINE))
             END IF
            END DO
           END DO
C
C CHECK FOR PILE OF POINTS HAVING FUNCTION VALUES THAT
C ARE IN THE RANGE OF EDIFF. IF MORE THAN 3 POINTS ARE IN
C THE PILE, SWITCH TO DERIVATIVE CHECK
C
           DO JLINE=2,ILINE
            IF (ABS(FUNCT(JLINE)-FUNCT(1)).GT.EDIFF) THEN
             NPILE=JLINE-1
             GOTO 55
            END IF
           END DO
           NPILE=ILINE
   55      CONTINUE
           IF (NPILE.GT.3) THEN
            IEMIN=.FALSE.
            CALL LOGCGR(77,2,ILGRET,0,0,0)
            DO JLINE=1,NPILE-1
             DO KLINE=JLINE+1,NPILE
              IF (ABS(DERIV(KLINE)).LT.ABS(DERIV(JLINE))) THEN
               CALL SWAP(GAMMA(KLINE),GAMMA(JLINE))
               CALL SWAP(FUNCT(KLINE),FUNCT(JLINE))
               CALL SWAP(DERIV(KLINE),DERIV(JLINE))
              END IF
             END DO
            END DO
           END IF
C
C GET INDICES OF THE THREE INTERESTING POINTS AND SORT THEM
C          
           IF (DERIV(1)*DERIV(2).LT.0.0) THEN
            IND=3
           ELSE
            DO JLINE=3,ILINE
             IF (DERIV(JLINE)*DERIV(1).LT.0.0) THEN
              IND=JLINE
              GOTO 105
             END IF
            END DO
  105       CONTINUE
            CALL SWAP(GAMMA(3),GAMMA(IND))
            CALL SWAP(FUNCT(3),FUNCT(IND))
            CALL SWAP(DERIV(3),DERIV(IND))
           END IF
           DO JLINE=1,2
            DO KLINE=JLINE+1,3
             IF (GAMMA(KLINE).LT.GAMMA(JLINE)) THEN
              CALL SWAP(GAMMA(KLINE),GAMMA(JLINE))
              CALL SWAP(FUNCT(KLINE),FUNCT(JLINE))
              CALL SWAP(DERIV(KLINE),DERIV(JLINE))
             END IF
            END DO
           END DO
C
C GET BRACKETING TRIPLE A,B,C
C
           AG=GAMMA(1)
           BG=GAMMA(2)
           CG=GAMMA(3)
           A1=DERIV(1)
           B1=DERIV(2)
           C1=DERIV(3)
           LINNEG=.FALSE.
           IF (B1.LT.0.0) LINNEG=.TRUE.
C
C CHECK FOR UNUSUAL BEHAVIOR OF FUNCTION
C IF ((A1 > 0) && (C1 < 0)) -> MAXIMUM IN (A,C) 
C                           -> ENLARGE INTERVAL
C
           IF ((A1.GT.0.0).AND.(C1.LT.0.0)) THEN
            IMODE=1
            CALL LOGCGR(77,3,ILGRET,0,0,0)
            GNEW=CG+DELTA
            IF (IEMIN) THEN
             IF (FUNCT(1).LT.FUNCT(3)) GNEW=AG-DELTA
            ELSE
             IF (ABS(A1).GT.ABS(C1)) GNEW=AG-DELTA
            END IF
            GOTO 200
C
C IF ((A1 > 0) && (C1 > 0)) -> MAXIMUM IN (A,B), MINIMUM IN (B,C) 
C                           -> LINEAR INTERPOLATION IN (B,C)
C
           ELSE IF ((A1.GT.0.0).AND.(C1.GT.0.0)) THEN
            CALL LOGCGR(77,4,ILGRET,0,0,0)
            GOTO 110
C
C IF ((A1 < 0) && (C1 < 0)) -> MINIMUM IN (A,B), MAXIMUM IN (B,C) 
C                           -> LINEAR INTERPOLATION IN (A,B)
           ELSE IF ((A1.LT.0.0).AND.(C1.LT.0.0)) THEN
            CALL LOGCGR(77,5,ILGRET,0,0,0)
            GOTO 120
           END IF
C
C THANK GOD, EVERYTHING IS OKAY
C TRY QUADRATIC INTERPOLATION
C
           G2=BG-AG
           G3=CG-AG
           D2=B1-A1
           D3=C1-A1
           DN=(G3-G2)*G3*G2
           IF (ABS(DN).LE.1D-20) GOTO 110
           DN=1.0/DN
           APOL=(D2*G3*G3-D3*G2*G2)*DN
           BPOL=(D3*G2-D2*G3)*DN
           IF (ABS(BPOL).LE.1D-10) GOTO 110
           BPOLR=1.0/BPOL
           DHLP= -0.5*APOL*BPOLR
           DRT= DHLP*DHLP-A1*BPOLR
           IF (DRT.LT.1D-20) GOTO 110
           DRT=SQRT(DRT)
           GNEW=DHLP+DRT
           IF ((2*BPOL*GNEW+APOL).LT.0.0) GNEW=DHLP-DRT
           GNEW=GNEW+AG
           IF ((GNEW.LT.AG).OR.(GNEW.GT.CG)) GOTO 110 
           IMODE=2
           GOTO 200
C
C SET UP VALUES FOR LINEAR INTERPOLATION IN INTERVAL (A,B)
C
  110      IF (B1.LT.0.0) THEN
            AG=BG
            BG=CG
            A1=B1
            B1=C1
           END IF
C
C TRY LINEAR INTERPOLATION 
C
  120      IF ((B1-A1).LE.1D-10) GOTO 140
           GNEW=(B1*AG-A1*BG)/(B1-A1)
           IF ((GNEW.LT.AG).OR.(GNEW.GT.CG)) GOTO 140 
           IMODE=3
           GOTO 200     
C
C BISECTION
C
  140      IMODE=4
           GNEW=0.5*(AG+BG)   
  200      GAMMA(ILINE+1)=GNEW
          END IF
         END IF
        END IF
       END IF
C
C WRITE RESULTS IN LOGFILE
C
       ILGRET=.FALSE.
       CALL LOGCGR(77,0,ILGRET,IMODE,ILINE,ISTEP)
C
C IF LINMAX EXCEEDED, RETURN BEST VALUE OF Y AND START ANEW
C
       IF(ILINE.GE.MAXLINE)THEN
        DO I=1,NPAR
         Y(I)=YBEST(I)
         G(I)=GBEST(I)
        END DO
        ENERGY=EBEST
        GOTO 40
       END IF
C
C UPDATE Y
C
       DO I=1,NPAR
        Y(I)=YOLD(I)+GAMMA(ILINE+1)*U(I)
       END DO
C
C RESTORE HISTORY FOR LINE MINIMIZATION
C
       REWIND(75)
       WRITE(75)ISTEP,ILINE,NPARRD,DELTA
       WRITE(75)(YOLD(I),I=1,NPAR)
       WRITE(75)(YBEST(I),I=1,NPAR)
       WRITE(75)(GOLD(I),I=1,NPAR)
       WRITE(75)(GBEST(I),I=1,NPAR)
       WRITE(75)(H(I),I=1,NPAR)
       WRITE(75)(U(I),I=1,NPAR)
       WRITE(75)(GAMMA(I),I=1,ILINE+1)
       WRITE(75)(FUNCT(I),I=1,ILINE)
       WRITE(75)(DERIV(I),I=1,ILINE)
       WRITE(75)EBEST,FBEST
       RETURN
       END
C
C**********************************************************************

        SUBROUTINE SWAP(X,Y)
        IMPLICIT NONE
        REAL*8 X,Y,TEMP
        TEMP=X
        X=Y
        Y=TEMP
        RETURN
        END

C**********************************************************************

        SUBROUTINE ISWAP(I,J)
        IMPLICIT NONE
        INTEGER I,J,TEMP
        TEMP=I
        I=J
        J=TEMP
        RETURN
        END
C
