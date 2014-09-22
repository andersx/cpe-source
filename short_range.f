c==============================================================================
c      evaluate short range expression i.e. sumR (gamma - 1/R)
c      
c      INPUT:
c      REAL*8    rh(3)       vector between rmu and rnu
c      REAL*8    basis(3,3)  basis of cell(unusual!!!:one line is a basis vector)
c      REAL*8    umu         hubbard parameter of orbital mu
c      REAL*8    unu         hubbard parameter of orbital nu
c      REAL*8    udermu      hubbard derivative of orbital mu
c      logical   xhgammahlp  .true. if h=exp(-((ui+uj)/2^zeta*r^2)) for gamma^h
c      REAL*8    zeta        parameter for gamma^h (see Gaus JCTC 2011)
c      REAL*8    tol         convergence tolerance (contribution of last shell)  
C      !!NOTE THAT umu,unu,udermu in this implementation is not orbital but atom specific!!
c
c      OUTPUT:
c      REAL*8    value       value of the short range sum (gamma)
c      REAL*8    valueder    value of the short range sum (Gamma=dgamma/dq)
c==============================================================================

       subroutine SHORTRANGE(rh,basis,umu,unu,udermu,xhgammahlp,
     &                       zeta,tol,value,valueder)
       IMPLICIT NONE

       REAL*8  rh(3),umu,unu,udermu,basis(3,3),zeta,tol,value,valueder
       INTEGER i,j,k,nreal,nmax,nmin 
       REAL*8  rvec(3),R(3),result,resultder,lastshell,tmp,norm,lastder
       REAL*8  tolder,gval,gder
       LOGICAL xhgammahlp
       character*1  sccmode
       common /scchelp/ sccmode

c
c      rh = rmu - rnu
c
       result = 0.0d0
       resultder = 0.0d0
       nmax = 50
       nmin = 3
       nreal = 0
       lastshell = tol+1.0d-8
c      /* sum over R until tolerance is reached */
c MG_UW1211: insert criteria for 3rd-order convergence
c            since in gammaall gder is always calculated as for DFTB3 
c            we just rise the tolder if DFTB2 is invoked 
       if (sccmode=="3") then
         lastder = tol+1.0d-8
         tolder  = tol
       else
         lastder = 0.0d0
         tolder  = tol*1.0d20 
       endif
!      /* sum over R until tolerance is reached */
       DO WHILE ((nreal .le. nmax) .and. ((dabs(lastshell) .gt. tol)
     &      .or.(dabs(lastder).gt.tolder)  .or. (nreal .le. nmin)  )) 
       lastshell = 0.0d0
       lastder   = 0.0d0
        DO i = -nreal,nreal
       	 DO j = -nreal,nreal
       	  DO k = -nreal,nreal
c            /*only R belonging to outer shells are new ones */
             IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j))  .or. 
     &	         (nreal.eq. abs(k)) ) THEN
                  
                  R(1)=i*basis(1,1)+j*basis(2,1)+k*basis(3,1) 
                  R(2)=i*basis(1,2)+j*basis(2,2)+k*basis(3,2) 
                  R(3)=i*basis(1,3)+j*basis(2,3)+k*basis(3,3) 

                  rvec(1) = rh(1) - R(1)  
                  rvec(2) = rh(2) - R(2)   
                  rvec(3) = rh(3) - R(3)  
              
                  norm   = dsqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

c              get value for gamma or gamma^h and Gamma
               call gammaall(norm,umu,unu,udermu,xhgammahlp,
     &                      zeta,gval,gder)

c              subtract long range part 1/R and multiply by Z(nu)
                           
               IF (norm .lt. 1.0d-6) THEN 
                  tmp =  gval
               ELSE
                  tmp =  ( gval  - 1.0d0/norm )
               ENDIF

               result = result + tmp
               resultder = resultder + gder
               ! for physical reasonable value gder converges much faster than gval, 
               ! however gval-1/norm converges faster than gder
               ! that is why additional criterium (lastder) is applied.
               lastshell = lastshell + tmp
               lastder   = lastder + gder
              
             END IF
            END DO
           END DO
          END DO
          nreal = nreal + 1
          END DO

        IF((dabs(lastshell) .gt. tol).or.(dabs(lastder).gt.tolder)) THEN
           STOP "tolerance in subroutine short not reached."
        END IF
        value = result 
        valueder = resultder
          

       END




c=============================================================================
c evaluate derivative of short range expression: sumR (d gamma/dR - (-1/R^2))
c      
c      INPUT:
c      REAL*8    rh(3)       vector between rmu and rnu
c      REAL*8    umu         hubbard parameter of orbital mu
c      REAL*8    unu         hubbard parameter of orbital nu
c      REAL*8    udermu      hubbard derivative of orbital mu
c      logical   xhgammahlp  .true. if h=exp(-((ui+uj)/2^zeta*r^2)) for gamma^h
c      REAL*8    zeta        parameter for gamma^h (see Gaus JCTC 2011)
c      REAL*8    basis(3,3)  basis of cell(unusual!!!:one line is a basis vector)
c      REAL*8    tol         convergence tolerance (contribution of last shell)  
c    
c      OUTPUT:
c      REAL*8    deriv(3)    derivative of the short range sum (gamma contr)
c      REAL*8    deriv3(3)   derivative of the short range sum (Gamma contr)
c==============================================================================

       subroutine SHORTRANGE1(rh,basis,umu,unu,udermu,xhgammahlp,zeta,
     &                        tol,deriv,deriv3)
       IMPLICIT NONE

       REAL*8 rh(3),umu,unu,udermu,basis(3,3),tol,deriv(3),deriv3(3)
       REAL*8 zeta
       INTEGER i,j,k,nreal,nmax,nmin
       REAL*8 rvec(3),R(3),lastshell,tmp,norm,lastder,tolder
       REAL*8 gdrv,gdrv3
       LOGICAL xhgammahlp
       character*1  sccmode
       common /scchelp/ sccmode

c
c      rh = rmu - rnu
c
       deriv(1) = 0.0
       deriv(2) = 0.0
       deriv(3) = 0.0
       deriv3(1) = 0.0d0
       deriv3(2) = 0.0d0
       deriv3(3) = 0.0d0
       nmax = 100
       nmin = 3
       nreal = 0
       
       lastshell = tol+1.0d-8
c MG_UW1211: insert criteria for 3rd-order convergence
c            since in gammaall gder is always calculated as for DFTB3 
c            we just rise the tolder if DFTB2 is invoked 
       if (sccmode=="3") then
         lastder = tol+1.0d-8
         tolder  = tol
       else
         lastder = 0.0d0
         tolder  = tol*1.0d20 
       endif
c      /* sum over R until tolerance is reached */
       DO WHILE ((nreal .le. nmax) .and. ((dabs(lastshell) .gt. tol)
     &      .or.(dabs(lastder).gt.tolder)  .or. (nreal .le. nmin)  )) 
       lastshell = 0.0d0
       lastder   = 0.0d0
        DO i = -nreal,nreal
         DO j = -nreal,nreal
          DO k = -nreal,nreal
c      /*only R belonging to outer shells are new ones */
       IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j))  .or. 
     &  (nreal.	eq. abs(k)) ) THEN

         R(1)=i*basis(1,1)+j*basis(2,1)+k*basis(3,1) 
         R(2)=i*basis(1,2)+j*basis(2,2)+k*basis(3,2) 
         R(3)=i*basis(1,3)+j*basis(2,3)+k*basis(3,3) 


         rvec(1) = rh(1) - R(1)
         rvec(2) = rh(2) - R(2)
         rvec(3) = rh(3) - R(3)

         norm   = dsqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

c       get derivative of gamma and Gamma

        call gammaall1(norm,umu,unu,udermu,xhgammahlp,zeta,gdrv,gdrv3)
                       
c       subtract long range -1/R^2
        IF (norm .lt. 1.0d-6) THEN
         tmp = gdrv
        ElSE
         tmp = gdrv + 1.0d0/(norm**2)
        ENDIF
                            
        deriv(1) = deriv(1) + tmp*rvec(1)/norm
        deriv(2) = deriv(2) + tmp*rvec(2)/norm
        deriv(3) = deriv(3) + tmp*rvec(3)/norm
        deriv3(1)= deriv3(1) + gdrv3*rvec(1)/norm
        deriv3(2)= deriv3(2) + gdrv3*rvec(2)/norm
        deriv3(3)= deriv3(3) + gdrv3*rvec(3)/norm

        lastshell = lastshell + tmp
        lastder   = lastder + gdrv3
        END IF
       END DO
       END DO
       END DO
       nreal = nreal + 1
       END DO
       
       IF((dabs(lastshell) .gt. tol).or.(dabs(lastder).gt.tolder)) THEN
        STOP "tolerance in subroutine short1 not reached."
       END IF

       END
