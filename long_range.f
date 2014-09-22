      subroutine initphi(basis)
      implicit none
      real*8 basis(3,3)
      real*8 recbasis(3,3)
      real*8 vol
      include 'phihelp.inc'
      integer i,j,k,l,ind1,ind2,ind3
      real*8 G(4)
      call REZVOL(basis,recbasis,vol)
      DO i=-nmax,nmax
       ind3=i+nmax+1
       DO j=-nmax,nmax
        ind2=j+nmax+1
        DO k=-nmax,nmax
         ind1=k+nmax+1
         G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
         G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
         G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)
         G(4)=G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
         do l=1,4
           gmat(l,ind1,ind2,ind3)=G(l)
         end do

         G(1)=-(i*basis(1,1)+j*basis(2,1)+k*basis(3,1))
         G(2)=-(i*basis(1,2)+j*basis(2,2)+k*basis(3,2))
         G(3)=-(i*basis(1,3)+j*basis(2,3)+k*basis(3,3))

         do l=1,3
           sumlat(l,ind1,ind2,ind3)=G(l)
         end do
         
        END DO
       END DO
      END DO
      return
      end
c  =====================================================================
c  evaluation of the potential phi 
c
c   phi(r) = 4*pi/Omega ( Summe{G neq 0} e^{-G^2/{4 alpha^2}}/{G^2} cos(G r)
c            +Summe{R, R neq r} (1-erf(alpha*|R-r|))/|R-r|
c            -pi/(Omega*alpha^2)
c
c
c   INPUT Parameter:
c   REAL*8 r(3)           position of evaluation of potential
c   REAL*8 basis(3,3)     basis of cell
c   REAL*8 recbasis(3,3)      basis of reciprocal cell
c   REAL*8 alpha          convergence parameter
c   REAL*8 vol            cell volume
c   REAL*8 tol            tolerance for convergence of "last shell"
c                       (can often be used as a criterion for global
c                        convergence) 
c   OUTPUT:
c   REAL*8 potential      value of Ewald potential
c  
c  ======================================================================

        subroutine phi(r,basis,recbasis,alpha,vol,tol,potential)
        IMPLICIT NONE
        external terfc
        REAL*8 terfc
        REAL*8 r(3), basis(3,3), recbasis(3,3), alpha, vol, tol
        REAL*8 potential
        REAL*8 reciprocal,rspace,cterm
        REAL*8 G(3),rh(3),help,norm,lastshell
        REAL*8 MPI
        INTEGER nrezi, nreal
        INTEGER i,j,k,l,ind1,ind2,ind3

        include 'phihelp.inc'

        MPI = 3.14159265358979323846d0

c       evaluate reciprocal space term ( sum over G <> 0) ...  
c       /* sum over G until tolerance is reached */
        nrezi = 1
        lastshell = tol+1.0d-8  
        reciprocal = 0.0d0
        DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or.
     &    (dabs(lastshell) .gt.  tol)))
          lastshell = 0.0d0  
          DO i=-nrezi,nrezi
           ind3=i+nmax+1
           DO j=-nrezi,nrezi
            ind2=j+nmax+1
            DO k=-nrezi,nrezi
c             /*only G belonging to outer shells are new ones */
              IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or.
     &        (nrezi. eq. abs(k)) ) THEN
              ind1=k+nmax+1
              do l=1,3
                g(l)=gmat(l,ind1,ind2,ind3)
              end do
              help = gmat(4,ind1,ind2,ind3)
              help = dexp(-help/(4.0d0*alpha*alpha))/help
              help = dcos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help
 
                  reciprocal = reciprocal + help
                  lastshell = lastshell + help/vol       
              END IF
             END DO
            END DO
          END DO
         nrezi = nrezi + 1
        END DO

c       stop if tolerance not reached
        IF ( dabs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi not reached in reciprocal space"     
        END IF

        reciprocal=(4.0d0*MPI*reciprocal)/vol


c       evaluate  real space term (sum over R)   
c       /* sum over R until tolerance is reached */
        rspace = 0.0d0
        nreal = 0
         lastshell = tol+1.0d-8  
        DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin) 
     &            .or. (dabs(lastshell) .gt.  tol)))
         lastshell = 0.0d0  
         DO i=-nreal,nreal
          ind3=i+nmax+1
          DO j=-nreal,nreal
          ind2=j+nmax+1
           DO k=-nreal,nreal
c            /*only R belonging to outer shells are new ones */
             IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or.
     &         (nreal. eq. abs(k)) ) THEN
                ind1=k+nmax+1
                rh(1)=r(1)+sumlat(1,ind1,ind2,ind3)
                rh(2)=r(2)+sumlat(2,ind1,ind2,ind3)
                rh(3)=r(3)+sumlat(3,ind1,ind2,ind3)
                norm=dsqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3))
                IF (norm .gt. 1.0d-20) THEN 
c               erfc=1-erf   
                  help   = terfc(alpha*norm)/norm
                  rspace = rspace + help
                  lastshell = lastshell + help
                ELSE 
                   lastshell = tol+1.0d-8  
              END IF
             END IF
            END DO
           END DO
          END DO
         nreal = nreal + 1
        END DO

c       stop if tolerance not reached
        IF ( dabs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi not reached in real space"     
        END IF


c       evaluate constant term pi/(Omega*alpha^2) 
        cterm = -MPI/(vol*alpha*alpha)

c       if r = 0 there is another constant to be added   
        IF ((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) .lt. 1.0d-20) THEN
            cterm = cterm -2.0d0*alpha/dsqrt(MPI)
        END IF

c       set the value of the potential
        potential = reciprocal + rspace + cterm

        END 


 
c  =====================================================================
c  evaluation of the derivative of the potential phi 
c
c   INPUT Parameter:
c   REAL*8 r(3)           position of evaluation of potential
c   REAL*8 basis(3,3)     basis of cell
c   REAL*8 recbasis(3,3)      basis of reciprocal cell
c   REAL*8 alpha          convergence parameter
c   REAL*8 vol            cell volume
c   REAL*8 tol            tolerance for convergence of "last shell"
c                       (can often be used as a criterion for global
c                        convergence) 
c   OUTPUT:
c   REAL*8 deriv(3)       derivative of ewlad potential
c  
c  ======================================================================

 
        subroutine phi1(r,basis,recbasis, alpha,vol,tol,deriv)

        IMPLICIT NONE
        external terfc
        REAL*8 terfc 
        REAL*8 r(3), basis(3,3), recbasis(3,3), alpha, vol, deriv(3)
        REAL*8 reciprocal(3),rspace(3),MPI 
        REAL*8 G(3),rh(3),norm,help,tol,lastshell 
        INTEGER i,j,k,l, nrezi, nreal,ind1,ind2,ind3

        include 'phihelp.inc'

        MPI =  3.14159265358979323846d0

        IF((r(1)*r(1) + r(2)*r(2) + r(3)*r(3)) .lt. 1.0d-20) THEN
          deriv(1) = 0.0d0
          deriv(2) = 0.0d0
          deriv(3) = 0.0d0
          return
        END IF
 
c       /* evaluate reciprocal space term (sum over G <> 0) ...  */
        nrezi = 1
        lastshell = tol+1.0d-8
        reciprocal(1) = 0.0d0
        reciprocal(2) = 0.0d0 
        reciprocal(3) = 0.0d0 
        DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or.
     &    (dabs(lastshell) .gt.  tol)))
          lastshell = 0.0d0
        DO i=-nrezi,nrezi
         ind3=i+nmax+1
         DO j=-nrezi,nrezi
          ind2=j+nmax+1
          DO k=-nrezi,nrezi
c             /*only G belonging to outer shells are new ones */
              IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or.
     &        (nrezi. eq. abs(k)) ) THEN
              ind1=k+nmax+1
              do l=1,3
                g(l)=gmat(l,ind1,ind2,ind3)
              end do
              help = gmat(4,ind1,ind2,ind3)
 
                  help=dexp(-help/(4.0d0*alpha*alpha))/help 
 
                  help=-dsin(G(1)*r(1)+G(2)*r(2)+G(3)*r(3))*help 
 
                  reciprocal(1)=help*G(1) + reciprocal(1)
                  reciprocal(2)=help*G(2) + reciprocal(2)
                  reciprocal(3)=help*G(3) + reciprocal(3)

                  lastshell = lastshell + help/vol
              ENDIF
          END DO
         END DO
        END DO
        nrezi = nrezi + 1
        END DO
 
c       stop if tolerance not reached
        IF ( dabs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi1 not reached in reciprocal space"
        END IF

        reciprocal(1)=(4.0d0*MPI*reciprocal(1))/vol 
        reciprocal(2)=(4.0d0*MPI*reciprocal(2))/vol 
        reciprocal(3)=(4.0d0*MPI*reciprocal(3))/vol 
 
 
c       /* evaluate  real space term (sum over R) */
c       /* sum over R until tolerance is reached */
        rspace(1) = 0.0  
        rspace(2) = 0.0  
        rspace(3) = 0.0 
        nreal = 0
        lastshell = tol+1.0d-8
        DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin)
     &            .or. (dabs(lastshell) .gt.  tol)))
        lastshell = 0.0
        DO i=-nreal,nreal
         ind3=i+nmax+1
         DO j=-nreal,nreal
          ind2=j+nmax+1
          DO k=-nreal,nreal
c            /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or.
     &         (nreal. eq. abs(k)) ) THEN
                ind1=k+nmax+1
                rh(1)=r(1)+sumlat(1,ind1,ind2,ind3)
                rh(2)=r(2)+sumlat(2,ind1,ind2,ind3)
                rh(3)=r(3)+sumlat(3,ind1,ind2,ind3)
 
            norm=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3)) 
 
            help = (-2.0d0/dsqrt(MPI)*dexp(-alpha*alpha*norm*norm)*
     &             alpha*norm - terfc(alpha*norm))/(norm*norm*norm) 
 
            rspace(1) = rh(1)*help + rspace(1)  
            rspace(2) = rh(2)*help + rspace(2)
            rspace(3) = rh(3)*help + rspace(3)
   
            lastshell = lastshell + help
            END IF
          END DO
         END DO
        END DO
       nreal = nreal + 1
      END DO

c       stop if tolerance not reached
        IF ( dabs(lastshell)  .gt. tol ) THEN
         STOP "tolerance in phi1 not reached in real space"
        END IF


c       /* add real and reciprocal parts */
        deriv(1) = rspace(1)  + reciprocal(1) 
        deriv(2) = rspace(2)  + reciprocal(2) 
        deriv(3) = rspace(3)  + reciprocal(3) 

        END 



c==============================================================================
c
c       evaluate the cross product of A x B
c
c==============================================================================

        subroutine CROSS( A, B, C) 
        IMPLICIT NONE
          REAL*8 A(3), B(3), C(3)

          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        END 


c       get reciprocal lattice vectors and volume of unit cell       

        subroutine REZVOL(basis,recbasis,vol) 
        IMPLICIT NONE
        external cross
        REAL*8  MPI
        REAL*8  basis(3,3), recbasis(3,3), vol
        REAL*8  hv1(3), hv2(3), hv3(3), hv4(3), fac
        INTEGER i
        MPI = 3.14159265358979323846d0


          DO i=1,3, 1
            hv1(i)=basis(1,i)
            hv2(i)=basis(2,i)
            hv3(i)=basis(3,i)
          END DO

          call CROSS(hv2,hv3,hv4)
          vol = dabs(basis(1,1)*hv4(1)+basis(1,2)*hv4(2)+
     &                basis(1,3)*hv4(3))
          fac = 2.0d0*MPI/vol

          recbasis(1,1)=hv4(1)*fac
          recbasis(1,2)=hv4(2)*fac
          recbasis(1,3)=hv4(3)*fac

          call CROSS(hv3,hv1,hv4)
          recbasis(2,1)=hv4(1)*fac 
          recbasis(2,2)=hv4(2)*fac 
          recbasis(2,3)=hv4(3)*fac

          call CROSS(hv1,hv2,hv4)
          recbasis(3,1)=hv4(1)*fac
          recbasis(3,2)=hv4(2)*fac
          recbasis(3,3)=hv4(3)*fac

        END 



c==============================================================================
c
c     Returns the (tabulated) complementary error function erfc(x) 
c     with fractional error everywhere less than 1.2 x 10^-7
c
c==============================================================================

        FUNCTION terfc(x)
        REAL*8 terfc,x
        REAL*8 t,z

        z = dabs(x)
        t = 1.0d0/(1.0d0+0.5d0*z)
        terfc = t*dexp(-z*z-1.26551223d0+t*(1.00002368d0+
     &   t*(0.37409196d0+t*(0.09678418d0+t*(-0.18628806d0+
     &   t*(0.27886807d0+t*(-1.13520398d0+t*(1.48851587d0+
     &   t*(-0.82215223d0+t*0.17087277d0)))))))))
        if (x .lt. 0.0d0) terfc=2.0d0-terfc

        RETURN
        END


c==============================================================================
c      get optimal alpha for Ewald potential phi
c
c      INPUT:
c      REAL*8 basis(3,3)     basis of lattice                     
c
c      RETURNS:
c      REAL*8                optimal alpha
c==============================================================================

       function getalpha(basis)
       IMPLICIT NONE
        REAL*8 basis(3,3)
        REAL*8 getalpha
        REAL*8 alpha, alphal, alphar
        INTEGER nopt
        REAL*8 recbasis(3,3), vol, tol
        REAL*8 G, R, help1, help2, help3
        EXTERNAL diffrecreal
        REAL*8 diffrecreal
       
        tol  = 1.0d-5

c       get reciprocal lattice vectors and cell volume

        CALL REZVOL(basis,recbasis,vol)

c       get sqnorm of smallest vector in reciprocal space 
        help1 = recbasis(1,1)**2+recbasis(1,2)**2+recbasis(1,3)**2
        help2 = recbasis(2,1)**2+recbasis(2,2)**2+recbasis(2,3)**2
        help3 = recbasis(3,1)**2+recbasis(3,2)**2+recbasis(3,3)**2
        G    = dsqrt(min(help1,help2,help3))
        
c       get norm of smallest vector in real space 
        help1 = basis(1,1)**2 + basis(1,2)**2 + basis(1,3)**2
        help2 = basis(2,1)**2 + basis(2,2)**2 + basis(2,3)**2
        help3 = basis(3,1)**2 + basis(3,2)**2 + basis(3,3)**2
        R    = dsqrt(min(help1,help2,help3))


c       optimise alpha

c       if (reciprocalspace decline - realspace decline) < 0 convergence
c       in real space too slow: increase alpha
c       set starting alphal
        alpha = 1.0d-5
         DO WHILE( diffrecreal(alpha,G,R,vol) .lt. tol )
          alphal = alpha
          alpha = alpha*2.0d0
         END DO

c       if (reciprocalspace decline - realspace decline) > 0 convergence
c       in reciprocal space too slow: decrease alpha
c       set starting alphar
        alpha = 1.0d+5
         DO WHILE( diffrecreal(alpha,G,R,vol) .gt. tol )
          alphar = alpha
          alpha = alpha / 2.0d0
         END DO

c       now find best value by refining the interval [alphal,alphar]
         alpha = (alphal + alphar)/2.0d0
         nopt  = 0
        DO WHILE ( dabs(diffrecreal(alpha,G,R,vol)) .gt. tol .and. 
     &             nopt .le. 20)
         IF ( diffrecreal(alpha,G,R,vol) .lt. tol ) THEN
          alphal = alpha
         END IF

         IF ( diffrecreal(alpha,G,R,vol) .gt. tol ) THEN
          alphar = alpha
         END IF
         alpha = (alphal + alphar)/2.0
         nopt  = nopt + 1
        END DO

        IF (nopt .gt. 20 ) THEN 
        alpha=dexp(-0.310104d0*dlog(vol)+0.786382d0)/2.0d0
         PRINT*, "WARNING: NO OPTIMISED ALPHA FOUND: "
         PRINT*, "STANDARD ALPHA USED. ALPHA SET TO", alpha 
        END IF
          
        getalpha = alpha

        END


c==============================================================================
c      
c      This subroutine returns the norm of the largest reciprocal space
c      vector and the norm of the largest real space vector needed for 
c      a converged Ewald summation
c
c      INPUT:
c      REAL*8  alpha         chosen convergence parameter
c      REAL*8  tol           disired tolerance of last contribution
c      REAL*8  vol           cell volume
c
c      OUTPUT:
c      REAL*8  Gmax          norm of largest vector in reciprocal space
c      REAL*8  Rmax          norm of largest vector in real space 
c==============================================================================


       subroutine getGRmax(alpha,tol,vol,Gmax,Rmax)
       IMPLICIT NONE
        REAL*8 alpha, Gmax, Rmax, tol, vol
        INTEGER nopt
        REAL*8 G, R, Gl, Gr, Rl, Rr
        EXTERNAL Gspace 
        REAL*8 Gspace 
        EXTERNAL Rspace 
        REAL*8 Rspace 
        EXTERNAL REZVOL
        


c       set starting Gl and Rl
        G = 1.0d-5
        R = 1.0d-5
         DO WHILE( Gspace(G,alpha,vol) .gt. tol + 1.0d-10)
          Gl = G
          G    = G*2.0d0
         END DO

         DO WHILE( Rspace(R,alpha) .gt. tol + 1.0d-10)
          Rl = R
          R    = R*2.0d0
         END DO

c       set starting Gr and Rr
        G = 1.0d+5
        R = 1.0d+5
         DO WHILE( Gspace(G,alpha,vol) .lt. tol - 1.0d-10)
          Gr = G
          G    = G/2.0d0
         END DO

         DO WHILE( Rspace(R,alpha) .lt. tol - 1.0d-10)
          Rr = R
          R    = R/2.0d0
         END DO



c       now find best value of G by refining the interval [Gl,Gr]
          G = (Gl + Gr)/2.0d0
         nopt  = 0
        DO WHILE ( Gspace(G,alpha,vol) .gt. tol .and. nopt .le. 20)
         IF ( Gspace(G,alpha,vol) .lt. (tol - 1.0d-10) ) THEN
          Gr = G
         END IF

         IF ( Gspace(alpha,G,vol) .gt. (tol + 1.0d-10) ) THEN
          Gl = G
         END IF

         G = (Gl + Gr)/2.0d0
         nopt  = nopt + 1
        END DO

        IF (nopt .ge. 20) THEN
          STOP "Gmax COULD NOT BE DETERMINED IN getGRmax"
        END IF

c       now find best value of R by refining the interval [Rl,Rr]
          R = (Rl + Rr)/2.0d0
         nopt  = 0
        DO WHILE ( Rspace(R,alpha) .gt. tol .and. nopt .le. 20)

         IF ( Rspace(R,alpha) .lt. (tol - 1.0d-10) ) THEN
          Rr = R
         END IF

         IF ( Rspace(alpha,R) .gt. (tol + 1.0d-10) ) THEN
          Rl = R
         END IF

         R = (Rl + Rr)/2.0d0
         nopt  = nopt + 1
        END DO

        IF (nopt .ge. 20) THEN
          STOP "Rmax COULD NOT BE DETERMINED IN getGRmax"
        END IF

        Gmax = G
        Rmax = R
        PRINT*,"GMAX:", GMAX, "REC:", Gspace(Gmax,alpha,vol)
        PRINT*,"RMAX:", RMAX, "REAL*8:", Rspace(Rmax,alpha)

        END



c==============================================================================
c      get differnce between decline in reciprocal and real space
c      this function is only used by function getalpha
c
c      INPUT:
c      REAL*8  alpha         convergence parameter
c      REAL*8  G             square of norm of smallest G
c      REAL*8  R             norm of smallest R
c      REAL*8  vol           cell volume
c
c      RETURNS:
c      REAL*8                difference between decline in reciprocal 
c                          space (rec(2G)-rec(3G)) and real space (real(2R) 
c                          - real(3R))
c==============================================================================

        function diffrecreal(alpha,G,R,vol)
        IMPLICIT NONE
        REAL*8 alpha, G, R, vol
        REAL*8 diffrecreal
        REAL*8 diffrec, diffreal
        EXTERNAL Gspace
        REAL*8 Gspace
        EXTERNAL Rspace
        REAL*8 Rspace
        REAL*8 MPI
        MPI = 3.14159265358979323846d0

c       make differences between decline at 2G and 3G / 2R and 3R
        diffrec = Gspace(2.0d0*G,alpha,vol) - Gspace(3.0d0*G,alpha,vol)
        diffreal= Rspace(2.0d0*R,alpha) - Rspace(3.0d0*R,alpha)

c       return difference between reciprocal and realspace decline
        diffrecreal        = diffrec - diffreal

        END


c==============================================================================
c       returns the "r independent" G space part of the Ewald sum
c      
c       INPUT:
c       REAL*8 G       norm of G
c       REAL*8 alpha   chosen convergence parameter
c       REAL*8 vol     cell volume
c
c       RETURNS:     
c       REAL*8         "r independent" G space part of the Ewald sum
c============================================================================== 

        function Gspace(G,alpha,vol)
 
        IMPLICIT NONE
        REAL*8 G, alpha, vol
        REAL*8 Gspace
        REAL*8 MPI
        MPI = 3.14159265358979323846d0

c       evaluate reciprocal space term at G
        Gspace = exp(-G**2/(4.0d0*alpha*alpha))/(G**2)
        Gspace = (4.0d0*MPI*Gspace)/vol

        END


c==============================================================================
c       returns the R space part of the Ewald sum
c      
c       INPUT:
c       REAL*8 R       norm of R
c       REAL*8 alpha   chosen convergence parameter
c
c       RETURNS:     
c       REAL*8         R space part of the Ewald sum
c============================================================================== 

        function Rspace(R,alpha)

        IMPLICIT NONE
        REAL*8 R, alpha
        REAL*8 Rspace
        EXTERNAL terfc
        REAL*8 terfc

c       evaluate real space term at R
        Rspace  = terfc(alpha*R)/R

        END
