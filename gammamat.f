c======================================================================
c   Build symmetric matrix gammamat containing Ewald potentials + short
c         range terms for periodic systems
c   Build non-symmetric matrix gammader
c
c   INPUT:
c   integer nn            number of atoms
c   real*8  x(3,NNDIM)    position of atoms
c   integer izp(NNDIM)    map from atoms to atom types
c   real*8  uhubb(MAXTYP,3) Hubbard parameters (l-dep)
c   real*8  uhder(MAXTYP,3) Hubbard derivatives dU/dq (l-dep)
c   real*8  zeta          parameter for gamma^h (see Gaus JCTC 2011)
c   real*8  basis(3,3)    basis of cell
c   logical period        .true. if periodic boundary conditions
c   logical izpxh(MAXTYP) .true. for atom types which need extra term
c                         in gamma^h if switched on
c
c   OUTPUT:
c   real*8 gammamat(*,*,3,3) matrix containing Ewald potentials + short 
c                        range terms  (for molecules the simple 
c                        gamma/gamma^h for all atom-pairings (ldep!))
c   real*8 gammader(*,*,3,3) matrix containing Gamma=dgamma/dq for DFTB3 
c                        for all atom-pairings (ldep!)
c                        (for periodic systems: only short range terms)
c
c   Note that code is made efficient (but still easily readable) 
c   for DFTB3, but allows also running DFTB2, therefore gammader 
c   is calculated by default in this function but of course may 
c   be zeroed or controlled by a subroutine calling get_gammamat 
c   or using the OUTPUT.
c
c======================================================================

       subroutine get_gammamat(nn,x,izp,uhubb,uhder,zeta,basis,period,
     &                     izpxh,ldep,lmax,gammamat,gammader)
       implicit none
       include 'maxima.inc'
       ! i/o variables
       integer  nn,izp(NNDIM),lmax(MAXTYP)
       logical  period,izpxh(MAXTYP),ldep
       real*8   x(3,NNDIM),uhubb(MAXTYP,3),uhder(MAXTYP,3)
       real*8   zeta,basis(3,3)
       real*8   gammamat(NNDIM,NNDIM,3,3),gammader(NNDIM,NNDIM,3,3)
       ! internal variables
       external getalpha
       integer  i,j,li,lj
       logical  xhgammahlp
       real*8   recbasis(3,3),vol,alpha,getalpha,tol,r(3),phivalue
       real*8   gval,gder,norm

       if (period) then

         ! get reciprocal lattice vectors and cell volume
         call rezvol(basis,recbasis,vol)
         ! choose good convergence parameter alpha
         alpha = getalpha(basis)
         ! set rolerance for convergence
         tol   = 1.0d-8
         do i=1,nn
           do j=1,nn
             r(1)=x(1,i)-x(1,j)
             r(2)=x(2,i)-x(2,j)
             r(3)=x(3,i)-x(3,j)
             call phi(r,basis,recbasis,alpha,vol,tol,phivalue)
             if ((izpxh(izp(i))).or.(izpxh(izp(j)))) then
               xhgammahlp=.true.
             else
               xhgammahlp=.false.
             endif
             if (ldep) then
                 do li=1,lmax(izp(i))
                 do lj=1,lmax(izp(j))
                   call shortrange(r,basis,uhubb(izp(i),li),
     &                             uhubb(izp(j),lj),uhder(izp(i),li),
     &                             xhgammahlp,zeta,tol,gval,gder)
                   gammamat(i,j,li,lj) = phivalue + gval
                   gammader(i,j,li,lj) = gder
                 enddo
                 enddo
             else !no ldep
               call shortrange(r,basis,uhubb(izp(i),1),uhubb(izp(j),1),
     &                uhder(izp(i),1),xhgammahlp,zeta,tol,gval,gder)
               gammamat(i,j,1,1) = phivalue + gval
               gammader(i,j,1,1) = gder
             endif
           enddo
         enddo

       else ! (no period)

         do i=1,nn
           do j=1,nn
             r(1)=x(1,i)-x(1,j)
             r(2)=x(2,i)-x(2,j)
             r(3)=x(3,i)-x(3,j)
             norm = dsqrt(r(1)**2+r(2)**2+r(3)**2)
             if ((izpxh(izp(i))).or.(izpxh(izp(j)))) then
               xhgammahlp=.true.
             else
               xhgammahlp=.false.
             endif
             if (ldep) then
               do li=1,lmax(izp(i))
               do lj=1,lmax(izp(j))
                 call gammaall(norm,uhubb(izp(i),li),uhubb(izp(j),lj),
     &                uhder(izp(i),li),xhgammahlp,zeta,gval,gder)
                 gammamat(i,j,li,lj) = gval
                 gammader(i,j,li,lj) = gder
               enddo
               enddo
             else
               call gammaall(norm,uhubb(izp(i),1),uhubb(izp(j),1),
     &              uhder(izp(i),1),xhgammahlp,zeta,gval,gder)
               gammamat(i,j,1,1) = gval
               gammader(i,j,1,1) = gder
             endif
           enddo
         enddo
 
       endif

       end

