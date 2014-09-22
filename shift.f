c=======================================================================
c      get qdiff and shift-vectors 
c      (necessary for Hubbard contribution to H matrix elements, total 
c       energy and forces)
c
c      INPUT:
c      integer     nn                     number of atoms (in one cell)
c      real*8      qmat(NNDIM)            mulliken charges 
c      real*8      qzero(MAXTYP,4)          neutral atomic charges (ldep+tot)
c      integer     izp(NNDIM)             map from atoms to atom types
c      real*8      gammamat(NNDIM,NNDIM,3,3)  gamma_{ij} for every atom 
c                                         pairing i and j (includes 
c                                         Ewald+shortrange if  
c                                         periodicity is switched on)
c      real*8      gammader(NNDIM,NNDIM,3,3)  Gamma_{ij)=dU_{ij}/dq_i for 
c                                         every atom pairing (ldep!)
c      character*1 sccmode                last term of DFT taylor series
c                                         which is included 
c                                         (e.g. 2=2nd order, 3=3rdorder)
c      logical     ldep                   true for l-dependent Hubbard
c      real*8      ql(3*NNDIM)            atom-shell-dependent charge (s,p,d)
c      logical     spin                   true for collinear spin-pola-
c                                         ization-formalism
c      real*8      qlup(3*NNDIM)          atom-shell-dependent (s,p,d)
c      real*8      qldown(3*NNDIM)        charges for up and down-electrons
c      real*8      wspin(MAXTYP,3,3)      spin-polarization constants
c      integer     indl(NNDIM+1)          index for l-dependent charges
c      integer     lmax(MAXTYP)           highest shell for element
c
c      OUTPUT: 
c      real*8  qdiff(NNDIM)   net charge of atoms
c      real*8  shift(NNDIM,3) shift(i) = \sum_{j}\Delta q_j \gamma_{ij}
c      real*8  shift3(NNDIM,3) shift3(i) = \sum_{j}\Delta q_i \Delta q_j
c                                          \Gamma_{ij}
c      real*8  shift3A(NNDIM,3) shift3A(i)= \sum_{j}\Delta q_j \Delta q_j
c                                          \Gamma_{ji}
c      real*8  spinshift(NNDIM,3)  ...(i,li)=sum_{l_i}
c
c  for l-dependence output is slightly different:
c      real*8  qdiff(NNDIM),qldiff(NNDIM) atomic and l-dep net charges
c      real*8  shift
c      real*8  shift3  
c      real*8  shift3A  
c      real*8  shift3B  

c======================================================================

       subroutine HAMILSHIFT(nn,qmat,qzero,izp,qdiff,gammamat,gammader,
     &               sccmode,shift,shift3,shift3A,shift3B,ldep,ql,
     &               spin,qlup,qldown,wspin,indl,lmax,spinshift)
       implicit none
       include 'maxima.inc'
       integer     nn,izp(NNDIM),indl(NNDIM+1),lmax(MAXTYP)
       character*1 sccmode
       logical     spin,ldep
       real*8      qmat(NNDIM),qzero(MAXTYP,4),qdiff(NNDIM)
       real*8      gammamat(NNDIM,NNDIM,3,3),gammader(NNDIM,NNDIM,3,3)
       real*8      shift(NNDIM,3),shift3(NNDIM,3),shift3A(NNDIM,3)
       real*8      qlup(3*NNDIM),qldown(3*NNDIM),spinshift(NNDIM,3)
       real*8      wspin(MAXTYP,3,3),ql(3*NNDIM),shift3B(NNDIM)
       ! internal variables
       integer     i,j,li,lj,indlj
       real*8      qldiff(NNDIM,3)

c  zero shift-vectors
       do i=1,nn
         shift3B(i) = 0.0d0
         do li=1,3
           shift(i,li)   = 0.0d0
           shift3(i,li)  = 0.0d0
           shift3A(i,li) = 0.0d0
           spinshift(i,li)=0.0d0
         enddo
       enddo

       if (ldep) then
c get qldiff 
         do i=1,nn
           qdiff(i) = qmat(i) - qzero(izp(i),4)
           do li=1,lmax(izp(i))
             qldiff(i,li) = ql(indl(i)+li) - qzero(izp(i),li)
           enddo
         enddo

c shift for 2nd order contributions
         do i=1,nn
         do j=1,nn
         do li=1,lmax(izp(i))
         do lj=1,lmax(izp(j))
           shift(i,li) = shift(i,li) + qldiff(j,lj)*gammamat(i,j,li,lj)
         enddo
         enddo
         enddo
         enddo

c shift for 3rd order contributions
         if (sccmode=="3") then
           do i=1,nn
           do li=1,lmax(izp(i))
           do j=1,nn
           do lj=1,lmax(izp(j))
            shift3(i,li) =shift3(i,li) +qldiff(j,lj)*gammader(i,j,li,lj)
            shift3A(i,li)=shift3A(i,li)+qldiff(j,lj)*gammader(j,i,lj,li)
     &                                 *qdiff(j)
            shift3B(i)   =shift3B(i)   +qldiff(j,lj)*gammader(i,j,li,lj)
     &                                 *qldiff(i,li)
           enddo
           enddo
            shift3(i,li)  = shift3(i,li) * qdiff(i)
           enddo
           enddo
         endif

c if not ldep
       else 
c get qdiff
       do i=1,nn
         qdiff(i) = qmat(i) - qzero(izp(i),4)
       enddo

c shift for 2nd order contributions
       do i=1,nn
         do j=1,nn
           shift(i,1) = shift(i,1) + qdiff(j)*gammamat(i,j,1,1)
         enddo
       enddo

c shift for 3rd order contributions
       if (sccmode=="3") then
         do i=1,nn
           do j=1,nn
          shift3(i,1) = shift3(i,1) +qdiff(j)*gammader(i,j,1,1)
          shift3A(i,1)= shift3A(i,1)+qdiff(j)*qdiff(j)*gammader(j,i,1,1)
           enddo
           shift3(i,1)  = shift3(i,1) * qdiff(i)
         enddo
       endif
c endif ldep
       endif

c shift for spin-polarized-formalism
       if (spin) then
         do i=1,nn
           do li=1,lmax(izp(i))
             do lj=1,lmax(izp(i))
               indlj=indl(i)+lj
               spinshift(i,li) = spinshift(i,li)
     &           +wspin(izp(i),li,lj)*(qlup(indlj)-qldown(indlj))
             enddo
           enddo
         enddo
       endif
       
       end


c=======================================================================
c get the gamma and Gamma contribution to the gradient
c -F_{kx}= 0.5d0*\Delta q_k\sum_{a!=k}\Delta q_a(dgamma_{ak}/dR_{kx}+
c          dgamma_{ka}/dR_{kx})+1/3 \Delta q_k\sum_{a!=k}\Delta q_a (
c          \Delta q_a dGamma_{ak}/dR_{kx}+\Delta q_k dGamma_{ak}/dR_{kx}
c          )
c For periodic systems: all expressions with gamma or Gamma also run
c over a sum of the cells (there is a cutoff)
c
c INPUT:
c integer  nn            number of atoms (in one cell)
c integer  nbeweg        number of movable atoms (in one cell)
c real*8   x(3,NNDIM)    coordinates
c real*8   izp(NNDIM)    map from atoms to atom types      
c real*8   uhubb(MAXTYP,3) Hubbard parameters
c real*8   uhder(MAXTYP,3) Hubbard derivatives
c real*8   basis         basis of cell
c logical  period        .true. if periodic boundary conditions
c logical  izpxh(MAXTYP) .true. for atom types which need extra term
c                         in gamma^h if switched on
c real*8   zeta          parameter for gamma^h (see Gaus JCTC 2011)
c real*8   qdiff(NNDIM)  atomic net charges (Mulliken)
c character*1 sccmode    last term of DFT taylor series which 
c                        is included (e.g. 2=2nd order, 3=3rdorder)
c
c OUTPUT:
c real*8   hgrad(3,NNDIM) gradient contribution 
c
c======================================================================

       subroutine GAMMAGRAD(nn,nbeweg,x,izp,uhubb,uhder,basis,period,
     &    izpxh,zeta,qdiff,sccmode,ldep,ql,qzero,lmax,indl,hgrad)
       implicit none
       include 'maxima.inc'
       integer     nn,nbeweg,izp(NNDIM),lmax(MAXTYP),indl(NNDIM+1)
       logical     period,izpxh(MAXTYP),ldep
       character*1 sccmode
       real*8      x(3,NNDIM),uhubb(MAXTYP,3),uhder(MAXTYP,3),basis(3,3)
       real*8      zeta,qdiff(NNDIM),ql(3*NNDIM),qzero(MAXTYP,4)
       real*8      hgrad(3,NNDIM)
       ! internal variables
       integer     ix,k,j,i,li,lk,lj
       logical     xhgammahlp
       external    getalpha
       real*8      recbasis(3,3),vol,tol,tmp(3),tmp3(3),r(3),getalpha
       real*8      alpha,long_deriv(3),short_deriv(3),norm
       real*8      dgdrkj,dgdr3kj,dgdrjk,dgdr3jk,hlp
       real*8      dgdr(NNDIM,NNDIM),dgdr3(NNDIM,NNDIM)
       real*8      SRdgdr(3,NNDIM,NNDIM),SRdgdr3(3,NNDIM,NNDIM)
       real*8      SRdgdrkj(3),SRdgdr3kj(3),SRdgdrjk(3),SRdgdr3jk(3)
       real*8      qldiff(NNDIM,3)

       do k=1,nn
         do ix=1,3
           hgrad(ix,k)=0.0d0
         enddo
       enddo

c first gradient for l-dependet Hubbards 
       ! not very efficient, but comprehensible formula hgrad=...
       if (ldep) then
c get qldiff
         do i=1,nn
           do li=1,lmax(izp(i))
             qldiff(i,li) = ql(indl(i)+li) - qzero(izp(i),li)
           enddo
         enddo
c gradient for periodic systems 
         if (period) then
           call rezvol(basis,recbasis,vol)
           alpha = getalpha(basis)
           tol = 1.0d-8
           do k=1,nbeweg
           do j=1,nn
           if (k.ne.j) then
             r(1)=x(1,k)-x(1,j)
             r(2)=x(2,k)-x(2,j)
             r(3)=x(3,k)-x(3,j)
             call phi1(r,basis,recbasis,alpha,vol,tol,long_deriv)
             if ((izpxh(izp(k))).or.(izpxh(izp(j)))) then
               xhgammahlp=.true.
             else
               xhgammahlp=.false.
             endif
             do lk=1,lmax(izp(k))
             do lj=1,lmax(izp(j))
               call shortrange1(r,basis,uhubb(izp(k),lk),
     &              uhubb(izp(j),lj),uhder(izp(k),lk),xhgammahlp,zeta,
     &              tol,SRdgdrkj,SRdgdr3kj)
               call shortrange1(r,basis,uhubb(izp(j),lj),
     &              uhubb(izp(k),lk),uhder(izp(j),lj),xhgammahlp,zeta,
     &              tol,SRdgdrjk,SRdgdr3jk)
               do ix=1,3
                 hgrad(ix,k) = hgrad(ix,k) + 0.5d0*qldiff(k,lk)
     &                  *qldiff(j,lj)*( SRdgdrjk(ix)+SRdgdrkj(ix)
     &                   +2.0d0*long_deriv(ix) )
                 if (sccmode=="3") then
                 hgrad(ix,k) = hgrad(ix,k) + qldiff(k,lk)/3.0d0
     &                  *qldiff(j,lj)*( qdiff(j)*SRdgdr3jk(ix)
     &                                 +qdiff(k)*SRdgdr3kj(ix) )
                 endif
               enddo
             enddo
             enddo
           endif
           enddo
           enddo

c gradient for non-periodic systems (cluster) and l-dependent Hubbard
         else 
           do k=1,nbeweg
           do j=1,nn
           if (j.ne.k) then
             r(1)=x(1,k)-x(1,j)
             r(2)=x(2,k)-x(2,j)
             r(3)=x(3,k)-x(3,j)
             norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
             if ((izpxh(izp(k))).or.(izpxh(izp(j)))) then
               xhgammahlp=.true.
             else
               xhgammahlp=.false.
             endif
             do lk=1,lmax(izp(k))
             do lj=1,lmax(izp(j))
               call gammaall1(norm,uhubb(izp(k),lk),uhubb(izp(j),lj),
     &              uhder(izp(k),lk),xhgammahlp,zeta,dgdrkj,dgdr3kj)
               call gammaall1(norm,uhubb(izp(j),lj),uhubb(izp(k),lk),
     &              uhder(izp(j),lj),xhgammahlp,zeta,dgdrjk,dgdr3jk)
               hlp = 0.5d0*qldiff(k,lk)*qldiff(j,lj)*(dgdrjk+dgdrkj)
               if (sccmode=="3") then
                 hlp = hlp + qldiff(k,lk)*qldiff(j,lj)/3.0d0
     &               *( qdiff(j)*dgdr3jk+qdiff(k)*dgdr3kj )
               endif
               do ix=1,3
                 hgrad(ix,k) = hgrad(ix,k) + hlp*r(ix)/norm
               enddo
             enddo
             enddo
           endif
           enddo
           enddo
            
         endif ! cluster/period


! now gradient for atomic Hubbards (original)
       else 
       if (period) then
         
         ! get dgamma_{kj}/dR_{kx} and dGamma_{kj}/dR_{kx} for periodic system
         ! Note that dgamma_{kj}/dR_{kx} = - dgamma_{kj}/dR_{jx} 
         ! Note that dGamma_{kj}/dR_{kx} = - dGamma_{kj}/dR_{jx} 
         ! Precalculation needs extra memory, however it should be faster because
         ! on-the-fly calculation would double the cost. 
         do k=1,nn
           do j=1,nn
           if (k.ne.j) then
             r(1)=x(1,k)-x(1,j)
             r(2)=x(2,k)-x(2,j)
             r(3)=x(3,k)-x(3,j)
             if ((izpxh(izp(k))).or.(izpxh(izp(j)))) then
               xhgammahlp=.true.
             else
               xhgammahlp=.false.
             endif
             call shortrange1(r,basis,uhubb(izp(k),1),uhubb(izp(j),1),
     &           uhder(izp(k),1),xhgammahlp,zeta,tol,SRdgdrkj,SRdgdr3kj)
             do ix=1,3
               SRdgdr(ix,k,j)  = SRdgdrkj(ix)
               SRdgdr3(ix,k,j) = SRdgdr3kj(ix)
             enddo
           endif
           enddo
         enddo
           
         ! calculate long-range contribution and built gammagradient 
         !contribution

         ! get reciprocal lattice vectors and cell volume
         call rezvol(basis,recbasis,vol)
         ! choose good convergence parameter alpha
         alpha = getalpha(basis)
         ! set tolerance for convergence
         tol = 1.0d-8
         do k=1,nbeweg
           do ix=1,3
             tmp(ix)=0.0d0
             tmp3(ix)=0.0d0
           enddo
           do j=1,nn
           if (k.ne.j) then
             r(1)=x(1,k)-x(1,j)
             r(2)=x(2,k)-x(2,j)
             r(3)=x(3,k)-x(3,j)
             call phi1(r,basis,recbasis,alpha,vol,tol,long_deriv)
             do ix=1,3
               tmp(ix)  = tmp(ix)  + qdiff(j)*( SRdgdr(ix,k,j)
     &                      -SRdgdr(ix,j,k)+2.0d0*long_deriv(ix)  )
               tmp3(ix) = tmp3(ix) + qdiff(j)*( qdiff(k)*SRdgdr3(ix,k,j)
     &                    -qdiff(j)*SRdgdr3(ix,j,k)     ) 
             enddo
           endif
           enddo
           if (sccmode=="3") then
             hgrad(1,k) = qdiff(k)*(0.5d0*tmp(1)+tmp3(1)/3.0d0)
             hgrad(2,k) = qdiff(k)*(0.5d0*tmp(2)+tmp3(2)/3.0d0)
             hgrad(3,k) = qdiff(k)*(0.5d0*tmp(3)+tmp3(3)/3.0d0)
           else
             hgrad(1,k) = qdiff(k)*0.5d0*tmp(1)
             hgrad(2,k) = qdiff(k)*0.5d0*tmp(2)
             hgrad(3,k) = qdiff(k)*0.5d0*tmp(3)
           endif
         enddo


       else ! no periodicity

         ! get dgamma/dr and dGamma/dr   (r=|R_j-R_k|)
         do k=1,nn
           do j=1,nn
             r(1)=x(1,k)-x(1,j)
             r(2)=x(2,k)-x(2,j)
             r(3)=x(3,k)-x(3,j)
             norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
             if ((izpxh(izp(k))).or.(izpxh(izp(j)))) then
               xhgammahlp=.true.
             else
               xhgammahlp=.false.
             endif
             call gammaall1(norm,uhubb(izp(k),1),uhubb(izp(j),1),
     &            uhder(izp(k),1),xhgammahlp,zeta,dgdrkj,dgdr3kj)
             dgdr(k,j)  = dgdrkj
             dgdr3(k,j) = dgdr3kj
           enddo
         enddo

         ! get dr/dR_{kx} and built gammagradient contribution
         do k=1,nbeweg
           do ix=1,3
             tmp(ix)=0.0d0
             tmp3(ix)=0.0d0
           enddo
           do j=1,nn
           if (j.ne.k) then
             r(1)=x(1,k)-x(1,j)
             r(2)=x(2,k)-x(2,j)
             r(3)=x(3,k)-x(3,j)
             norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
             do ix=1,3
               tmp(ix) = tmp(ix) + qdiff(j) * ( dgdr(j,k) + dgdr(k,j) )
     &                                           * r(ix)/norm
               tmp3(ix) = tmp3(ix) + qdiff(j)*( qdiff(j)*dgdr3(j,k) +
     &                    qdiff(k)*dgdr3(k,j) ) * r(ix)/norm
             enddo
           endif
           enddo
           if (sccmode=="3") then
             hgrad(1,k) = qdiff(k)*(0.5d0*tmp(1)+tmp3(1)/3.0d0)
             hgrad(2,k) = qdiff(k)*(0.5d0*tmp(2)+tmp3(2)/3.0d0)
             hgrad(3,k) = qdiff(k)*(0.5d0*tmp(3)+tmp3(3)/3.0d0)
           else
             hgrad(1,k) = qdiff(k)*0.5d0*tmp(1)
             hgrad(2,k) = qdiff(k)*0.5d0*tmp(2)
             hgrad(3,k) = qdiff(k)*0.5d0*tmp(3)
           endif
         enddo

       endif
       endif !no ldep

       end

