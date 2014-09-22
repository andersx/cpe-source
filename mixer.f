c purpose: organize charges to be fed to broyden-mixer

       subroutine mixer(niter,almix,nn,qmold,qmat,spin,ndim,lldim,ldep,
     &     lmax,izp,indl,qlold,ql,qlupold,qldownold,qlup,qldown)
       implicit none
       include 'maxima.inc'
       logical spin,ldep
       integer niter,nn,ndim,lldim,lmax(MAXTYP),izp(MAXTYP)
       integer indl(NNDIM+1)
       real*8  almix,qmold(NNDIM),qmat(NNDIM),ql(3*NNDIM),qlold(3*NNDIM)
       real*8  qlup(3*NNDIM),qldown(3*NNDIM),qlupold(3*NNDIM)
       real*8  qldownold(3*NNDIM)

       real*8  qmix(6*NNDIM),qmixold(6*NNDIM)
       integer i,li,mi,mu,l

c spinpolarized (change all q_al) 
       if (spin) then
         do i=1,lldim
           qmix(i)       = qlup(i)
           qmix(i+lldim) = qldown(i)
           qmixold(i)       = qlupold(i)
           qmixold(i+lldim) = qldownold(i)
         enddo 
         call broyden(niter,almix,2*lldim,qmixold,qmix)
         do i=1,lldim
           qlup(i)   = qmixold(i)
           qldown(i) = qmixold(i+lldim)
           qlupold(i)   = qmixold(i)
           qldownold(i) = qmixold(i+lldim)
           ql(i)     = qlup(i)+qldown(i)
         enddo
         do i = 1,nn
           qmat(i)=0.0d0
           do li = 1,lmax(izp(i))
             qmat(i) = qmat(i)+qlup(indl(i)+li)+qldown(indl(i)+li)
           enddo
         enddo
c spinunpolarized
       else
c and l-dependence Hubbard
         if (ldep) then
           call broyden(niter,almix,lldim,qlold,ql)
           do i=1,lldim
             ql(i) = qlold(i)
           enddo
           do i = 1,nn
             qmat(i)=0.0d0
             do li = 1,lmax(izp(i))
               qmat(i) = qmat(i)+ql(indl(i)+li)
             enddo
           enddo
c no l-dependent Hubbard (standard DFTB2/3)
         else
           call broyden(niter,almix,nn,qmold,qmat)
           do i = 1,nn
             qmat(i)= qmold(i)
           end do
         endif
       endif

       end
