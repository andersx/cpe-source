       subroutine externalchgrad(nn,x,izp,gr,qmat)
       implicit REAL*8 (A-H,O-Z)
       include 'maxima.inc'
       integer n,nE,izp(NNDIM),nEiter
       real*8 x(3,NNDIM)
c       real*8 xE(3,NNDIM),ZE(NNDIM)
c      real*8 xE(3,4*NNDIM),ZE(4*NNDIM)
       real*8 xE(3,EXTDIM),ZE(EXTDIM)
       real*8 dif(3),r,r2,gr(3,NNDIM),dgr
       real*8  qmat(NNDIM)
       character*2 EXT
       logical ldep
       common /mcharge/ qzero(MAXTYP,4), uhubb(MAXTYP,3), ldep
       common /extchr/ xE, ZE, nE,  nEiter,icycle, EXT
c       const=1/0.529177
       do j = 1,nn
        do i=1,3
         gr(i,j) = 0.0 
        enddo
       enddo
       if (EXT.eq.'CH') then
       do j = 1,nn  
        izpj = izp(j)
        do k =  1,nE
        r2 = 0.0
        do i = 1,3 
         dif(i) = x(i,j) - xE(i,k)
         r2 = r2 + dif(i)**2
c         write(*,*) x(i,j), xE(i,k)
        enddo 
         r=sqrt(r2)
c         write(*,*) r
c         uhub=uhubb(izp(j),1)
c        gamma= gamE(r2,uhubb(izpj),uhubb(izpj))
c          gamma = 1.0/sqrt(r2 + (0.5/uhub + 0.5/uhub)**2)
         gamma=1/r
        do l=1,3 
         dgr = dif(l)*(gamma**3)*
     c   (qmat(j)-qzero(izpj,4))*ZE(k)
         gr(l,j) = gr(l,j) + dgr
        enddo 
       enddo 
      enddo 
      endif
c   external field
       if (EXT.eq.'EF') then
        if (icycle .eq. nEiter ) then
        write(*,'(A20,3F12.6)') 'external field E= ',(ZE(i),i=1,3)
        endif
        if (icycle .ge. nEiter ) then
        do j = 1,nn
         izpj = izp(j)
          do i=1,3 
           dgr= -(qmat(j)-qzero(izpj,4))*ZE(i) 
           gr(i,j) = gr(i,j) + dgr
          enddo
         enddo
         endif
       endif
      end
