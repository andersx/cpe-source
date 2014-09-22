       subroutine externalshift(nn,x,izp,shiftE)
       implicit REAL*8 (A-H,O-Z)
       include 'maxima.inc'
       integer n,nE,izp(NNDIM),nEiter
       real*8 x(3,NNDIM)
c      real*8 xE(3,4*NNDIM),ZE(4*NNDIM)
       real*8 xE(3,EXTDIM),ZE(EXTDIM)
c       real*8 xE(3,NNDIM),ZE(NNDIM)
       real*8 shiftE(NNDIM),dif(3),r,r2
       character*2 EXT
       logical ldep
       common /mcharge/ qzero(MAXTYP,4), uhubb(MAXTYP,3), ldep
       common /extchr/ xE, ZE, nE,  nEiter,icycle, EXT

c       const=1/0.529177
       if (EXT.eq.'CH') then
        do j=1,nn 
        shiftE(j) = 0.0
        do k=1,nE
         r2=0.0
         do i=1,3 
          dif(i) = x(i,j) - xE(i,k)
          r2=r2 + dif(i)**2
         enddo 
c          uhub=uhubb(izp(j))
c          gamma = 1.0/sqrt(r2 + (0.5/uhub + 0.5/uhub)**2)
         gamma =  1.0/sqrt(r2)
         shiftE(j) = shiftE(j) + gamma*ZE(k)
        enddo 
        enddo 
       endif 
c end CH!, now field EF
       if (EXT.eq.'EF') then
       if (icycle .ge. nEiter ) then
       do i=1,nn
        shiftE(i) = 0.0
        do j=1,3
         shiftE(i) = shiftE(i) + ZE(j)*(x(j,i)-xE(1,j))
        enddo 
       enddo
       endif
       endif 
      end
