      subroutine dispersionread(nn,ntype,izp,x)
      include 'maxima.inc'
      integer nei(NNDIM),izp(NNDIM)
      real*8 Ni0(NTYPE),Ni(NNDIM)
      real*8 x(3,NNDIM),dif(3),r2
      real*8 h1(MAXTYP,4),h2(MAXTYP,4),Edis
      real*8 R0vdw(MAXTYP,4),C6(NNDIM,NNDIM)
      real*8 Rvdw(NNDIM,NNDIM)
      real*8 hh1(NNDIM),hh2(NNDIM) 
      real*8 C,A,B,r0,rv,scale
      logical read
      common /disper/ Edis, dispers
      common /dispertmp/ A,B,C,r0,rv,C6,Rvdw
      read = .false.
      open(54,file="DISPERSION.INP",status="unknown")
c      read(54,*) A,B,C,rv,r0,scale
c     parameters as in Elstner et al., JCP 2000
      A=7.0d0
      B=4.0d0
      C=-3.0d0
      rv=0.0d0
      r0=0.0d0
      scale=1.0d0
c
      if (scale .le. 0.0) then
      write(*,*) 'London type dispersion, pol and dis in eV,scale='
     c ,scale
      else
       write(*,*)  'Slater-Kirkwood dispersion switched on:'
      endif

      do i=1,ntype 
      if ( scale .ge.0.0 ) then
c  $1-$4 = static atomic polarizability in Angstrom^3 for element with
c          0,1,2 or 3 Hydrogen neighbors
c  $5-$8 = atomic vanderWaals radius in Angstrom
c  $9    = slater-kirkwood effective number of electrons 
c          0.8 for H and 1.17+0.33*(numberofvalenceelectrons) for others
       read(54,*,end=10) (h1(i,j),j=1,4),(h2(i,j),j=1,4),Ni0(i)
       else
       read(54,*,end=10) (h1(i,j),j=1,4),(h2(i,j),j=1,4)
       endif
      enddo 
       write(*,*) 
     c '  read C6, Rvdw (eV,A) from DISPERSION.INP'
      open(15,file='DISP.CHEK')
      write(*,*) '  careful, parameter determined from # of H atoms'

      goto 20
10    continue
       read = .true.
       open(16,file='DISP.INP')
         write(*,*) '  DISPERSION.INP empty, read from DISP.INP' 
20    continue
c if we read from DISPERSION.INP:
c  determine hybridization from number of Hydrogens
c wir starten mit 1 Nachbarn. Dann zaehlen wir die
c Anzahl der Wasserstoffe durch: fuer jeden Wasserstoff gibt es
c einen Nachbarn mehr. Die werte fuer C und N unterscheiden sich nur in den
c Waserstoffen. N mit 2H ist anders als nur N oder N mit 1H
c C ist nur fuer 3H unterschiedlich
c
c  determine parameters for all atoms
      do i=1,nn
       if (read)  then
         read(16,*)  hh1(i),hh2(i),Ni(i) 
       else 
         nei(i) = 1
          do j=1,nn
           r2 = 0.0d0
           if ( j.ne.i) then
            do n = 1,3
             dif(n) = x(n,i) - x(n,j)
             r2 = r2 + dif(n)**2
            enddo
            r = dsqrt(r2)*0.529177d0
            if (r.le.1.2) then
             nei(i) = nei(i) +1
            endif
           endif
          enddo !j=1,nn
         if (nei(i).gt.4) then
           hh1(i)=h1(izp(i),4)
           hh2(i)=h2(izp(i),4)
         else  
           hh1(i)=h1(izp(i),nei(i))
           hh2(i)=h2(izp(i),nei(i))
         endif
         Ni(i) = Ni0(izp(i))
         write(15,'(3F12.6)')  hh1(i),hh2(i),Ni(i)
       endif
   

c   check values
       if ((hh1(i).eq. 0.0).or.(hh2(i).eq.0.0).or.(Ni(i).eq.0.0)) then
        write(*,*) 'a parameter is 0.0 for atom',i,izp(i),nei(i)
        stop
       endif
c
      enddo !i=1,nn
c set up mixed coefficients 
c mixing from Halgren JACS 1992 114 p7827
       if (.not. read)  then
       write(15,*) ' --------------'
       write(15,*) ' I J  typeI typeJ C6 R NeiI NeiJ'
       endif
      do i=1,nn
       do j=1,i
        if (scale .le. 0.0) then
         C6(i,j) = -scale*1.5d0*hh1(i)*hh1(j)*
     c   hh2(i)*hh2(j)/
     c   ( hh2(i)+hh2(j) )
        else
cccc  17.532 conversion from eV in a.u. for polarizability
c   0.5975 conversion from [H au**6] to [eV A**6]
c  total 17.532*0.5975 = 10.476
         C6(i,j) = scale*1.5d0*10.476d0*hh1(i)*hh1(j)/
     c   ( dsqrt(hh1(i)/Ni(i) ) + 
     c    dsqrt( hh1(j)/Ni(j) ) ) 
        Rvdw(i,j)=(hh2(i)**3 + hh2(j)**3)/
     c    (hh2(i)**2 + hh2(j)**2 )
        Rvdw(j,i) =  Rvdw(i,j)
        endif
        C6(j,i) = C6(i,j) 
       if (.not. read)  then
        write(15,'(4I4,2F12.6,2I4)') i,j,izp(i),izp(j),
     c   C6(i,j),Rvdw(i,j),nei(i),nei(j)
       endif
       enddo
      enddo
c      write(*,*) 'end reading disper'
      end

