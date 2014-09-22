c---------------------------------------------------------------------      
c      Do the dirty work of reading from whichever file
c---------------------------------------------------------------------
      subroutine input(fileno,atoms,period,types,coords,vects,atnames)
c
c     fileno      = number of unit for input
c     atoms       = number of atoms
c     period      = logical (T) supercell, (F) cluster
c     types ()    = atom types
c     coord (,3)  = coordinates of atoms in Angstrom
c     vects (3,3) = supercell unit vectors

      
      implicit none
      include 'maxima.inc'
      
      integer fileno,atoms,types(NNDIM),i,j,xmo,iop,ientrn
      real*8 coords(3,NNDIM),vects(3,3),origen(3),xhlp(3),reada
      logical period
      character*70 atnames
      character*2 atyp(maxtyp),tmpchr
      character*(mxline) line
      character*(maxldt) str(mxentr)     
      integer nco,nc(NNDIM,2)
      real*8 r0yz,r0x,kspr
c      ME: for constraints like in CHARMM
      integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
      real*8 rccval,kccval,rccinc
c      common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc

c
      common /correct/ r0yz , r0x ,kspr, nco, nc
      common /atosym/ xmo,atyp
c     ME initialize
      lcc=0
 
      do i=1,NNDIM
        types(i)=0
        do j=1,3
          coords(j,i)=0.0d0
        enddo
      enddo
      do i=1,3
        do j=1,3
          vects(i,j)=0.0d0
        enddo
      enddo

 10   read(fileno,'(a)',end=100) line 
      call parse(line,mxentr,mxline,str,ientrn,iop)
      if(iop.eq.1) goto 10 
      atoms  =nint(reada(str(1),1))
      tmpchr =str(2)
c     get xmo, every xmo iteration a xmol file is written
      if(ientrn.lt.3) then
        xmo=100000
      else 
        xmo=nint(reada(str(3),1))
        if(xmo.eq.0) xmo=100000
      endif

 20   read(fileno,'(a)',end=100) line
      atnames =line
      call parse(line,mxentr,mxline,str,ientrn,iop)
      if(iop.eq.1) goto 20 
      if(ientrn.gt.maxtyp) then
        write(*,'(1x,a,i3,a)') 'Number of atomtypes exceeds',maxtyp,' !'
        stop
      endif
      do i=1,ientrn 
        atyp(i)=str(i)
      enddo 

      if( (tmpchr.eq.'s').or.(tmpchr.eq.'S').or.(tmpchr.eq.'SM') ) then
       period = .TRUE.
      else
       period = .FALSE.
      endif
      do i=1,atoms
       read(fileno,*)j,types(i),(xhlp(j),j=1,3)
       do j=1,3
        coords(j,i)=xhlp(j)/0.529177d0
       enddo
      enddo
      
      if(period) then
       read(fileno,*) (origen(j),j=1,3)
      
       do i=1,3
         read(fileno,*)(xhlp(j),j=1,3)
         do j=1,3
          vects(i,j)=xhlp(j)/0.529177d0
         enddo
       enddo
      endif
c     ME: constraints like in CHARMM
      if ((tmpchr.eq.'M').or.(tmpchr.eq.'SM')) then
      read(fileno,*) lcc,kccval,rccval,rccinc
      read(fileno,*) (jcc(i),atcc1(i),atcc2(i),i=1,lcc)
      rccval = rccval/0.529177d0
      rccinc = rccinc/0.529177d0 
      kccval = kccval*0.529177d0
      write(*,*) 'const with:(N,[AA],[H/AA] ',
     c   lcc,rccval,kccval
      endif

ccc add something emprically 
      if (tmpchr.eq.'K') then
       read(fileno,*) nco
       do i=1,nco
         read(fileno,*) (nc(i,j),j=1,2)
       enddo
      endif
ccc add restoring potential:
      if (tmpchr.eq.'R') then
      write(*,*) 'nco,xy-zylinder,z-Zylinder,spring,[kcal]'
        read(fileno,*) nco,r0yz,r0x,kspr
        write(*,*) 'restraining potential used for first ', -nco, 
     c   ' atoms'
      endif
cccccccccccc
      return

 100  print *,'Unexpected end of structure file !'
      stop

      end
