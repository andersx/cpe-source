# HIer fuegen wir den zmatrix-input ein, mit einfrieren von
#  koordinaten das wird dann mode 9


#      PROGRAM DYLCAO
#     ================
#
#     Copyright 1991 by Peter Blaudeck, Dirk Porezag 
#
# *********************************************************************
#
#     PROGRAM CHARACTERISTICS 
#     -----------------------
#
# DYLCAO calculates the dynamics of various systems within a 
# two-centre SETB formalism 
#
# *********************************************************************
#
# maxima definitions
#
implicit REAL*8 (A-H,O-Z)
include 'maxima.h'
#
integer izp(NNDIM), lmax(MAXTYP),irun,maxrun
integer dim(MAXTYP,MAXTYP), conat(NNDIM+1)
integer numint(MAXTYP,MAXTYP)
integer nexc,nrem,bfgserror
real*8 x(3*NNDIM), xparam(3*NNDIM), xold(3*NNDIM), gr(3*NNDIM), norm
real*8 skhtab(10,MAXTAB,MAXTYP,MAXTYP), skstab(10,MAXTAB,MAXTYP,MAXTYP)
real*8 skself(3,MAXTYP), sr(MAXTYP,MAXTYP), nel, convec(3,NNDIM)
real*8 coeff(6,MAXINT,MAXTYP,MAXTYP),xr(2,MAXINT,MAXTYP,MAXTYP)
real*8 efkt(3,MAXTYP,MAXTYP),cutoff(MAXTYP,MAXTYP)
real*8 xm(MAXTYP), xmrc(MAXTYP), rmax(MAXTYP,MAXTYP)
real*8 ar(10,MAXTYP,MAXTYP),qmat(NNDIM),slkcutoff
real*8 espin(MAXTYP),ekin
real*8 xhlp(3*NNDIM)
character*65 infile, outfile, name, tmpchr, namev
character*60 smode(16)
character*70 atnames
character*1  sccmode,sccmodeS
logical mdtemp, mdfree, stdc, cgrmin,bfgsmin, atomic, evector,bcalcvib
logical converge, mdcorr, dorelax, doecalc
logical chrr,constr,ex
logical lbfgs,lzmat,lgrad,readco
logical lexc
# K.W. 2014-04-14
logical excited
external skhpar, skspar
#CHARMM
# real*8 xE(3,4*NNDIM),ZE(4*NNDIM)
real*8 xE(3,EXTDIM),ZE(EXTDIM)
#integer nE,nEiter
character*2 EXT
#      ME: for constraints like in CHARMM
      integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
      real*8 rccval,kccval,rccin
character*5 dummystr
logical xhgamma,izpxh(MAXTYP),fudge,dummyl,ldep
real*8  kl1,a0fudge,bfudge,cfudge
real    cputime_start,cputime_stop
integer clock_count_start,clock_count_stop,clock_rate,clock_max
# MG:120412  for collinear spin-formalism (JPCA2007,111,5622) (s,p,d-orbital)
real*8  nelunp,wspin(MAXTYP,3,3),nelup,neldown
real*8  ql(3*NNDIM),qlup(3*NNDIM),qldown(3*NNDIM)
logical spin
integer lcount
# K.W. 2013-06-20 CM3 calculation
logical docm3
real*8  qcm3(NNDIM)
# K.W. 2014-04-14 CPE calculation
logical docpe

#DISPERSION
# integer subsys(NNDIM)
# real*8 pol(MAXTYP,4),ion(MAXTYP,4),Edis,C6(NNDIM)
real*8 Edis
logical dispers
common /disper/  Edis, dispers
#END DISPERSION
#CHARMM
#common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc
common /extchr/ xE, ZE, nE, nEiter ,icycle, EXT
common /sktab/  skhtab,skstab,skself,sr,dim
common /spin/   espin,nelup,neldown,wspin,spin
common /spltab/ coeff,xr,efkt,cutoff,numint
common /concgr/ convec, constr, conat
common /fudgev/ fudge,dummyl,a0fudge,bfudge,cfudge
#
###### 
# These data are read in by gettab. 
# Every pair of atoms requires a table containing the Slater-Koster 
# data. The first line consists of the distance stepwidth and the
# number of distance values. If the file describes homonuclear
# interactions, an additional second line contains the energies
# of the valence d, p and s electrons (diagonal elements of the
# hamiltonian). 
# All other lines contain 10 hamiltonian and 10 overlap numbers
# in the order dds ddp ddd pds pdp pps ppp sds sps sss.
# The first of these lines (r = 0) contains:
# element 01:    mass of the atom (0 for heteroatomic files)
# element 02-10: dummy-variables
# All energies in H = Hartrees = 4.358d-11 g*cm**2/s**2
# All masses in nucleon masses = 1.6603d-24 g
#
common /lmax/ lmax
common /izp/ izp, nel, nbeweg
common /rdaten/ ar
common /machine/ racc
# in commonbox are things like boxsiz,...
include 'commonbox.h'
common /mcharge/ qzero(MAXTYP,4), uhubb(MAXTYP,3), ldep
common /uhq/ uhder(MAXTYP,3)
common /scchelp/ sccmode
common /hxmod/   kl1,xhgamma,izpxh

call cpu_time(cputime_start)
call system_clock ( clock_count_start, clock_rate, clock_max )

data smode /'MD relaxation with heat reservoir   ',
            'free MD without external constraints',
            'steepest descent relaxation         ',
            'conjugate gradient relaxation       ',
            'bfgs relaxation',
            'Mulliken charge and atomic energy calculation',
            'Option number 6 + prints out the eigenvectors',
            'BFGS relaxation in internal coordinates',
      	    'Option number 6 + Z-matrix input',
            '',
            '',
            '',
            'conjugate gradient relaxation + writing phonir input',
            '',
            'BFGS relaxation + writing phonir input',
            'Energy evalulation with fixed HOMO-LUMO occupation' /
#
lcc=0 #????
# print *,'** dftb (version 09.03.2011) **'
print *, '** DFTB (version 09.05.2014) ** '
print *, '** pimped with CM3, CPE, fixed occupation number ** '
#
# Get the machine accuracy
#
racc = 1.0d0
while ((1.0d0+racc) > 1.0d0) racc = racc*0.5d0
racc = 2.0d0*racc 
#
# Enter relaxation mode and fmax and scftol. Valid relaxation modes are
# 1 ... MD with scaling of velocities according to temperature
# 2 ... MD without scaling of velocities according to temperature
# 3 ... Steepest descent (velocities are set to zero after each step)
# 4 ... Conjugate gradient relaxation
# 5 ... bfgs relaxation
# 6 ... Mulliken analysis and atomic energy calculation
# 7 ... Option number 6 + prints out the eigenvectors
# 8 ... BFGS relaxation in internal coordinates
# 9 ... Option number 6 + Z-matrix input
# 10... Excitation spectrum via linear response --> not available in this version
# 11... Option 1, read old coordinates
# 12... Option 2, read old coordinates
# 13... Option 4 + print input for phonir 
#       (calculate energy,force,dipole,... to use for build the HESSIAN
#       with finite differences)
# 15... Option 5 + print input for phonir 
#       (calculate energy,force,dipole,... to use for build the HESSIAN
#       with finite differences)
# K.W.: 2014-04-14
# 16... energy calculation of HOMO to LUMO single excited state
#
# If the total force acting on each atom in the structure is smaller 
# than fmax, the conjugate gradient or steepest descent routine is 
# converged and the program terminates.
# scftol is the convergence criterion for the charge SCF convergence.
#
#
mdtemp = .false. ; mdfree = .false. ; stdc = .false. ; cgrmin = .false.
atomic = .false. ; evector = .false.; constr =.false.; lbfgs = .false.
lzmat = .false. ; lgrad = .true. ; lexc = .false. ; bfgsmin = .false.
readco = .false. ; fudge = .false.; bcalcvib = .false.; 
xhgamma = .false. ; excited = .false.
docm3   = .true. ; docpe   = .false. 
bfgserror=0; ldep=.false.; spin=.false.
print *,'enter mode, fmax, scfmode, scftol, read charge,dispers,EXT'
read *,mode,fmax,sccmode,scftol,chrr,dispers,EXT
print *,'mode, fmax, scfmode, scftol, read charge,dispers,EXT set to'
print *,mode,fmax,sccmode,scftol,chrr,dispers,EXT
if(sccmode == "1") scftol=-abs(scftol)
else scftol=dabs(scftol)
if (sccmode == "1")     print *,'+++ DFTB1 (PorezagPRB1995)'
else if(sccmode == "2") print *,'+++ DFTB2, formerly called SCC-DFTB (ElstnerPRB1998)'
else if(sccmode == "3") print *,'+++ DFTB3 (GausJCTC2011)'
else if(sccmode == "C") print *,'+++ DFTB3 (GausJCTC2011) with CPE correction (KaminskiJPCA2012)'
else if(sccmode == "S") print *,'detailed specifications: enter order taylor series (1,2,3) , gamma^H(T,F) , fudge(T,F), l-depHubb(T,F), nelunp(real)'
else {
  print *,'scc-mode "',sccmode,'" is not available'
  stop
}
if (sccmode == "C") {
  # The 'C' stands for CPE correction
  sccmode = "3" 
  docpe = .true.
  docm3 = .false.
}
if (sccmode=="3") xhgamma=.true.
if(sccmode=="S"){
  read *,sccmodeS,xhgamma,fudge,ldep,nelunp
  if      (sccmodeS=="1") {
    sccmode="1"
    print *,'Scc=off (PorezagPRB1995)'
  }
  else if (sccmodeS=="2") {
    sccmode="2"
    print *,'2ndOrderScc=on (ElstnerPRB1998)'
  }
  else if (sccmodeS=="3") {
    sccmode="3"
    print *,'3rdOrderScc=on'
  }
  else {
    print *,'no such sccmodeS found, exit program'
    stop
  }

  if(xhgamma){
    xhgamma=.true.
    print *,'XHgamma=on'
  }
  else {
    xhgamma=.false.
    print *,'XHgamma=off'
  }

  if(fudge){
    print *,'fudge-term=on (YangJCTC2008)'
    print *,'enter  "fudge",a0fudge,bfudge,cfudge'
    read *,dummystr,a0fudge,bfudge,cfudge
    if (dummystr != 'fudge'){
      print *,'  no fudge variables specified, exit program'
      stop
    }
    print *,'  a0fudge = ',a0fudge,'  bfudge = ',bfudge,'  cfudge = ',cfudge
  }
  else print *,'fudge=off'
  if (ldep)  print *,'l-shell dependent Hubbards=on'
  else       print *,'l-shell dependent Hubbards=off'
  if ((ldep).and.(fudge)) {
    print *,'l-shell dependent Hubbard not implemented for fudge-term (YangJCTC2008)'
    stop
  }
  if (nelunp>1.0d-12){
    if(sccmode=="1") {
      print *,'ERROR: spin-polarization formalism not available in non-SCC mode'
      stop
    } 
    spin=.true.
    print *,'Collinear spin-formalism will be applied (JPCA2007,111,5622), number of unpaired electrons: ',nelunp
  } 
} # endif sccmode=="S"


#if (mode <= 0 ){
# mode = -mode
# constr=.true.
#}
if (mode == 11 ){
 readco = .true.
 mode = mode - 10
}
if (mode == 12 ){
 readco = .true.
 mode = mode - 10
}
if (mode == 1) mdtemp = .true.
else if (mode == 2) mdfree = .true.
else if (mode == 3) stdc = .true. 
else if (mode == 4) {
 cgrmin = .true.
 # evector = .true.
}
else if (mode == 5) {
 bfgsmin = .true.
 # evector = .true.
}
else if (mode == 6) {
 atomic = .true.
 stdc   = .true.
}
else if (mode == 7) {
 atomic  = .true.
 stdc    = .true.
 evector = .true.
 }
else if (mode == 8) lbfgs = .true.
else if (mode == 9) {
 lbfgs = .true.
 lzmat = .true.
 }
else if (mode == 10) {
 print *,'mode 10 not available in this dftb-version, stop program'
 goto 9
 lgrad = .false.
 lexc  = .true.
}
else if (mode == 13) {
  cgrmin   = .true.
  bcalcvib = .true.
}
else if (mode == 15) {
  bfgsmin  = .true.
  bcalcvib = .true.
}
else if (mode == 16) {
  excited     = .true.
  print *, 'energy evaluation with fixed HOMO/LUMO occupation number'
}
else {
  print *,'invalid relaxation mode' 
  goto 9 
}
write(*,'(1x,2a)') 'relaxation mode: ',smode(mode)
#
# Read in name of structure file
#
print *,'enter filename for input structure'
read *,infile
# 
# All lengths in a.u. = 5.292d-9 cm
#
print *,'infile :',infile
open (1,file=infile,status='old'); rewind 1

write(*,*) 'Input total charge'
read(*,*) nel
write(*,*) 'Total charge set to ',nel
write(*,*) 'Number of moveable atoms'
read(*,*) nbeweg
if(nbeweg==0) {
  write(*,*) 'No constraints set.'
}
else {
  write(*,*) 'Number of moveable atoms set to ',nbeweg
}
if(nbeweg<0) {
 constr = .true.
 write(*,*) nbeweg,'atoms to be constraint:'
# read(*,*) conat(1)
 conat(1) = -nbeweg
 nbeweg = 0
 write(*,*) 'Enter ',conat(1),' times: atom-number constraint-vector (x y z) :'
 do i=1,conat(1) {
  read(*,*) conat(i+1),(convec(j,i),j=1,3)
#put other constaraints: fix x,y,z Koordinate individually!
# also changed end of eglcao.f
   norm = 1.0d0
#  norm=convec(1,i)*convec(1,i)+convec(2,i)*convec(2,i)+convec(3,i)*convec(3,i)
  norm=sqrt(norm)
  convec(1,i)=convec(1,i)/norm
  convec(2,i)=convec(2,i)/norm
  convec(3,i)=convec(3,i)/norm
 }
}

#
# Read in coordinates and box parameters
#
#if(lzmat) {
# call zinput(1,nn,period,izp,xparam)
#}
#else {
 call input(1,nn,period,izp,x,boxsiz,atnames)
#}
if (nn > NNDIM) {
  print *,'dftb: nn > ', NNDIM
  goto 9
}
if(nbeweg==0){
nbeweg=nn
}
n3 = 3*nn
nbew3 = 3*nbeweg

#
# calculate Inverse of Boxmatrix
#
if(period) {
 call inversebox
}

close (1)
#
# Check if ntype in limits of MAXTYP 
#
ntype = 1
do n = 1,nn {
  if (izp(n) > ntype) ntype = izp(n)
}
if (ntype > MAXTYP) {
  print *,'dftb: ntype >', MAXTYP
  goto 9
}
#
# Read name of output file
#
write (*,*) 'enter filename for output structure'
read *,outfile
write(*,*) 'outfile: ',outfile

#
# Read in whether atom shall have a virtual atom at the same place,
# i.e. wether you want an additional s (or p) orbital or not
#

n3 = 3*nn
nbew3 = 3*nbeweg

#
# Read in maximal angular momentum for each type
# 1,2,3  for  s,p,d, respectively
#
write (*,*) 'enter ',ntype,' * lmax (for each atom type and sorted as in ',infile,')'
read *,(lmax(i),i = 1,ntype)
write (*,*) 'lmax set to: ',(lmax(i),i = 1,ntype)
#
# Read Hubbard derivative parameters and zeta for gamma^h
#
  if (ldep) {
    print *,'enter ',ntype,' lines (one for each atom type sorted as in ',infile,') of Hubbard derivatives Udd, Udp, Uds'
    do i=1,ntype {
      read *,(uhder(i,4-j),j=1,3)
      print *,(uhder(i,4-j),j=1,3)
    }
    if (sccmode!="3") print *,'Hubbard derivatives not used for this sccmode'
    }
  else { # no ldep
    write (*,*) 'enter ',ntype,' * dU/dq (for each atom type and sorted as in ',infile,')'
    read *,(uhder(i,1),i=1,ntype)
    if (sccmode=="3") write (*,*) 'dU/dq set to: ',(uhder(i,1),i = 1,ntype)
    else {
      write(*,*) 'values for dU/dq not used for this sccmode' 
      do i=1,ntype { uhder(i,1)=0.0d0 }
    }
  }
  write(*,*) 'enter zeta-parameter for gamma^h'
  read *,kl1
  if (xhgamma) write (*,*) 'zeta for gamma^h set to: ',kl1
  else         write (*,*) 'zeta not used for regular gamma function'
#
# Read wspin matrix for each atom type: Wss, Wsp, Wps, Wpp, Wsd, Wpd, Wdd, Wds, Wdp
#
  if (spin){
    print *,'enter ',ntype,' lines (one for each atom type sorted as in ',infile,') of atomic spin constants'
    print *,'Wss Wsp Wps Wpp Wsd Wpd Wdd Wds Wdp'
    do i=1,ntype {
      read  *,wspin(i,1,1),wspin(i,1,2),wspin(i,2,1),wspin(i,2,2),wspin(i,1,3),wspin(i,2,3),wspin(i,3,3),wspin(i,3,1),wspin(i,3,2)
      write (*,'(10F13.8)') wspin(i,1,1),wspin(i,1,2),wspin(i,2,1),wspin(i,2,2),wspin(i,1,3),wspin(i,2,3),wspin(i,3,3),wspin(i,3,1),wspin(i,3,2)
    }
  }
  
#
# Read in Slater-Koster-table 
#
call gettab(ntype)
#
# nel is the net charge
# nel=1, one electron missing
# nel =-1 one electron more
 nel = -nel
 do i = 1,nn {
   nel = nel+qzero(izp(i),4)
 }
write(*,*) 'total number of electrons:', nel
if (spin){
  nelup = 0.5d0 * (nel + nelunp)
  neldown = 0.5d0 * (nel - nelunp)
  write(*,*) 'total number of up-electrons:', nelup
  write(*,*) 'total number of down-electrons:', neldown
}
#
# set charges
if(sccmode .ne. "1") {
if(!chrr) {
 lcount=0
 do i = 1,nn {
   qmat(i) = qzero(izp(i),4)
   do j=1,lmax(izp(i)){
     lcount=lcount+1
     ql(lcount)     = qzero(izp(i),j)
     qlup(lcount)   = qzero(izp(i),j)/2.0d0
     qldown(lcount) = qzero(izp(i),j)/2.0d0
   }
 }
}
else {
 if(spin){
   write(*,*) "ERROR: read-charges and combination with spin-polarization-formalism not implemented"
   stop
 }
 if(ldep){
   write(*,*) "ERROR: read-charges and combination with l-dependent Hubbards not implemented"
   stop
 }
 INQUIRE(file="CHR.DAT",exist=ex)
 if (ex) {
 open(37,file="CHR.DAT",status="unknown")
 do i=1,5 {
  read(37,*) tmpchr
 }
 do i=1,nn {
  read(37,*) j,qmat(i)
  if ((qmat(i)-qzero(izp(i),4))>1.0) {
   write(*,*) 'warning: check charge on atom',i
  }
 }
 write(*,*) "read initial charges from CHR.DAT"
 close(37)
 } else
 do i = 1,nn {
   qmat(i) = qzero(izp(i),4)
 }

}
}

#
# Get and check masses, data for repulsive part and next neighbours 
#
do i = 1,ntype {
  ar(1,i,i)=skhtab(1,1,i,i)
  xm(i) = ar(1,i,i)
  if ((! cgrmin)&&(!bfgsmin)) {
    if (xm(i) < 0.1d0) {
      write (*,*) 'mass of sort ',i,' is too small: ',xm(i)
      goto 9
    }
    else xmrc(i) = 1.0d0/(1822.887428d0*xm(i))
  }
  # get izpxh(MAXTYP): criteria for which elements gamma^h will be used
  if (xhgamma && xm(i)<3.5d0) {
    izpxh(i)=.true.
    print *,'+++ gamma^h used for combinations with atom type ',i,' (criteria: mass < 3.5 u)'
  }
  else izpxh(i)=.false.

  do j = 1,ntype {
    do k = 2,10 {
      ar(k,i,j) = skhtab(k,1,i,j)
    }
    rmax(i,j) = skstab(1,1,i,j)
  }
}

#
#DISPERSION
if(dispers) {
 call dispersionread(nn,ntype,izp,x)
} 
#END DISPERSION
#CHARMM
if (EXT=='CH') { 
 open(98,file='EXTCHARGES.INP')
 write(*,*) 'EXTCHARGES.INP'
 write(*,*) 'number of external charges'
 read(98,*) nE
 write(*,*) nE,' times x,y,z,charge'
 do i=1,nE {
   read(98,*) (xE(j,i),j=1,3),ZE(i)
   do j=1,3 {
    xE(j,i)=xE(j,i)/0.529177d0
   }
 }      
}
if (EXT=='EF') {
# same variables a for extcharges!
#nEiter: ionic step, at which Field
# is switched on
 open(98,file='EXTFIELD.DAT')
 write(*,*) 'EXTFIELD.DAT'
 read(98,*) nE,nEiter
 write(*,*) '3 , icycle: when to switch field=',nEiter
 do i=1,nE {
   read(98,*) xE(1,i),ZE(i)
   write(*,'(A20,2F12.6)') 'x0 , Ex = ',xE(1,i),ZE(i)
    xE(1,i)=xE(1,i)/0.529177d0
 }  
}
#END CHARMM
# Check boxsize according to the cutoff in SLK-data:
#    the shortest distance of two vertices has to be twice the shortest
#    cutoff
#

  slkcutoff=0.0d0
  do i = 1,ntype {
    do j = 1,ntype {
     yhlp=sr(i,j)*dim(i,j)+0.3d0
     if(yhlp>slkcutoff) slkcutoff=yhlp
   }
  }
  skcut2=slkcutoff*slkcutoff

if (period) { 
  call gamma_summind(slkcutoff)
  write(*,*) 'Number of lattice sums for Gammapoint matrix:', nlat(1), nlat(2), nlat(3)
  call initphi(boxsiz)
}
# Initialize random generator
# Set up initial velocities and old coordinates
# Take care: frozen coordinates must have xold(i) = x(i)
#
inr1 = 1802; inr2 = 9373; call rmarin(inr1,inr2)
if (!readco) {
do i = 1,n3 {
  xold(i) = x(i)
}
if (mdtemp) {
  do i = 1,nbew3 {
    xold(i) = x(i) - 0.1d0*(ranmar()-0.5d0)
  }
}
}
if (readco) {
open(15,file='VELOC.DAT')
write(*,*) 'read old coordinates'
 do i=1,nn {
 read(15,*) xold(3*i-2),xold(3*i-1),xold(3*i)
 }
 do i =1,n3 {
  xold(i)=xold(i)/0.529177d0
 }
 close (15)
}

#
# Setup for relaxation loop
#
icycle = 0; maxrun = 0; irun = 1; converge = .false.
e = 0.0d0; eel = 0.0d0; dxsmx = 0.0d0; fsmx = 0.0d0; 
#
# Start relaxation loop 
# Multiple input lines are only allowed if (mdtemp), they do not make
# sense in any other case. There would be a problem with xold if the
# time stepwidth is changed since xold is with respect to the old time
# stepwidth. However, since vnorm is always called in the first step,
# xold will be set to a reasonable value. 
#
# BFGS optimization in internal coordinates starts here
if (lbfgs) {
#19 continue
#  call bfgs(n3,telec,x,xparam,scftol,atomic,e,eel,gr,niter_
#       ,evector,qmat,fmax,izp,lzmat,lexc,nexc,nrem)
#  converge = .true.
#  icycle=  icycle+1 
#  goto 8
write(*,*) 'this mode is not available'
stop
}
# BFGS optimization in xyz coordinates
221 continue
if (bfgsmin) {
   write (*,*) 'enter deltat, tatom, telec, wvscale, maxrun'
   read (*,*,end = 9) deltat, tatom, telec, wvscale, maxrun
   write (*,*) 'set to: ',deltat, tatom, telec, wvscale, maxrun
   
   if (deltat < 1.0d-10 ) maxrun=0

   ### prepare 
   icycle = icycle+1
   if (icycle==1) write(*,12)
   # calculate atomization energy
   geseatom=0.0d0
   do i=1,nn {
    geseatom = geseatom+espin(izp(i))
   }

   # BFGS calls eglcao(n3,telec,x,scftol,atomic,e,eel,gr,niter,evector,qmat,ql,qlup,qldown,lgrad,lexc,nexc,nrem)
   call bfgsgeom(n3,telec,x,scftol,atomic,e,eel,gr,niter,evector,qmat,_
     lgrad,lexc,nexc,nrem,geseatom,icycle,nbeweg,izp,xm,maxrun,fmax,name,_
     namev,period,ntype,boxsiz,atnames,mode,outfile,qzero,converge,_
     bfgserror,ql,qlup,qldown,lmax,docm3,qcm3)
   if (bfgserror<3) {
     goto 8
   }
   else {
     write(*,*) "Now switch to Conjugate Gradient minimization."
     cgrmin=.true.
     icycle=icycle-1
     maxrunhlp=maxrun
     maxrun=0
   }
}
# for all other, go into this loop
repeat {
 
  if (irun > maxrun) {
   icycle = icycle+1
    irun = 1
   if ((lcc!=0)&&(maxrun!=0) ) {
    rccval = rccval + rccinc
    converge = .false.
    if (cgrmin) call cgrad(3,nbew3,fmax,e,x,gr)
   }
#
# Close cgrad files and leave
#
    if ((irun > 1) && (! mdtemp) && (! mdfree) && (lcc==0)) {   
      icgr = 1
      if (cgrmin) call cgrad(icgr,nbew3,fmax,e,x,gr)
      goto 9
     }
#
#
# Read in deltat, tatom, telec, wvscale, maxrun
# deltat is the time stepwidth in 2.4189d-17 s. tatom is the temperature
# of the atomar system, telec the the temperature of the electronic
# system, both in K. wvscale is the probability to rescale the
# velocities after each step according to the temperature. It may be
# seen as a simple coupling parameter between the system and a
# heat reservoir. maxrun is maximum number of steps to perform
# in that cycle.
#
    dorelax = .true. ; doecalc = .true.
    if (bfgsmin) goto 222 
    write (*,*) 'enter deltat, tatom, telec, wvscale, maxrun'
    read (*,*,end = 9) deltat, tatom, telec, wvscale, maxrun
    write (*,*) 'set to: ',deltat, tatom, telec, wvscale, maxrun
222 continue

    if(lexc) {
      read(*,*) nexc,nrem
    }
    if (mdfree && (icycle !=1)) write(*,*) 'careful: same timestep?'
    if(atomic) {
     deltat = 0.0d0
     maxrun = 0
     }
    delth2 = deltat**2
    if(telec > 5.0 && maxrun>1) {
      write(*,*) 'Warning: Relaxation with finite electron temperature may not'
      write(*,*) '         lead to the right result!'
    }
#
# Check if deltat is greater than zero 
#
    if (deltat < 1.0d-10) {
      dorelax = .false. ;
      if (maxrun > 0) {
        maxrun = 0
        print *,'maxrun has been set to 0'
      }
    }
    else rcdt = 1.0d0/deltat
    if (maxrun < 0) doecalc = .false.
# FOR MD! change if you read velocities or
# do multiple input lines for MD: then do not
# use mdcorrect every icycle
#     if (mdtemp || mdfree) mdcorr = .true. ; else mdcorr = .false.
     if (mdtemp) mdcorr = .true. ; else mdcorr = .false.
     if (mdfree && (icycle==1) && (!readco)) mdcorr = .true. 
        
    if (icycle == 1 ) {
      if (bfgsmin) open(95,file='ENERGY.TMP',form='formatted',status='unknown',access='append')
      else         open(95,file='ENERGY.TMP',form='formatted',status='unknown')
      if ( mode==1 || mode==2) {
        write (*,11)
      }
      else {
        write (*,12)
      }
    }
  }
  if (bfgserror>2) {
    bfgsmin=.false.
    maxrun=maxrunhlp
  }
#
# Done with setup
#
  if (dorelax && (! cgrmin)) {
#
# Mass centre --> 0,0,0
#  
#   if (! constr ){
    if (EXT=='NO') {
      if (deltat > 1.0d-10) call xnorm(nn,x,xold,xm,izp)
    } 
#   }
#
# All atoms back into the supercell
#
    if (period) {
      
      do i = 1,nn {
      call coordback(x(3*i-2),x(3*i-1),x(3*i),xold(3*i-2),xold(3*i-1),xold(3*i))
      }
    }
    
    
  }
#
# Set total impulse = 0, scale velocities according to temperature
# Velocities must always be scaled in the first step
#
  if (dorelax && mdtemp) {
    if (irun == 1) write(*,*) 'rescale velocities'
    if (irun == 1) call vnorm(nbeweg,x,xold,tatom,deltat,xm,izp)
#    if ((irun == 1) && (icycle == 1))  _
#     call vnorm(nbeweg,x,xold,tatom,deltat,xm,izp)
    else if (ranmar() < wvscale) {
    write(*,*) 'ranmar() < wvscale'
     mdcorr = .true.
     call vnorm(nbeweg,x,xold,tatom,deltat,xm,izp) 
    }
  }
#
# If (! doecalc), skip energy evaluation and geometry update
#
  if (! doecalc) {goto 7}

#
# Do the work (total energy, gradient)
#
  call eglcao(n3,telec,x,scftol,atomic,e,eel,gr,niter,_
              evector,qmat,ql,qlup,qldown,lgrad,lexc,nexc,nrem,docm3,qcm3,xm,excited,docpe)
 

### numerical forces for debugging:                                             
#  e0      = e
#  deltagr = 0.00001d0
#  do i=1,n3{
#    xhlp(i)=x(i)
#  }
#  do i=1,n3{
#    lcount=0
#    do k = 1,nn {
#      qmat(k) = qzero(izp(k),4)
#      do j=1,lmax(izp(k)){
#        lcount=lcount+1
#        ql(lcount)     = qzero(izp(k),j)
#        qlup(lcount)   = qzero(izp(k),j)/2.0d0
#        qldown(lcount) = qzero(izp(k),j)/2.0d0
#      }
#    }
#    xhlp(i)=x(i)-deltagr
#    call eglcao(n3,telec,xhlp,scftol,atomic,em,eel,gr,niter,evector,qmat,_
#       ql,qlup,qldown,lgrad,lexc,nexc,nrem,qcm3,xm)
#    xhlp(i)=x(i)+deltagr
#    call eglcao(n3,telec,xhlp,scftol,atomic,ep,eel,gr,niter,evector,qmat,_
#       ql,qlup,qldown,lgrad,lexc,nexc,nrem,qcm3,xm)
#    gr(i)=(ep-em)/(2.0d0*deltagr)
#    xhlp(i)=x(i)
#  }
### end numerical forces for debugging

### print out trajectory - for easier debugging
#     namev=name
#     write (name,'(a,i3.3)') outfile(:index(outfile,' ')-1), irun
#     call outcoord(name, namev, nn, period, ntype, izp, _
#       x, xold, boxsiz, atnames,mode)

#
# Output of necessary imnformation (forces, masses, etc.) for vibrational modes
#
#
# if (atomic) {
 call outforces(nbeweg,izp,gr,xm)
# }

#
# Calculate fsmx, the maximal force acting on a single atom
#
  fsmx = 0.0
  do i = 1,nbeweg {
    ftmp = 0.0
    do j = -2,0 {
      ftmp = ftmp + gr(3*i+j)**2
    }
    ftmp = sqrt(ftmp)
    fsmx = max(fsmx,ftmp)
  }
#
# Check if conjugate gradient or stdc method is converged
# Setup icgr
# If (converge && cgrmin), delete file cgrad 
#
  if (stdc || cgrmin) {
    if (fsmx < fmax) converge = .true.
  }
  icgr = 0
  if (converge) icgr = 2
  if (converge && cgrmin) call cgrad(icgr,nbew3,fmax,e,x,gr)
#
# If (converge || (! dorelax)), skip geometry update
#
  if (converge || (! dorelax)) goto 7
#
# Geometry update (conjugate gradient, stdc or verlet)
# Hint: When using the Verlet algorithm and velocities at the same 
# time, we have a problem. The time dependence of x can be written as
#
# (1)  x(t+dt) = x(t) + v(t)*dt + 0.5*dt*dt*f(t)/m
#
# When we use Verlet, we approximate v(t) pretty reasonable as
#
# (2)  v(t) = (x(t+dt)-x(t-dt))/(2*dt)
# 
# Substituting (2) in (1) leads to the Verlet formula
# 
# (3)  x(t+dt) = 2*x(t) - x(t-dt) + dt*dt*f(t)/m
#
# This geometry update is due to v(t) and f(t). That means, we have the
# coordinates for x(t+dt), but not the velocities. But we desperately 
# need these velocities to determine the temperature of the system.
# Thus we have to make a compromise: use v(t+dt) = (x(t+dt) - x(t))/dt 
# to determine the temperature and to rescale the velocities.
# There is one more thing to take care of: There are a lot of cases when
# we exactly know v(t): either in the very first step of the MD simulation
# or after we did a rescaling of the velocities according to the current
# temperature. In either case, we have determined the velocity of the 
# atoms using the formula v(t) = (x(t)-x(t-dt))/dt, which is in
# contrast to (2). However, after calling eglcao we know f(t) and
# are able to correct for this error by resetting x(t-dt) according to
#
# (4)  x(t-dt) = x(t-dt) + 0.5*f(t)/m
#
# This is what is done in the next loop. It is executed if (irun == 1) 
# or if vscale has been called.
# Since the gradient of the frozen coordinates has already been set to
# zero, they will not be changed.
#
  if (mdcorr  ) {
    mdcorr = .false.
     write(*,*) 'mdcorrect: correct velocities'
    do i = 1,nbew3 {
     xold(i) = xold(i) + 0.5d0*delth2*gr(i)*xmrc(izp((i+2)/3))
    }
  }
#
# Now start with the actual geometry update
# dxsmx is the maximal displacement of a single atom
#
  ekin = 0.0d0
  dxsmx = 0.0d0
  if (cgrmin || stdc) {
    if (cgrmin) call cgrad(icgr,nbew3,fmax,e,x,gr)
    else {
      do i = 1,nbew3 {
      x(i) = xold(i) - delth2*gr(i)*xmrc(izp((i+2)/3))
      }
    }
    do i = 1,nbeweg {
      dxtmp = 0.0d0
      do j = -2,0 {
        dxtmp = dxtmp + (x(3*i+j)-xold(3*i+j))**2
        xold(3*i+j) = x(3*i+j)
      }
    dxtmp = sqrt(dxtmp)
    dxsmx = max(dxsmx,dxtmp)
    }
  }
  else if (mdtemp || mdfree) { 
    do i = 1,nbew3 { 
      xv = 2*x(i) - xold(i) - delth2*gr(i)*xmrc(izp((i+2)/3)) 
      ekin = ekin + 0.5d0*(xv-xold(i))**2/ _
        (4.0d0*delth2*xmrc(izp((i+2)/3)) )
      xold(i) = x(i)
      x(i) = xv
    }
    do i = 1,nbeweg {
      dxtmp = 0.0d0
      do j = -2,0 {
        dxtmp = dxtmp + (x(3*i+j)-xold(3*i+j))**2
      }
      dxtmp = sqrt(dxtmp)
      dxsmx = max(dxsmx,dxtmp)
    }
  }
#
# create neighbour table, get number of atoms without
# neighbours
#
7 nfree = 0
  
#
# normalize coordinates if one cycle is complete
#
  if (converge || (irun >= maxrun)) {
    if (dorelax && (!cgrmin) ) {
      if (EXT=='NO') {
        call xnorm(nn,x,xold,xm,izp)
      }
    }
  }
# calculate atomization energy
  geseatom=0.0d0
  do i=1,nn {
    geseatom = geseatom+espin(izp(i))
  }
#
#
# output to stdout
#
  if ( mode==1 || mode ==2) {
  write (*,13) icycle,irun,niter,e,eel,(e-geseatom)*627.5095d0,Edis,_
    ekin,ekin+e
  write (95,13) icycle,irun,niter,e,eel,(e-geseatom)*627.5095d0,Edis,_
    ekin,ekin+e
  }
  else  {
  write (*,13) icycle,irun,niter,e,eel,(e-geseatom)*627.5095d0,Edis,dxsmx, _
    fsmx
  write (95,13) icycle,irun,niter,e,eel,(e-geseatom)*627.5095d0,Edis,dxsmx, _
    fsmx
 }
if ( (converge||(irun==maxrun)) && (lcc!=0) ) {
    write(*,'(A10,2I5,6F12.4)') ' constr',icycle,irun,  _
    rccval*0.529177d0,e,e*627.5095d0,(e-geseatom)*627.5095d0,(e-eel)*627.5095d0,fsmx 
    write(95,'(A10,2I5,6F12.4)') ' constr',icycle,irun,  _
    rccval*0.529177d0,e,e*627.5095d0,(e-geseatom)*627.5095d0,(e-eel)*627.5095d0,fsmx 
    irun = maxrun+1
   }
#  close(95)
#
# output to structure outputfile
#
  irun = irun+1
#8 continue  
  write (name,'(a,i3.3)') outfile(:index(outfile,' ')-1), icycle
  if ( mode==1 || mode ==2 ) {
  write (namev,'(a,i3.3,a)') outfile(:index(outfile,' ')-1), icycle,'v'
  }
#  if ( mode == 3 ) {
#  write (name,'(a,i3.3)') outfile(:index(outfile,' ')-1), irun
#  }
#  if(lzmat) {
#   call zoutput(name,izp,xparam)
#  }
#  else {
   call outcoord(name, namev, nn, period, ntype, izp, _
     x, xold, boxsiz, atnames,mode)
#  }
8 continue 
# K.W. 2013-07-11
  if (converge && (lcc==0) ) goto 9
#  if (lbfgs && (lcc!=0)) goto 19
  if (bfgsmin) goto 221
#
} 
# end irun loop
#
# end of relaxation loop
#
9 continue

### calculate single point energies for all strucs where one atom is 
### moved by +or-0.001 AA in one direction 
### finally output files for phonir-programm
if (bcalcvib){
  open (1,file='POSITION',status='replace'); close(1)
  open (2,file='FORCE',status='replace')   ; close(2)
  open (3,file='DIPOLE',status='replace')  ; close(3)
  open (4,file='ENERGY',status='replace')  ; close(4)
# K.W. 2013-06-20 
  open (5,file='CM3DIPOLE',status='replace'); close(5)

  call outvib(.true.,nn,atnames,x,gr,xm,izp,e,qzero,qmat)
  # at1 x y z, at2 x y z,...
  do i=1,n3{ 
    xhlp(i)=x(i)
  }
  do i=1,n3{ 
    lcount=0
    do k = 1,nn {
      qmat(k) = qzero(izp(k),4)
      do j=1,lmax(izp(k)){
        lcount=lcount+1
        ql(lcount)     = qzero(izp(k),j)
        qlup(lcount)   = qzero(izp(k),j)/2.0d0
        qldown(lcount) = qzero(izp(k),j)/2.0d0
      }
    }
    do j=1,2 {
      xhlp(i)=xhlp(i)+0.001d0/0.529177d0*((j-2)*2+1)
      call eglcao(n3,telec,xhlp,scftol,atomic,e,eel,gr,niter,_
                  evector,qmat,ql,qlup,qldown,lgrad,lexc,nexc,nrem,docm3,qcm3,xm)
      call outvib(.false.,nn,atnames,xhlp,gr,xm,izp,e,qzero,qmat,docm3,qcm3)
      xhlp(i)=x(i)
    }
  }
  # restore all files that are written out by eglcao before bcalcvib
  lcount=0
  do k = 1,nn {
    qmat(k) = qzero(izp(k),4)
    do j=1,lmax(izp(k)){
      lcount=lcount+1
      ql(lcount)     = qzero(izp(k),j)
      qlup(lcount)   = qzero(izp(k),j)/2.0d0
      qldown(lcount) = qzero(izp(k),j)/2.0d0
    }
  }
  call eglcao(n3,telec,x,scftol,atomic,e,eel,gr,niter,_
              evector,qmat,ql,qlup,qldown,lgrad,lexc,nexc,nrem,qmat,docm3,qcm3,xm)

}


  
#
# format specifications
#
11 format(/2x,'icycle iter niter',6x,'e(total)',5x,'e(bandstr)', _
           5x,'e(bind[kcal])',4x,'Edis',8x,_
            'deltax(max)',6x,'f(max)',8x,_
            'ekin',9x,'etot'/,1x,105('='))
12 format(/2x,'icycle iter niter',6x,'e(total)',5x,'e(bandstr)',_
           5x,'e(bind[kcal])',4x,'Edis',8x,'deltax(max)',6x,'f(max)'/,_
            1x,105('='))
#12 format(/2x,'icycle iter niter',6x,'e(total)',5x,'e(bandstr)',_
#                5x,'e(bind[kcal])',4x,'Edis',8x,'deltax(max)',6x,'f(max)'/,_
#		                 1x,135('='))
#            'ekin',9x,'etot'/,1x,135('='))
#13 format(i5,1x,i5,' / ',i2,3(1x,f14.6),1x,f14.6,4(1x,f14.6))
#13  format(i5,1x,i5,' / ',i2,3(1x,f14.6),1x,f14.6,2(1x,f14.6))
13  format(i5,1x,i5,' / ',i2,3(1x,f14.6),1x,f14.6,2(1x,f14.6))

call cpu_time(cputime_stop)
call system_clock ( clock_count_stop, clock_rate, clock_max )
print *,' '
print *,'CPU time in seconds:   ',cputime_stop - cputime_start
print *,'real time in seconds: ',(clock_count_stop - clock_count_start)/real(clock_rate)
print *,' '
print *,'***** end of dftb *****'
end

