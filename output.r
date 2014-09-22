#     SUBROUTINE OUT_EIGENVECTORS
#     =================
#
#     Copyright 1997 by Adam Gali
#
# *********************************************************************
#
#     PROGRAM CHARACTERISTICS
#     -----------------------
#
# print eigenvectors prints out the eigenvalues and the contributed    
# eigenvectors in ascending order in file 'EVC.DAT'. 
# No effects to the calculation!!
#
# PARAMETERS:
# a      r() the matrix of eigenvectors
# ev     r() the matrix of eigenvalues
# occ    r() vector of occupation number
# ind    i() index of orbitals possessed by one atom
# nn       i number of atoms
#
# *********************************************************************
#
subroutine  outeigenvectors(a,ev,occ,ind,nn)
#################################
#
# maxima definitions
#
implicit REAL*8 (A-H,O-Z)
include 'maxima.h'
#
integer ind(NNDIM+1), izp(NNDIM), nn
integer i,j,k,l
integer neig,filled,hfilled
#  number of eigenstates, filled levels and half occupied states

real*8 a(MDIM,MDIM), ev(MDIM), occ(MDIM), nel
character*4 orbital(9)

common /izp/ izp, nel, nbeweg
#              atom types, number of electrons, number of movable atoms

data orbital/'S   ',
             'Px  ',
             'Py  ',
             'Pz  ',
             'Dxx ',
             'Dxy ',
             'Dyy ',
             'Dyz ',
             'Dzz '/


#The last index points to the number of last eigenstates
#

neig = ind(nn+1)

filled = 0
hfilled = 0
do i = 1,neig {
  if(occ(i) > 1) {
    filled = filled+1
  }
  else {
    if(occ(i) > 0) hfilled = hfilled+1
  }
}

open (1,file='EVC.DAT',status='unknown'); rewind 1
write (1,'(20x,a)') 'THE ONE-ELECTRON EIGENVALUES AND EIGENVECTORS'
write (1,*)
write (1,*)
write (1,*) 'NUMBER OF FILLED LEVELS: ',filled
if(hfilled > 0) write (1,*) 'NUMBER OF HALF OCCUPIED STATES: ',hfilled
write (1,*)
write (1,*)

do i=1,neig{
  write(1,13) 'THE ',i,'. EIGENVALUE:',ev(i),' H',ev(i)*27.2116,' eV'
  write(1,*) 'Atom No.   Atom Type'
  do j=1,nn{
    write(1,*) '  ',j,'         ',izp(j)
    k=0
    do l=ind(j)+1, ind(j+1){
      k=k+1
      write(1,'(a,f8.5)') orbital(k),a(l,i)
    }
 }
 write(1,*)
} 
  
endfile 1; close (1)

# format specification
#
13 format(1x,a,i4,a,f15.8,a,4x,f12.5,a)

END



subroutine outspec(nn,nmov,ndim,ind,lmax,izp,ev,occ,efermi, _ 
               qmat,qmulli,qtot,dipol,dipabs, _
               spin,evdown,occdown,qmatup,qmatdown,qmulliup,_ 
               qmullidown,qtotup,qtotdown,efermiup,efermidown,_
               qzero,docm3,qcm3,qcm3tot,cm3dipol,cm3dipabs,cbnd,_ 
               order,number1,number2,bord,_
               docpe,E_cpe,C,cpedipabs,cpedipol,grad_cpe,shiftCPE)

# maxima definitions
 
  implicit REAL*8 (A-H,O-Z)
  include 'maxima.h'
#
  integer nn,ndim,ind(NNDIM+1),lmax(MAXTYP),izp(NNDIM)
  real*8 ev(MDIM),occ(MDIM),efermi,qmat(NNDIM),qmulli(MDIM)
  real*8 qtot,dipol(3),dipabs
  logical spin
  real*8 evdown(MDIM),occdown(MDIM),qtotup,qtotdown
  real*8 qmulliup(MDIM),qmullidown(MDIM),qmatup(NNDIM),qmatdown(NNDIM)
  integer i
# K.W. 2013-06-20 additional stuff for CM3 charge and bond order output
  logical docm3
  real*8 qzero(MAXTYP,4)
  real*8 qcm3(NNDIM),cm3dipol(3),cm3dipabs, qcm3tot
  real*8 order(NNDIM),bord(NNDIM,NNDIM)
  integer number1(NNDIM),number2(NNDIM),cbnd
# K.W. 2014-05-01 additional stuff for CPE calculation 
  logical docpe
  real*8 E_cpe,C(3*NNDIM),cpedipabs,cpedipol(3), shiftCPE(NNDIM), grad_cpe(3,NNDIM)

  open (1,file='SPE.DAT',status='unknown'); rewind 1
  if (spin) {
    write(1,*) 'spin up'
    write (1,'(f20.12)') efermiup
  } 
  else write (1,'(f20.12)') efermi
  write (1,'(2(1x,f20.12))') (ev(i),occ(i),i = 1,ndim)
  if (spin) {
    write(1,*) 'spin down'
    write (1,'(f20.12)') efermidown
    write (1,'(2(1x,f20.12))') (evdown(i),occdown(i),i = 1,ndim)
  }
  endfile 1; close (1)

  open (1,file='CHR.DAT',status='unknown'); rewind 1
  write (1,'(a)') '#MULLIKEN CHARGE'
  write (1,'(a)') '#CHARGE PER ATOM / PER ATOMIC STATES'
  write (1,'(a)') '#ATOM No.'
  write (1,'(a,2x,i4)') '#',nmov
  write (1,'(a,2x,i4)') '#',nn
  if (spin) {
    write(1,*) 'spin up'
    do n = 1,nn {
      ind1 = ind(n)+1; ind2 = ind1+lmax(izp(n))**2-1
      write (1,98) n,qmatup(n),(qmulliup(i), i= ind1,ind2)
    }
    write(1,*) 'spin down'
    do n = 1,nn {
      ind1 = ind(n)+1; ind2 = ind1+lmax(izp(n))**2-1
      write (1,98) n,qmatdown(n),(qmullidown(i), i= ind1,ind2)
    }
  }
  else{
    do n = 1,nn {
      ind1 = ind(n)+1; ind2 = ind1+lmax(izp(n))**2-1
      write (1,98) n,qmat(n),(qmulli(i), i= ind1,ind2)
    }
  }
  endfile 1; close (1)
  
98  format(1x,i4,11(1x,f12.6))
# K.W. 2013-06-20 additional file 'CHARGES.DAT'
  if (docm3) {
    open (1,file='CHARGES.DAT',status='unknown'); rewind 1
    write (1,'(a)') '# ATOMNR. ATOM-TYPE  MULLIKEN   CM3'
    do i = 1,nn {
      write (1,99) i,izp(i), qzero(izp(i),4)-qmat(i), qzero(izp(i),4)-qcm3(i)
    }
    endfile 1; close(1)
  }
99  format(1x,i4,1x,i4,2f12.6)

# K.W. 2014-05-01 additional file 'CPE.DAT'
  if (docpe) {
    open (1,file='CPE.DAT', status='unknown'); rewind 1
    write (1,'(a)') 'Output from CPE correction'
    write (1,'(a)') ' '
    write (1,'(a)') 'Coefficients for CPE Slater functions'
    write (1,'(a)') '====================================='
    do i = 1, 3*nn {
      write (1,'(f12.6)') C(i)
    }
    write (1,'(a)') ' '
    write (1,'(a)') 'CPE contribution to the Hamiltonian'
    write (1,'(a)') '====================================='
    do i = 1, nn {
      write (1,'(f12.6)') shiftCPE(i)
    }
    write (1,'(a)') ' '
    write (1,'(a)') 'CPE contribution to total energy'
    write (1,'(a)') '================================'
    write (1,'((a),f9.6)') 'E_CPE: ', E_cpe
    write (1,'(a)') ' '
    write (1,'(a)') 'CPE contribution to gradient'
    write (1,'(a)') '================================'
    do j = 1, nn {
      write (1,'(i5,3(1x,e18.6))') j,(grad_cpe(i,j), i = 1,3)
    }
    write (1,'(a)') ' '
    write (1,'(a)') 'CPE corrected dipole moment'
    write (1,'(a)') '==========================='
    write (1, '(1x,3(f12.4))') (cpedipol(i),i=1,3)
    write (1,'(a)') ' '
    write (1,'(a)') '                       norm'
    write (1,'(a)') '================================'
    write (1,'(1x,f9.6)') cpedipabs
    endfile 1; close(1)
  }
  
  open (1,file='REST.DAT',status='unknown'); rewind 1
  if (spin) {
    write (1,'(a)') 'Total charge, total spin-up charge, total spin-down charge'
    write (1,'(a)') '============'
    write (1,'(3(1x,f20.12))') qtot, qtotup, qtotdown
  }
  else {
    write (1,'(a)') 'Total charge'
    write (1,'(a)') '============'
    write (1,'(1x,f20.12)') qtot
  }
  write (1,'(a)') ' '
  write (1,'(a)') 'Mulliken dipole moment [Debye]'
  write (1,'(a)') '=============================='
  write (1,'(3(1x,f20.12))') (dipol(i), i= 1,3)
  write (1,'(a)') ' '
  write (1,'(a)') 'Norm of dipole moment [Debye]'
  write (1,'(a)') '============================='
  write (1,'(1x,f20.12)') dipabs
  write (1,'(a)') 
# K.W. 2013-06-20 additional content of CM3 dipole moment
  if (docm3) {
    write (1,'(a)') ' '
    write (1,'(a)') 'Total charge'
    write (1,'(a)') '============'
    write (1,'(1x,f20.12)') qcm3tot
    write (1,'(a)') ' '
    write (1,'(a)') 'CM3 dipole moment [Debye]'
    write (1,'(a)') '=============================='
    write (1,'(3(1x,f20.12))') (cm3dipol(i), i= 1,3)
    write (1,'(a)') ' '
    write (1,'(a)') 'Norm of CM3 dipole moment [Debye]'
    write (1,'(a)') '============================='
    write (1,'(1x,f20.12)') cm3dipabs
    write (1,'(a)') 
  }
# K.W. 2014-05-01 additional contant of CPE dipole moment
  if (docpe) {
    write (1,'(a)') 'CPE dipole moment [Debye]'
    write (1,'(a)') '=============================='
    write (1,'(3(1x,f20.12))') (cpedipol(i), i= 1,3)
    write (1,'(a)') ' '
    write (1,'(a)') 'Norm of CPE dipole moment [Debye]'
    write (1,'(a)') '============================='
    write (1,'(1x,f20.12)') cpedipabs
    write (1,'(a)') 
  }
  
  endfile 1; close (1)
# K.W. 2013-06-20 additional file for bond order output
  if (docm3) {
    open (1, file='BORD.DAT', status='unknown'); rewind 1
    write (1, '(a,3x,a,3x,a)') 'atom1','atom2','bond order'
    write (1, '(a)') ' '
    do j = 1,cbnd {
      write (1,'(2(i4,4x),7x,f5.3)') number1(j),number2(j),order(j)
    }
    write (1,'(a)') ' '
    write (1,'(a,3x,a,3x,a)') 'Atom', 'Bond order matrix'
    write (1,'(a)') ' '
    do j = 1,nn {
      write (1,'(i4,1x,100(1x,f6.3))') j, (bord(i,j),i=1,nn)
    }
    endfile 1; close(1)
  }

END

#
# Output of coordinates , velocities and locked atoms
#
  SUBROUTINE outcoord(name, namev, nn, period, ntype, izp, _
    x, xold, boxsiz, atnames,mode)
#
# maxima definitions
#
   implicit REAL*8 (A-H,O-Z)
   include 'maxima.h'
#
  integer mode
  character*65 name,namev
  character*70 atnames
  character tmpchr
  logical period
  integer nn, ntype, izp(NNDIM)
  real*8 x(3,NNDIM),boxsiz(3,3),xold(3*NNDIM)
  integer i,j
  real*8 xhlp(3)
#      ME: for constraints like in CHARMM
      integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
      real*8 rccval,kccval,rccinc
#      common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc
#

#   
  open (1,file=name,status='unknown'); rewind 1
  n3 = 3*nn
  if ( mode ==1 || mode ==2 ) {
  open(10,file=namev,status='unknown'); rewind 10
   do i=1,nn {
   write(10,'(3F12.6)') xold(3*i-2)*0.529177d0,_
     xold(3*i-1)*0.529177,xold(3*i)*0.529177d0
 
   }
   close (10)
  }
  if(period) {
   tmpchr='S' 
   if(lcc!=0) {
   tmpchr='SM'
   }
   }
  else {
   tmpchr='C'
#  if(lcc!=0) {
#   tmpchr='M'
#  }
  }
  
  write(1,'(I4,1X,A1)') nn,tmpchr
  write(1,'(A)') atnames
  do i=1,nn {
    do j=1,3 {
     xhlp(j)=x(j,i)*0.529177d0
    }
    write(1,14) i,izp(i),(xhlp(j),j=1,3)
  }
  
  if(period) {
   write(1,15) 0.0,0.0,0.0
   do i=1,3 {
     do j=1,3 {
      xhlp(j)=boxsiz(i,j)*0.529177d0
     }
     write(1,15)(xhlp(j),j=1,3)
   }
  }
#     ME: constraints like in CHARMM
#      if (lcc!=0) {
#    write(1,'(I3,3F12.6)') lcc,kccval/0.529177,rccval*0.529177,rccinc*0.529177
#    write(1,'(100I3)')(jcc(i),atcc1(i),atcc2(i),i=1,lcc)
#      }

     
14  format(I4,1X,I3,3(1X,F12.6))
15  format(3(1X,F12.6))  

  endfile 1; close (1)
  
END

SUBROUTINE outforces(nbeweg,izp,gr,xm)
implicit REAL*8 (A-H,O-Z)
integer nbeweg,izp(*)
real*8 gr(*),xm(*)

integer i,j
#
# Output of energy, old coordinates and forces for vibrational modes
#
  open (2,file='FRC.DAT',status='unknown'); rewind 2 
  write(2,'(A)') '# INTERATOMIC FORCES AT ATOM N0. I'
  write(2,'(A)') '# ATOM N0.'
  write(2,'(A)') '# FORCE_X, FORCE_Y, FORCE_Z'
  write(2,'(A)') '#'
  write(2,'(A)') '#'
  do i = 1,nbeweg {
    write (2,'(i5,3(1x,e19.10))') i,(-gr(3*i+j-3), j = 1,3)
  }
  write(2,'(A)') 'rms/d:   0.0000   0.0000'
  close (2)
#
#
# Output of number of atoms and masses for vibrational modes
#
  open (2,file='MAS.DAT',status='unknown'); rewind 2
  write (2,'(1x,i5)') nbeweg
  do i = 1,nbeweg {
    write (2,'(1x,f24.15)') xm(izp(i))
  }
  close (2)
END




SUBROUTINE outeglcao(nn,nmov,ndim,ind,lmax,izp,ev,occ,efermi, _ 
                    qmat,qmulli,qtot)
#
# maxima definitions
#
   implicit REAL*8 (A-H,O-Z)
   include 'maxima.h'
#
  integer nn,ndim,ind(NNDIM+1),lmax(MAXTYP),izp(NNDIM)
  real*8 ev(MDIM),occ(MDIM),efermi,qmat(NNDIM),qmulli(MDIM)
  real*8 qtot

  integer i

  open (1,file='eglcao.out',status='unknown'); rewind 1
  write (1,'(3(1x,f20.12))') (ev(i),occ(i),occ(i),i = 1,ndim)
  write (1,'(f20.12)') efermi
  write (1,'(a)') ' '
  
  write (1,'(a)') 'Mulliken charges'
  write (1,'(a)') '================'
  do n = 1,nn {
      ind1 = ind(n)+1; ind2 = ind1+lmax(izp(n))**2-1
      write (1,72) 'Atom ',n,': ',qmat(n)
      write (1,73) (qmulli(i), i= ind1,ind2)
  }
  write (1,'(a)') ' '
  write (1,'(a)') 'Total charge'
  write (1,'(a)') '============'
  write (1,'(1x,f20.12)') qtot
  write (1,'(a)') ' '
  endfile 1; close (1)
72  format(a,i4,a,f20.12)
73  format(3(1x,f20.12))
END
