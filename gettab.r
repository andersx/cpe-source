subroutine gettab(ntype)
implicit REAL*8 (A-H,O-Z)
include 'maxima.h'
real*8 skhtab(10,MAXTAB,MAXTYP,MAXTYP),skstab(10,MAXTAB,MAXTYP,MAXTYP)
real*8 espin(MAXTYP),wspin(MAXTYP,3,3)
real*8 skself(3,MAXTYP),dr(MAXTYP,MAXTYP)
real*8 coeff(6,MAXINT,MAXTYP,MAXTYP),xr(2,MAXINT,MAXTYP,MAXTYP)
real*8 efkt(3,MAXTYP,MAXTYP),cutoff(MAXTYP,MAXTYP)
integer numint(MAXTYP,MAXTYP)
integer dim(MAXTYP,MAXTYP),ppp,lmax(MAXTYP),nelup,neldown
character*128 skfile
character*128 chdummy
logical spin,ldep
common /sktab/ skhtab,skstab,skself,dr,dim
common /spin/ espin,nelup,neldown,wspin,spin
common /spltab/ coeff,xr,efkt,cutoff,numint
common /mcharge/ qzero(MAXTYP,4), uhubb(MAXTYP,3), ldep
common /lmax/ lmax

do i = 1,ntype {
 qzero(i,4) = 0.0d0
 do k=1,3 {
   qzero(i,k) = 0.0d0
   uhubb(i,k) = 0.0d0
 }
  do j = 1,ntype {
    print *,'enter file name for sk-data of pair ',i,j
    read *,skfile
    open (3,file=skfile,status='unknown')
    rewind 3
    if (i != j) {
    read (3,*) dr(i,j),dim(i,j)
    }
    if (i == j) { 
#     read (3,*) dr(i,j),dim(i,j),lmax(i)
      read (3,*) dr(i,j),dim(i,j)
      read (3,*) (skself(l,i),l = 1,3),espin(i), _
                 (uhubb(i,4-l), l=1,3), (qzero(i,4-l), l=1,3)
      do k = 1,3 {
       qzero(i,4)   = qzero(i,4) + qzero(i,k)
      }
      do l=1,3 {
       espin(i) = espin(i) + skself(l,i)*qzero(i,4-l)
      }
    }
    do k = 1,dim(i,j) {
      read (3,*) (skhtab(l,k,i,j),l = 1,10),_
                 (skstab(l,k,i,j),l = 1,10)
    }
    while (.true.) {
     read(3,'(A)') chdummy
     if (chdummy=='Spline') break
    }   
    read(3,*) numint(i,j),cutoff(i,j)
    if(numint(i,j)>MAXINT) {
     print *,'Too many intervalls!'
     goto 99
    }
    read(3,*) (efkt(ppp,i,j),ppp=1,3)

    do ppp=1,numint(i,j) {
      if(ppp<numint(i,j)) {
       read (3,*) xr(1,ppp,i,j),xr(2,ppp,i,j),coeff(1,ppp,i,j),coeff(2,ppp,i,j),coeff(3,ppp,i,j),coeff(4,ppp,i,j)
      }
      else {
       read (3,*) xr(1,ppp,i,j),xr(2,ppp,i,j),coeff(1,ppp,i,j),coeff(2,ppp,i,j),coeff(3,ppp,i,j),coeff(4,ppp,i,j),coeff(5,ppp,i,j),coeff(6,ppp,i,j)
      }
    }

    if(xr(2,numint(i,j),i,j)!=cutoff(i,j)) {
     print *,'Error in Datafile'
     goto 99
    }
    print *,'skfile for pair ',i,j,' :',skfile
    close (3)
  }
}
99 continue
end

