#
# Calculates Hamilton and Overlap matrix
# for Atoms i and j, the results are stored in ham,over
#
SUBROUTINE slkmatrices(i,j,xat,ham,over)
implicit none
include 'maxima.h'
integer i,j, izp(NNDIM), nbeweg, l1,l2,nu,nv,nw, _
        ind1,ind2,ind3
real*8 xat(3,*), nel
real*8 ham(LDIM,LDIM),over(LDIM,LDIM),ham_tmp(LDIM,LDIM),over_tmp(LDIM,LDIM)
include 'commonbox.h'
include 'phihelp.inc'
common /izp/ izp, nel, nbeweg
external skspar,skhpar

integer izpi,izpj
real*8 dif(3),dif0(3)

dif0(1)=xat(1,j)-xat(1,i)
dif0(2)=xat(2,j)-xat(2,i)
dif0(3)=xat(3,j)-xat(3,i)
izpi=izp(i)
izpj=izp(j)
DO l2 = 1,LDIM {
 DO l1 = 1,LDIM {
   ham(l1,l2)  = 0.0d0
   over(l1,l2) = 0.0d0
 }
}

if(period) {
 
 DO nu =-nlat(1),nlat(1) {
  ind3=nu+nmax+1
  DO nv = -nlat(2),nlat(2) {
   ind2=nv+nmax+1
   DO nw = -nlat(3),nlat(3) {
    ind1=nw+nmax+1
    dif(1) = dif0(1)-sumlat(1,ind1,ind2,ind3)
    dif(2) = dif0(2)-sumlat(2,ind1,ind2,ind3)
    dif(3) = dif0(3)-sumlat(3,ind1,ind2,ind3)
    call SLKODE(dif,skcut2,izpi,izpj,ham_tmp,skhpar)
    call SLKODE(dif,skcut2,izpi,izpj,over_tmp,skspar)
    if (izpi == -1969 ) {
      izpi=izp(i)
    } else {
      DO l2 = 1,LDIM {
       DO l1 = 1,LDIM {
        ham(l1,l2)  = ham(l1,l2)  + ham_tmp(l1,l2)
        over(l1,l2) = over(l1,l2) + over_tmp(l1,l2)
       }
      }
    }
   
   }
  }
 }

}
else {
 call SLKODE(dif0,skcut2,izpi,izpj,ham,skhpar)
 call SLKODE(dif0,skcut2,izpi,izpj,over,skspar)
}

END


SUBROUTINE SLKODE(DUM,skcut2,I,J,EM,iovpar)
implicit REAL*8 (A-H,O-Z)
include 'maxima.h'
real*8 DUM(3),EM(LDIM,LDIM)
real*8 X(6),X2(6),dummy(LDIM,LDIM)
external iovpar
common /lmax/ lmax(MAXTYP)
R2=0.0
DO L=1,3 {
 X(L)=DUM(L)
 X2(L)=X(L)*X(L)
 R2=R2+X2(L)
}
#!!!!!!!!!!!!!!!!!!!!!!

if (R2>=skcut2) {
  I=-1969
  RETURN
}
#!!!!!!!!!!!!!!!!!!!!!
IF(R2 >= 1.0E-8) {
 R2I=1.0/R2
 RI=SQRT(R2I)
 DO L=1,3 {
  X(L)=X(L)*RI
  X(L+3)=X(L)
  X2(L)=X2(L)*R2I
  X2(L+3)=X2(L)
 }
 maxmax = max(lmax(i),lmax(j))
 minmax = min(lmax(i),lmax(j))
#
# s Interaction
#
 call skss(x,x2,i,j,r2,iovpar,em,LDIM)
 if (maxmax <= 1) return
#
# p Interaction
#
 if (minmax >= 2) {
  call skpp(x,x2,i,j,r2,iovpar,em(2,2),LDIM)
  call sksp(x,x2,i,j,r2,iovpar,em(1,2),em(2,1),LDIM)
  if (i != j)
   call sksp(x,x2,j,i,r2,iovpar,dummy,em(2,1),LDIM)
 }
 else if (lmax(j) >= 2)
  call sksp(x,x2,i,j,r2,iovpar,em(1,2),em(2,1),LDIM)
 else
  call sksp(x,x2,j,i,r2,iovpar,dummy,em(2,1),LDIM)
 if (maxmax <= 2) return
#
# d Interaction
#
 if (minmax == 3) {
  call skdd(x,x2,i,j,r2,iovpar,em(5,5),LDIM)
  call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),LDIM)
  call skpd(x,x2,i,j,r2,iovpar,em(2,5),em(5,2),LDIM)
  if (i != j) {
   call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),LDIM)
   call skpd(x,x2,j,i,r2,iovpar,dummy,em(5,2),LDIM)
  }
 }
 else if (lmax(i) == 1)
  call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),LDIM)
 else if (lmax(i) == 2) {
  call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),LDIM)
  call skpd(x,x2,i,j,r2,iovpar,em(2,5),em(5,2),LDIM)
 }
 else if (lmax(j) == 1)
  call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),LDIM)
 else {
  call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),LDIM)
  call skpd(x,x2,j,i,r2,iovpar,dummy,em(5,2),LDIM)
 }
}
else {
 do k = 1,LDIM {
  do l = 1,LDIM {
   em(k,l) = 0.
  }
 }
 if (i!=j) return
 call selfs(i,j,r2,iovpar,em,LDIM)
 if (lmax(i) <= 1) return
 call selfp(i,j,r2,iovpar,em(2,2),LDIM)
 if (lmax(i) <= 2) return
 call selfd(i,j,r2,iovpar,em(5,5),LDIM)
}
end
