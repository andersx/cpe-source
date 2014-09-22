subroutine usualgrd(nn,nbeweg,ndim,izp,lmax,ind,period,doscf, _
                    x,ev,a,b,occ,shift,shift3,shift3A,shift3B,dacc, _
                    boxsiz,xinvbox,grad,ldep,spin,spinfac,spinshift)
#
#
#implicit REAL*8 (A-H,O-Z)
implicit none
include 'maxima.h'
#
real*8 deltax
parameter ( deltax = 0.01d0 )

#
# input parameter
#
logical period,doscf,spin,ldep
integer nn,nbeweg,ndim,izp(*),lmax(*),ind(*)
real*8 x(3,*),ev(*),a(MDIM,MDIM),b(MDIM,MDIM),occ(*),grad(3,*)
real*8 dacc,boxsiz(3,3),xinvbox(3,3),shift(NNDIM,3),shift3(NNDIM,3),shift3A(NNDIM,3),shift3B(NNDIM)
real*8 spinfac,spinshift(NNDIM,3)

#
# Parameter used local
#
integer m,n,i,j,k,lj,lk,mj,mk
integer mu,nu,izpj,izpk,indj,indk
real*8 au(LDIM,LDIM),bu(LDIM,LDIM),ocmcc,dgrh,rcdx,xhelp,dgrs,dgr
real*8 auh(LDIM,LDIM),buh(LDIM,LDIM),dgr3,dgrspin

character*1  sccmode
common /scchelp/ sccmode


rcdx=1.0d0/deltax

do i=1,nn {
  do k=1,3 {
    grad(k,i) = 0.0d0
  }
}
dgr=0.0d0
dgr3=0.0d0

# setup of density and energy weighted matrix 
#
do n = 1,ndim {
  do m = 1,ndim {
    b(m,n) = 0.0d0
  }
}

do i = 1,ndim {
  if (occ(i) < dacc) break
  do m = 1,ndim {
    do n = 1,m-1 {
      ocmcc = occ(i)*a(m,i)*a(n,i)
      b(n,m) = b(n,m) + ocmcc*ev(i)   ### upper triangle
      b(m,n) = b(m,n) + ocmcc         ### lower triangle
    }
  }
}
do n = 1,ndim {
  do m = 1,ndim {
    if (abs(b(m,n)) < dacc) b(m,n) = 0.0d0
  }
}

#
# gradient: normal contribution
#
# for the atom on position l" the force contribution of the hamilton matrix 
# can be written as (The sum over the i for the different c_i can be 
# evaluated earlier and is included in the coefficients C):
#
# f =  sum_{l,l'} sum_{mu,nu} C_{l,mu}  C_{l,nu} (d H_{l,l'}^{mu,nu}/ d R_{l"})
#
#  l,l' are the atom positions mu,nu the orbitals
#  H_{l,l'}^{mu nu} = < phi_mu(r-R_{l}) | H | phi_nu(r-R_{l'}) >
#
# the term (d H_{l,l'}^{mu,nu}/ d R_{l"}) is not zero only if (l or l' are
# equal l") and (l not equal l')
#
# we have also (d H_{l,l"}^{mu,nu}/ d R_{l"}) = (d H_{l",l}^{nu,mu}/ d R_{l"})
#
# since H_{l,l"}^{mu,nu} = H_{l",l}^{nu,mu}
# 
# so that we can write
#
# f =sum_(l neq l") sum_{mu,nu} C_{l,mu}C_{l,nu} 2*(d H_{l,l"}^{mu nu}/d R_{l"})
#
# If you want to sum over either index l and l', you have once to replace
# dif by -dif in the call of slkode! 
#

do j = 1,nbeweg {       ### for every movable atom
  indj=ind(j)
  izpj=izp(j)
  do k = 1,nn {         ### for every atom that acts on the moveable atom
    if(k!=j) {
    indk=ind(k)
    izpk=izp(k)
    
      do i = 1,3 {          ### for every spatial component


# hamilton and overlap matrix contribution
#
        xhelp  = x(i,j)
        x(i,j) = xhelp + deltax       
        call slkmatrices(k,j,x,au,bu)
        x(i,j) = xhelp - deltax
        call slkmatrices(k,j,x,auh,buh)
        x(i,j) = xhelp
        
#
#   use sumation over angular momentum and magnetic quantum numbers
#   since shift is actually defined for the orbitals 
#
        do lj = 1,lmax(izpj) {
         do mj=1,2*lj-1 {
          n  = (lj-1)**2 + mj
          nu = n + indj
          do lk = 1,lmax(izpk) {
            do mk = 1,2*lk-1 {
              m  = (lk-1)**2 + mk
              mu =  m + indk

#
# dgrh = 2 * ( d H_{k,j}^{m,n}/ d R_{j} )
#
              dgrh = (au(m,n)-auh(m,n))*rcdx
#
# dgrs = - 2 * ( d S_{k,j}^{m,n}/ d R_{j} )
#
              dgrs = -(bu(m,n)-buh(m,n))*rcdx
#
# dgr  =  ( d S_{k,j}^{m,n}/ d R_{j} ) * (shift(k)+shift(j))
# dgr3 =  ( d S_{k,j)^{m,n}/ d R_{j} ) * (2*shift3(k)+shift3A(k))/3.0 (for ldep different) 
#
              if (ldep){            
                if (doscf)  dgr = -0.5d0*dgrs*(shift(k,lk)+shift(j,lj))
                if (sccmode=="3") dgr3=-0.5d0*dgrs/3.0d0* _
                    (shift3(k,lk)+shift3A(k,lk)+shift3(j,lj)+shift3A(j,lj)+shift3B(k)+shift3B(j))
              }
              else {
                if (doscf)  dgr = -0.5d0*dgrs*(shift(k,1)+shift(j,1))
                if (sccmode=="3") dgr3=-0.5d0*dgrs/3.0d0* _
                    (2.0d0*shift3(k,1)+shift3A(k,1)+2.0d0*shift3(j,1)+shift3A(j,1))
              }
              if (spin)   dgrspin = -0.5d0*dgrs*(spinshift(k,lk)+spinshift(j,lj))
#
# only lower triangle contains sum_i n(i) c_{mu,i} c_{nu,i}
# only upper triangle contains sum_i epsilon(i) n(i) c_{mu,i} c_{nu,i}
#
              if(mu>nu) {
               dgrh = dgrh * b(mu,nu)
               dgrs = dgrs * b(nu,mu)
               if(doscf)         dgr = dgr * b(mu,nu)
               if (sccmode=="3") dgr3= dgr3* b(mu,nu)
               if(spin)       dgrspin= dgrspin*b(mu,nu)
               }
              else {
               dgrh = dgrh * b(nu,mu)
               dgrs = dgrs * b(mu,nu)
               if(doscf) dgr = dgr * b(nu,mu)
               if (sccmode=="3") dgr3=dgr3*b(nu,mu)
               if(spin)       dgrspin= dgrspin*b(nu,mu)
              }
              grad(i,j) = grad(i,j) + dgrh + dgrs
              if(doscf)        grad(i,j) = grad(i,j) + dgr
              if(sccmode=="3") grad(i,j) = grad(i,j) + dgr3 
              if(spin)         grad(i,j) = grad(i,j) + spinfac*dgrspin
            }
          }
         }
        }
      
      }
    }
  }
}

end
