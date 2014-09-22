   subroutine MULLIKEN(nn,x,izp,qmat,qzero,qmulli,ql,qtot,ndim,dacc,occ,_
                       a,overl,lmax,ind,dipol,dipabs)
#
# maxima definitions
#
   implicit REAL*8 (A-H,O-Z)
   include 'maxima.h'
#
      integer nn,ndim,lmax(MAXTYP),izp(NNDIM),ind(NNDIM+1)
      real*8 dacc,occ(MDIM),a(MDIM,MDIM),overl(MDIM,MDIM)
      real*8 x(3,*),qzero(MAXTYP,4),qtot,qmulli(MDIM)
      real*8 qmat(NNDIM),dipol(3),dipabs,afoo(MDIM,MDIM)
      real*8 ql(3*NNDIM)

      integer i,j,lj,m,n,izpj,adim,lcount
      real*8 sum,qhelp,conv

      do i = 1,ndim {
        if (occ(i) < dacc) break      
      }
      adim = i-1
      call dsymm('L','U',ndim,adim,1.0d0,overl,MDIM,a,MDIM,0.0d0,afoo,MDIM)
      do n = 1,ndim {
        qmulli(n) = 0.0d0
        do i = 1,adim {
          qmulli(n) = qmulli(n) + occ(i)*afoo(n,i)*a(n,i)
        }
      }      

      lcount=0    
      do j = 1,nn {
        qtot = 0.0d0; 
        do lj = 1,lmax(izp(j)) {
          jofn = ind(j)+(lj-1)**2; qhelp = 0.0d0
          do mj = 1,2*lj-1 {
           qhelp = qhelp + qmulli(jofn+mj)
          }
          lcount=lcount+1
          ql(lcount) = qhelp
          qtot = qtot + qhelp
        }
        qmat(j) = qtot
      }

# calculation of total charge 
#
    qtot = 0.0d0
    do j = 1,nn {
        qtot = qtot+qmat(j)
    }

# calculation of dipol moment
#
    do i = 1,3 {
      dipol(i) = 0.d0
      do j = 1,nn {
        izpj = izp(j)
        qhelp = qzero(izpj,4) - qmat(j)
	dipol(i) = dipol(i) + qhelp*x(i,j)
      }
# conversion to debye
      dipol(i) = dipol(i)*2.541765d0
    }

# norm of dipol moment
#
    dipabs = 0.d0
    do i = 1,3 {
      dipabs = dipabs + dipol(i)**2
    }
    dipabs = sqrt(dipabs)

end




