      integer nmax,nmin
      parameter (nmax = 10)
      parameter (nmin = 2)
      real*8 gmat,sumlat
      common /comgmat/gmat(4,2*nmax+1,2*nmax+1,2*nmax+1)
      common /comrmat/sumlat(3,2*nmax+1,2*nmax+1,2*nmax+1)
