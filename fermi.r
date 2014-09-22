 subroutine FERMI(nel,telec,ndim,ev,dacc,occ,efermi,occperorb)
# MG: occperorb: maximal occupation number per orbital =1 for spinpolarized/ =2 for spinunpolarized calculations
#
# maxima definitions
#
   implicit REAL*8 (A-H,O-Z)
   include 'maxima.h'
#
   integer ndim
   real*8 nel,telec,efermi,dacc,ev(MDIM),occ(MDIM),occperorb
#
# Boltzmann constant in H/K
#
parameter(ckbol = 3.16679e-6)
#
# degeneracy tolerance for T = 0 
#
parameter(degtol = 1.0e-4) 


   real*8 beta,etol,occdg,chleft,ef1,ef2,ceps,eeps,charge,fac
   integer i,nef1,nef2,nup,ndown,nocc2,ndeg,istart,iend,hlp
   logical tzero
common /machine/ racc
#              machine accuracy

# start defining occupation numbers and their derivatives 
#
  do i = 1,ndim {
    occ(i) = 0.0d0
  }
  if (nel < 1.0e-5) goto 10
  if (nel > int(occperorb)*ndim) {
    print *,'too many electrons', nel,int(occperorb)*ndim
    stop
  }

# find energy interval for Fermi energy
#
  if (telec > 5.0d0) {
    beta = 1.0d0/(ckbol*telec) 
    etol = ckbol*telec*(log(beta)-log(racc)); tzero = .false.
  }
  else {
    etol = degtol; tzero = .true.
  }

# find energy range for Fermi energy
#
  hlp = int(occperorb)
  if(nel>int(nel)) {
   nef1 = (nel+hlp)/hlp; nef2 = (nel+hlp)/hlp
  }
  else {
   nef1 = (nel+hlp-1)/hlp; nef2 = (nel+hlp)/hlp
   # nef1 = (nel+1)/2; nef2 = (nel+2)/2
  }
  efermi = 0.5d0*(ev(nef1)+ev(nef2))
  nup = nef1; ndown = nef1;
  while (nup < ndim) {
    if (abs(ev(nup+1)-efermi) > etol) break;
    nup = nup+1
  }
  while (ndown > 0) {
    if (abs(ev(ndown)-efermi) > etol) break;
    ndown = ndown-1
  }
  ndeg = nup-ndown; nocc2 = ndown;
  do i = 1,nocc2 {
    occ(i) = occperorb
  }
  if (ndeg == 0) goto 10

# for T == 0 occupy orbitals as usual
#
  if (tzero) {
    occdg = ndeg; occdg = (nel-occperorb*nocc2)/occdg
    do i = nocc2+1,nocc2+ndeg {
      occ(i) = occdg
    }
  }

# for finite T, use Fermi distribution
#
  else {
    chleft = nel-occperorb*nocc2
    istart = nocc2+1; iend = istart+ndeg-1
    if (ndeg == 1) {
      occ(istart) = chleft
      goto 10
    }

# bracket and iterate Fermi energy by bisection
#
    ef1 = efermi-etol-degtol; ef2 = efermi+etol+degtol; 
    ceps = dacc*chleft; eeps = dacc*max(abs(ef1),abs(ef2))
    repeat {
      efermi = 0.5d0*(ef1+ef2); charge = 0.0d0;
      do i = istart,iend {
        occ(i) = occperorb/(1.0d0+dexp(beta*(ev(i)-efermi)))
        charge = charge+occ(i)
      }
      if (charge > chleft) ef2 = efermi; else ef1 = efermi
    }
    until ((abs(charge-chleft) < ceps) || (abs(ef1-ef2) < eeps))

# rescale occ(i) if the accuracy was not good enough 
#
    if (abs(charge-chleft) > ceps) {
      fac = chleft/charge;
      do i = istart,iend {
        occ(i) = occ(i)*fac
      } 
    }

  }
10 continue
end
