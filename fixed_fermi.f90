! This subroutine gets called, when we want to evaluate the energy of a system
! with a fixed occupation number in the HOMO and LUMO, e.g. to calculate
! the IP of an excited state. Requires the presence of a small file 'occ.dat' 
! with 2 lines. E.g.
! HOMO: 2.0000
! LUMO: 0.0000
! -- K.W. 2014-04-26

subroutine fixed_fermi(nel,ndim,occ,efermi,occperorb)
    implicit none
    include 'maxima.h'
    !<- IN
    integer,intent(in) :: ndim
    double precision, intent(in)  :: nel,occperorb
    ! OUT -> 
    double precision, dimension(MDIM), intent(out) :: occ
    double precision, intent(out) :: efermi

    ! LOCAL
    double precision  :: homo, lumo
    character(len=10) :: orbital
    integer :: i_homo,i
    
    ! Boltzmann in Hartree/Kelvin
    double precision, parameter :: ckbol = 3.16679e-6
    ! degeneracy tolerance for T = 0
    double precision, parameter :: degtol = 1.0e-04

    ! read in desired configuration for HOMO/LUMO
    open (unit=77, file='occ.dat', status='OLD')
      read (77,*) orbital, homo
      read (77,*) orbital, lumo
    close (77) 

    ! do not calculate efermi here, does it matter somewhere?
    efermi = 0.0d0

    ! some system checks
    if (nel < 1.0e-5) then
        stop 'Too few electrons!'
    endif

    if (nel > int(occperorb) * ndim) then
        stop 'too many electrons!'
    endif

    ! number of the homo
    i_homo = int(nel) / 2

    do i = 1,ndim 
      if (i < i_homo) then
          occ(i) = 2.0d0
      else if (i == i_homo) then
          occ(i) = homo
      else if (i == i_homo+1) then
          occ(i) = lumo
      else 
          occ(i) = 0.0d0
      end if
    end do
end subroutine fixed_fermi
