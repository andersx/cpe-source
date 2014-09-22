! *********************************************************
! calculation of CM3 charges according to: 
! Kalinowski et al. JPC A 108 (2004), 2545-2549.      
! correction parameters fitted for DFTB3 taken from:
! Kaminski et al. JPC A 116 (2012), 11927-11937.
! Steve Kaminski did the original implementation into DFTB
! 2013-06-19
! *********************************************************

subroutine cm3charges(nn,izp,qcm3,qmat,ndim,x,qzero, &
        cm3dipol,cm3dipabs,ind,occ,dacc,a,overl,cbnd, &
        number1,number2,order,xm,bord,qcm3tot,lmax)

        implicit none
        include 'maxima.h'
        
        ! -> INPUT
        integer,dimension(NNDIM),intent(in) :: izp
        integer,dimension(NNDIM+1),intent(in) :: ind
        integer,dimension(MAXTYP),intent(in) :: lmax
        integer,intent(in) :: nn,ndim
        double precision,dimension(NNDIM),intent(in) :: qmat
        double precision,dimension(MDIM),intent(in) :: occ
        double precision,dimension(MAXTYP),intent(in) :: xm
        double precision,dimension(MAXTYP,4),intent(in) :: qzero
        double precision,dimension(3,nn),intent(in) :: x
        double precision,dimension(MDIM,MDIM), intent(in) :: a,overl
        double precision,intent(in) :: dacc

        ! <- OUTPUT
        integer,dimension(NNDIM),intent(out) :: number1,number2
        integer,intent(out) :: cbnd
        double precision,dimension(NNDIM),intent(out) :: qcm3,order
        double precision,dimension(NNDIM,NNDIM),intent(out) :: bord
        double precision,dimension(3),intent(out) :: cm3dipol
        double precision,intent(out) :: qcm3tot,cm3dipabs 

        !  LOCAL
        integer :: i,j,k,m,n,test,ind_i,ind_j,l_i,l_j
        double precision,dimension(NNDIM,NNDIM) :: corrfact
        double precision,dimension(MDIM,MDIM)   :: dens,bondorder_per_orbital,ps
        double precision :: corrhelp

        ! correction parameter, depend linear on bond order 
        !("V" in Kaminski et al., Table 1)
        double precision, parameter :: D_NC = -0.066148350d0
        double precision, parameter :: D_NO =  0.000021086d0
        double precision, parameter :: D_NH = -0.116000680d0 
        double precision, parameter :: D_CO =  0.052336804d0
        double precision, parameter :: D_CH = -0.036270959d0
        double precision, parameter :: D_OH =  0.000006479d0 
        double precision, parameter :: D_SC =  0.000024998d0 
        double precision, parameter :: D_SH =  0.029239804d0 
        double precision, parameter :: D_SN =  0.142501660d0 
        double precision, parameter :: D_SP =  0.018375336d0 
        double precision, parameter :: D_SO =  0.080901696d0 
        double precision, parameter :: D_PC =  0.083956737d0 
        double precision, parameter :: D_PH = -0.032020411d0  
        double precision, parameter :: D_PN =  0.065625396d0  
        double precision, parameter :: D_PO =  0.101318670d0  
        ! correction parameter, depends quadratically on bond border
        !("W" in Kaminski et al., Table 1)
        double precision, parameter :: C_CO = -0.034477992d0
        double precision, parameter :: C_SP = -0.004322256d0
        double precision, parameter :: C_PO =  0.006754956d0
        double precision, parameter :: C_NO = -0.027890339d0

        ! Initializations
        ps    = 0.0d0
        bondorder_per_orbital = 0.0d0
        bord  = 0.0d0
        dens  = 0.0d0
        qcm3  = 0.0d0
        cbnd  = 0

        ! build density matrix
        do i = 1,ndim
            ! sum over occupied orbitals only
            if (occ(i) < dacc) then
                exit
            else
                do m = 1,ndim
                    do n = 1,ndim
                        dens(m,n) = dens(m,n) + occ(i) * a(m,i) * a(n,i)
                    end do
                end do
            end if
        end do

        ! Multiply density matrix with overlap matrix
        do i = 1,ndim
            do j = 1,ndim
                do k = 1,ndim
                    ps(i,j) = ps(i,j) + dens(i,k) * overl(k,j)
                end do
            end do
        end do

        ! building the bondorder matrix, which is PS * PS^T
        !bondorder_per_orbital = ps * transpose(ps)
        ! -> the explicit double loop is about 2x faster than
        ! the intrinsic statement above 
        do i = 1,ndim
            do j = 1,ndim
                bondorder_per_orbital(i,j) = ps(i,j) * ps(j,i)
            end do
        end do
       

        ! build the bond order matrix (per atoms) from the bondorder_per_orbital matrix (for orbitals)
        do i = 1,nn
            do j = 1,nn
                ! loop over orbitals for each atom
                do l_i = 1, lmax(izp(i))**2
                    ind_i = ind(i) + l_i
                    do l_j = 1, lmax(izp(j))**2
                        ind_j = ind(j) + l_j
                        bord(i,j) = bord(i,j) + bondorder_per_orbital(ind_i,ind_j)
                    end do
                end do

                if ( (abs(bord(i,j)) >= 0.75) .and. (i /= j) ) then
                    cbnd = cbnd + 1
                    order(cbnd) = bord(i,j) 
                    number1(cbnd) = i
                    number2(cbnd) = j

                    ! now assign the fit parameters based on the mass (is there
                    ! no better criteria in stand alone DFTB???
                    test = 1000 * nint(xm(izp(i))) + nint(xm(izp(j)))
                    select case(test)
                        ! C-H bond
                        case(12001)
                            corrfact(i,j) = -D_CH * bord(i,j)
                        ! H-C bond
                        case(1012) 
                            corrfact(i,j) =  D_CH * bord(i,j)
                        ! O-C bond
                        case(16012) 
                            corrfact(i,j) = -D_CO * bord(i,j) &
                                           - C_CO * bord(i,j)**2
                        ! C-O bond
                        case(12016) 
                            corrfact(i,j) =  D_CO * bord(i,j) &
                                           + C_CO * bord(i,j)**2
                        ! O-H bond
                        case(16001) 
                            corrfact(i,j) = -D_OH * bord(i,j)
                        ! H-O bond
                        case(1016)
                            corrfact(i,j) =  D_OH * bord(i,j)
                        ! N-H bond
                        case(14001)
                            corrfact(i,j) = -D_NH * bord(i,j)
                        ! H-N bond
                        case(1014)
                            corrfact(i,j) =  D_NH * bord(i,j)
                        ! N-C bond
                        case(14012)
                            corrfact(i,j) = -D_NC * bord(i,j)
                        ! C-N bond
                        case(12014)
                            corrfact(i,j) =  D_NC * bord(i,j)
                        ! N-O bond
                        case(14016)
                            corrfact(i,j) = -D_NO * bord(i,j) &  
                                            -C_NO * bord(i,j)**2
                        ! O-N bond
                        case(16014)
                            corrfact(i,j) =  D_NO * bord(i,j) &
                                           + C_NO * bord(i,j)**2 
                        ! S-O bond
                        case(32016)
                            corrfact(i,j) = -D_SO * bord(i,j) 
                        ! O-S bond
                        case(16032)
                            corrfact(i,j) =  D_SO * bord(i,j)
                        ! S-N bond
                        case(32014)
                            corrfact(i,j) = -D_SN * bord(i,j)
                        ! N-S bond
                        case(14032)
                            corrfact(i,j) =  D_SN * bord(i,j)
                        ! S-C bond
                        case(32012)
                            corrfact(i,j) = -D_SC * bord(i,j) 
                        ! C-S bond
                        case(12032)
                            corrfact(i,j) =  D_SC * bord(i,j)
                        ! S-H bond
                        case(32001)
                            corrfact(i,j) = -D_SH * bord(i,j)
                        ! H-S bond
                        case(1032)
                            corrfact(i,j) =  D_SH * bord(i,j)
                        ! P-O bond
                        case(31016)
                            corrfact(i,j) = -D_PO * bord(i,j) &
                                            -C_PO * bord(i,j)**2
                        ! O-P bond
                        case(16031)
                            corrfact(i,j) =  D_PO * bord(i,j) & 
                                           + C_PO * bord(i,j)**2
                        ! P-N bond
                        case(31014)
                            corrfact(i,j) = -D_PN * bord(i,j)
                        ! N-P bond
                        case(14031)
                            corrfact(i,j) =  D_PN * bord(i,j)
                        ! P-C bond
                        case(31012)
                            corrfact(i,j) = -D_PC * bord(i,j)
                        ! C-P bond
                        case(12031)
                            corrfact(i,j) =  D_PC * bord(i,j)
                        ! P-H bond
                        case(31001)
                            corrfact(i,j) = -D_PH * bord(i,j)
                        ! H-P bond
                        case(1031)
                            corrfact(i,j) =  D_PH * bord(i,j)
                        ! S-P bond (very common bonding motive ;)
                        case(32031)
                            corrfact(i,j) = -D_SP * bord(i,j) &
                                           - C_SP * bord(i,j)**2
                        ! P-S bond 
                        case(31032)
                            corrfact(i,j) =  D_SP * bord(i,j) &
                                           + C_SP * bord(i,j)**2
                    end select
                    corrhelp = corrhelp + corrfact(i,j)
                end if
          end do
          qcm3(i) = qmat(i) + corrhelp 
          corrhelp = 0.0
    end do

    ! remove diagonal elements of bond order matrix
    forall (i = 1:nn)
        bord(i,i) = 0.0d0
    end forall

    ! calculation of CM3 dipol moment
    cm3dipol = 0.0d0
    do j = 1,nn
       cm3dipol = cm3dipol + (qzero(izp(j),4)- qcm3(j)) * x(:,j)
    end do

    !atomic units to Debye
    cm3dipol = cm3dipol * 2.541765d0

    ! norm of dipol moment
    cm3dipabs = norm2(cm3dipol)
    !cm3dipabs = sqrt( cm3dipol(1)**2 + cm3dipol(2)**2 + cm3dipol(3)**2 )

    ! total charge
    qcm3tot = sum(qcm3)

end subroutine cm3charges
