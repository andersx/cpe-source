! ************************************************************************
! COMPUTE CONTRIBUTION OF CPE RESPONSE TO THE DFTB HAMILTONIAN
! implemented by S. Kaminski JPC A. 2012 (116), 9131-9141.
!
! There are two types of Coulomb type integrals to be evaluated: 
! 'M' is over 1 CPE basis function and 1 DFTB3 basis function
! 'N' is over 2 CPE basis functions
! by DFTB3 basis function I mean the Slater function for the charge
! density fluctiations \delta\rho; CPE basis functions are Gaussians, 
! which are used to represent the CPE response density (eq.12, eq. 17).
! For each atom, there are 3 CPE basis functions, which are similar to 
! p-orbitals (one can think of them as polarization functions). 
! ***********************************************************************
! K.W. 2014-05-01 

subroutine cpe_scf(nn,qmat,qzero,izp,x,xm,E_cpe,shiftCPE,C,dipol,cpedipol, & 
                   cpedipabs,uhder,uhubb,CPEPar)
  implicit none
  include 'maxima.h'

  !-> INPUT 
  integer, intent(in) :: nn
  integer, dimension(NNDIM), intent(in) :: izp
  double precision, dimension(NNDIM), intent(in) :: qmat
  double precision, dimension(3), intent(in) :: dipol
  double precision, dimension(3,NNDIM), intent(in) :: x
  double precision, dimension(MAXTYP), intent(in) ::  xm
  double precision, dimension(MAXTYP,3), intent(in) :: uhder, uhubb
  double precision, dimension(MAXTYP,4), intent(in) :: qzero

  ! OUTPUT -> 
  double precision, intent(out) :: E_cpe,cpedipabs      ! response energy contribution
  double precision, dimension(3*nn), intent(out) :: C   ! response dipole coefficients (zxy)
  double precision, dimension(nn), intent(out)   :: shiftCPE
  double precision, dimension(3), intent(out)    :: cpedipol
  double precision, dimension(NNDIM,4),intent(out) :: CPEPar

  ! LOCAL
  integer :: i,j,k,l,i_k,i_l,i1,i2,j1,j2
  double precision, dimension(3) :: Rab,Ga,dGa,Gb,dGb,Gc,dGc,Gd,dGd,Y
  double precision, dimension(3,3) :: P, Pa, Pb, dYda
  double precision, dimension(3,1) :: Pc, Pd, P2
  double precision, dimension(nn) :: tau, Q_qm
  double precision, dimension(3*nn,3*nn) :: A, dAdQa, dAdQb, AI
  double precision, dimension(3*nn,nn) :: Mmat, dMdQb
  double precision, dimension(3*nn) :: dMdQa, M, zw
  double precision :: s,f,Za,Zb,Zc,Zd,r2,r,b1,b2,br1,br2,c1,c2,j00,d00dr,d00dr2,j10
  double precision :: dedr,rlo,rhi,j11,e1,e2,dedb,dbdQa,dbdQb,d00db,dedQa,dedQb,de2db2
  double precision :: db2dQd,db2dQc,d10dQd,d10dQc,dedrdb,d00dr2db,d00drdb,d10db,d10dQa
  double precision :: d10dQb,d11dQa,d11dQb,dZadQa,dZbdQb,dZcdQc,dZddQd
  double precision, dimension(3), parameter :: XSZ = (/ 0.4068833884920483, &
                                                        0.09179674341953627, &
                                                        3.515225231758639 /)
  double precision, dimension(3), parameter :: XSC = (/ 0.3620527755096057, &
                                                        0.6262050528632612, &
                                                        0.01174217162750757 /)
  double precision, parameter :: SQRT_PI = 1.77245385
  double precision, parameter :: TO_DEBYE = 2.541765d0

  ! ********** assign CPE fit parameters (from Steve's paper) ****** 
  ! I think he used DFTB3, but with 'mio' parameters. So, I am not 
  ! sure, if one should re-fit these parameters !?

  do i = 1,nn
    Q_qm(i) = qzero(izp(i),4) - qmat(i)
    ! exponent of the DFTB Slater functions 
    tau(i) = (16.0d0/5.0d0) * (uhubb(izp(i),1) + uhder(izp(i),1) * Q_qm(i))

    select case ( nint(xm(izp(i))) )
      case (1) ! hydrogen
        CPEPar(i,1) = 4.10743330     ! Z
        CPEPar(i,2) = 4.37992780     ! B
        CPEPar(i,3) = 0.65017096     ! r_l
        CPEPar(i,4) = 0.89976532     ! r_u
        
      case (12) ! carbon
        CPEPar(i,1) = 2.12941380
        CPEPar(i,2) = 0.46271552
        CPEPar(i,3) = 1.58230960
        CPEPar(i,4) = 2.60816630

      case (16) ! oxygen
        CPEPar(i,1) = 4.59123620
        CPEPar(i,2) = 1.05271210
        CPEPar(i,3) = 4.14445140
        CPEPar(i,4) = 4.50938330

      case (14) ! nitrogen
        CPEPar(i,1) = 2.58954660
        CPEPar(i,2) = 0.50147210
        CPEPar(i,3) = 2.99983790
        CPEPar(i,4) = 3.12407570

      case (32) ! sulfur
        CPEPar(i,1) = 2.26117890
        CPEPar(i,2) = 0.74985889
        CPEPar(i,3) = 3.40497140
        CPEPar(i,4) = 624.334400

      case (31) ! phosphorous
        CPEPar(i,1) = 36.9483900
        CPEPar(i,2) = 105.240460
        CPEPar(i,3) = 74.7108830
        CPEPar(i,4) = 1998.27850
        
      case default
        print *, 'element ', i, ' is not recognized!'
        stop
    end select
  end do

  ! initialize some matrices involved in the integral evaluation
  A     = 0.d0
  AI    = 0.d0
  Mmat  = 0.d0
  dMdQa = 0.d0
  dMdQb = 0.d0
  dAdQa = 0.d0
  dAdQb = 0.d0

  ! evaluate the second order matrix (see paper for reference)
  do i = 1,nn
    i1 = (i-1)*3 + 1
    i2 = i1 + 2
    Za = CPEPar(i,1) * exp( CPEPar(i,2) * Q_qm(i) )
    dZadQa = CPEPar(i,2) * Za

    do j = 1,nn
      j1 = (j-1)*3 + 1
      j2 = j1 + 2

      Zb = CPEPar(j,1) * exp( CPEPar(j,2) * Q_qm(j) )
      dZbdQb = CPEPar(j,2) * Zb

      Rab = x(:,i) - x(:,j)

      ! SSJPP_IntQ (I have no idea, what that means, I just copied it)
      Ga = XSZ * Za**2
      Gb = XSZ * Zb**2
      dGa = XSZ * 2*Za * dZadQa
      dGb = XSZ * 2*Zb * dZbdQb

      ! GGJPP_IntQ
      P = 0.d0
      Pa = 0.d0
      Pb = 0.d0

      ! r2 = Rab(1)*Rab(1) + Rab(2)*Rab(2) + Rab(3)*Rab(3)
      r2 = dot_product(Rab,Rab)

      if (r2 < 1.e-25) then
        e1 = 0.d0
        dedQa = 0.d0
        dedQb = 0.d0

        do k = 1,3
          do l = 1,3
            b1 = sqrt( Ga(k) * Gb(l) / ( Ga(k) + Gb(l) ) )
            dbdQa = dGa(k) * (Gb(l) / (Ga(k) + Gb(l)))**2 / (2.0d0*b1)
            dbdQb = dGb(l) * (Ga(k) / (Ga(k) + Gb(l)))**2 / (2.0d0*b1)
            c1 = XSC(k) * XSC(l)
            e1 = e1 + c1 * 4.0d0 * b1**3 / ( 3.0d0 * SQRT_PI )
            dedb = c1 * 12.0d0 * b1**2 / ( 3.0d0 * SQRT_PI )
            dedQa = dedQa + dedb * dbdQa
            dedQb = dedQb + dedb * dbdQb
          end do
        end do

        do k = 1,3
          P(k,k) = e1
          Pa(k,k) = dedQa
          Pb(k,k) = dedQb
        end do
      else

        r = sqrt(r2)
        Y = Rab/r

        do k = 1,3
          do l = 1,3
           dYda(k,l) = -Y(k)*Y(l)/r
          end do
          dYda(k,k) = dYda(k,k) + 1.d0/r
        end do

        j10 = 0.d0
        j11 = 0.d0
        d10dQa = 0.d0
        d10dQb = 0.d0
        d11dQa = 0.d0
        d11dQb = 0.d0

        do k = 1,3
          do l = 1,3
            b1      = sqrt( Ga(k) * Gb(l) / (Ga(k) + Gb(l)) )
            dbdQa   = dGa(k) * (Gb(l) / (Ga(k) + Gb(l)))**2 / (2.0d0*b1)
            dbdQb   = dGb(l) * (Ga(k) / (Ga(k) + Gb(l)))**2 / (2.0d0*b1)
            br1     = b1 * r
            c1      = XSC(k) * XSC(l)
            e1      = (2.d0/SQRT_PI) * b1 * exp(-br1*br1)
            dedr    = -2*b1*br1 * e1
            dedb    = (2.0d0/SQRT_PI) * (1.0d0-2.0d0*br1*br1) * exp(-br1*br1)
            dedrdb  = -4*br1*e1 - 2*b1*br1*dedb
            j00     = erf(br1)/r
            d00dr   = e1/r - j00/r
            d00dr2  = dedr/r - e1/r2 - d00dr/r + j00/r2
            d00db   = (2.0d0/SQRT_PI) * exp(-br1*br1)
            d00drdb = dedb/r - d00db/r
            d00dr2db= dedrdb/r - dedb/r2 - d00drdb/r + d00db/r2
            j10     = j10 + c1 * d00dr
            j11     = j11 + c1 * d00dr2

            d10dQa  = d10dQa + c1 *  d00drdb * dbdQa
            d10dQb  = d10dQb + c1 *  d00drdb * dbdQb
            d11dQa  = d11dQa + c1 * d00dr2db * dbdQa
            d11dQb  = d11dQb + c1 * d00dr2db * dbdQb
          end do
        end do

        do k = 1,3 
          i_k = mod(k+1,3)+1
          do l = 1,3
            i_l = mod(l+1,3)+1

            P(k,l)  =  j11 * Y(i_k) * (-Y(i_l)) - j10 * dYda(i_k,i_l)
            Pa(k,l) =  d11dQa * Y(i_k) * (-Y(i_l)) - d10dQa * dYda(i_k,i_l)
            Pb(k,l) =  d11dQb * Y(i_k) * (-Y(i_l)) - d10dQb * dYda(i_k,i_l)

          end do
        end do
      end if

      A(i1:i2, j1:j2)     = P
      dAdQa(i1:i2, j1:j2) = Pa
      dAdQb(i1:i2, j1:j2) = Pb

    end do
  end do

  ! first order matrix
  do i = 1,nn
    i1 = (i-1)*3 + 1
    i2 = i1 + 2

    Zc = CPEPar(i,1) * exp(CPEPar(i,2) * (Q_qm(i)))
    dZcdQc = CPEPar(i,2) * Zc

    do j = 1,nn
      Rab = x(:,i) - x(:,j)
      Zd  = tau(j)
      dZddQd = (16.0d0/5.0d0)*(uhder(izp(j),1))

      ! SSJPS_IntQ
      Gc  = XSZ * Zc**2
      Gd  = XSZ * Zd**2
      dGc = XSZ * 2*Zc * dZcdQc
      dGd = XSZ * 2*Zd * dZddQd

      ! GGJPS_IntQ
      P2 = 0.d0
      Pc = 0.d0
      Pd = 0.d0

      r2  = dot_product(Rab,Rab)

      if (r2 < 1.e-25) then
        P2 = 0.d0
        Pc = 0.d0
        Pd = 0.d0
      else
        r = sqrt(r2)
        Y = Rab / r
        j10 = 0.d0
        d10dQc = 0.0d0
        d10dQd = 0.0d0

        do k = 1,3
          do l = 1,3
           b2      = sqrt( Gc(k) * Gd(l) / (Gc(k) + Gd(l)))
           db2dQc  = dGc(k) * (Gd(l) / (Gc(k) + Gd(l)))**2 / (2.0d0*b2)
           db2dQd  = dGd(l) * (Gc(k) / (Gc(k) + Gd(l)))**2 / (2.0d0*b2)
           br2     = b2 * r
           c2      = XSC(k) * XSC(l)
           e2      = (2.d0/SQRT_PI) * b2 * exp(-br2*br2)
           j00     = erf(br2) / r
           d00dr   = e2/r - j00/r
           d00db   = (2.0d0/SQRT_PI) * exp(-br2*br2)
           de2db2  = (2.0d0/SQRT_PI) * (1.0d0 - 2.0d0*br2*br2) * exp(-br2*br2)
           d00drdb = de2db2/r - d00db/r
           j10     = j10 + c2 * d00dr 
           d10db   = c2 * d00drdb 
           d10dQc  = d10dQc + d10db * db2dQc 
           d10dQd  = d10dQd + d10db * db2dQd 
          end do
        end do

        P2(1,1) = P2(1,1) + j10 * Y(3)
        P2(2,1) = P2(2,1) + j10 * Y(1)
        P2(3,1) = P2(3,1) + j10 * Y(2)
        Pc(1,1) = Pc(1,1) + d10dQc * Y(3)
        Pc(2,1) = Pc(2,1) + d10dQc * Y(1)
        Pc(3,1) = Pc(3,1) + d10dQc * Y(2)
        Pd(1,1) = Pd(1,1) + d10dQd * Y(3)
        Pd(2,1) = Pd(2,1) + d10dQd * Y(1)
        Pd(3,1) = Pd(3,1) + d10dQd * Y(2)
      endif

      rlo = CPEPar(i,3) + CPEPar(j,3)
      rhi = CPEPar(i,4) + CPEPar(j,4)

      ! Switching function
      if (r .ge. rhi) then
        s = 0.d0
      else if (r .le. rlo) then
        s = 1.d0
      else
        f = (rhi-r) / (rhi-rlo)
        s = 10.d0*f**3 - 15.d0*f**4 + 6.d0*f**5   
      end if
      s = 1.d0 - s
            
      P2 = s * P2
      Pc = s * Pc
      Pd = s * Pd
        
      Mmat(i1:i2, j)  = Mmat(i1:i2, j)  + P2(1:3,1)
      dMdQa(i1:i2)    = dMdQa(i1:i2)    + Pc(1:3,1) * Q_qm(j)
      dMdQb(i1:i2, j) = dMdQb(i1:i2, j) + Pd(1:3,1) * Q_qm(j)

    end do
  end do

  call invert_eta_cpe(3*nn, A, AI)

  E_cpe = 0.d0
  shiftCPE = 0.d0
  C     = 0.d0

  M = matmul(Mmat,Q_qm)
  C = -matmul(AI,M) 
  E_cpe = dot_product(C,M) + 0.5d0 * dot_product(C,matmul(A,C)) 

  do i = 1, nn
    i1 = (i-1)*3 + 1
    i2 = i1 + 2

    shiftCPE(i) = shiftCPE(i) + dot_product(C,Mmat(:,i)) + dot_product(C(i1:i2),dMdQa(i1:i2))

    do j = 1,nn
      shiftCPE(j) = shiftCPE(j) + dot_product(C(i1:i2),dMdQb(i1:i2,j))
    end do
  end do

  do i = 1,nn
    i1 = (i-1)*3 + 1
    i2 = i1 + 2

    do j = 1,nn
      j1 = (j-1)*3 + 1
      j2 = j1 + 2

      shiftCPE(i) = shiftCPE(i) + 0.50d0 * dot_product( C(i1:i2), matmul(dAdQa(i1:i2,j1:j2), C(j1:j2)) )
      shiftCPE(j) = shiftCPE(j) + 0.50d0 * dot_product( C(i1:i2), matmul(dAdQb(i1:i2,j1:j2), C(j1:j2)) )
    end do
  end do

  ! corrected dipole moment
  zw = 0.0d0

  do i = 1,nn
    j = (i-1)*3
    zw(1) = zw(1) + C(j+1)
    zw(2) = zw(2) + C(j+2)
    zw(3) = zw(3) + C(j+3)
  end do

  cpedipol = 0.0d0
  cpedipol(1) = dipol(1) + (zw(2) * TO_DEBYE)
  cpedipol(2) = dipol(2) + (zw(3) * TO_DEBYE)
  cpedipol(3) = dipol(3) + (zw(1) * TO_DEBYE)

  ! norm of the dipole moment
  cpedipabs = norm2(cpedipol)
end subroutine cpe_scf
