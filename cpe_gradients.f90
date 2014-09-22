! compute the contribution of CPE to the gradient
! re-implementation of Steve Kaminski's code
! K.W. 2014-05-02

subroutine cpe_grad(nn,izp,x,qmat,qzero,CPEPar,C,grad_cpe, uhubb, uhder)
  implicit none
  include 'maxima.h'

  ! -> INPUT
  integer, dimension(NNDIM), intent(in) :: izp
  integer, intent(in) :: nn
  double precision, dimension(NNDIM), intent(in) :: qmat
  double precision, dimension(MAXTYP,4), intent(in) :: qzero
  double precision, dimension(MAXTYP,3), intent(in) :: uhubb, uhder
  double precision, dimension(3*nn), intent(in) :: C 
  double precision, dimension(3,NNDIM), intent(in) :: x
  double precision, dimension(NNDIM,4), intent(in) :: CPEPar

  ! OUTPUT -> 
  double precision, dimension(3,NNDIM), intent(out) :: grad_cpe

  ! LOCAL
  integer :: i,j,k,l,m,i_k,i_l,i1,i2,j1,j2

  double precision, parameter :: SQRT_PI = 1.7724538
  double precision, dimension(3), parameter :: XSZ = (/ 0.4068833884920483, &
                                                        0.09179674341953627, & 
                                                        3.515225231758639 /)
  double precision, dimension(3), parameter :: XSC = (/ 0.3620527755096057, &
                                                        0.6262050528632612, &
                                                        0.01174217162750757 /)
  double precision, dimension(3) ::  Y, dEdRa, Rab, Ga, Gb
  double precision, dimension(3,1) :: P2
  double precision, dimension(3,3) :: dYda, P
  double precision, dimension(3,3,1) :: dP2
  double precision, dimension(3,3,3) :: dYda2(3,3,3), dP
  double precision, dimension(nn) :: Q_qm, Zeta_qm
  double precision :: j00,d00dr,d00dr2,d00dr3,j10,j11,j21,r2,r3,s,f,sf,ds,dsf
  double precision :: dedr,dedr2,dedr4,e1,e2,c1,c2,b1,b2,br1,br2,r,rlo,rhi,Za,Zb

  grad_cpe = 0.0d0

  do i = 1,nn
    Q_qm(i) = qzero(izp(i),4) - qmat(i)
  end do

  do i = 1,nn
    zeta_qm(i)  = (16.0d0/5.0d0) *(uhubb(izp(i),1) + uhder(izp(i),1) * Q_qm(i))
  end do

  ! 2nd order gradient
  do i = 1,nn
    i1 = (i-1)*3 + 1
    i2 = i1 + 2
    Za = CPEPar(i,1) * exp(CPEPar(i,2) * Q_qm(i))

    do j = 1,nn
      j1 = (j-1)*3 + 1
      j2 = j1 + 2
      Zb = CPEPar(j,1) * exp(CPEPar(j,2) * Q_qm(j))

      Rab = x(:,i) -x(:,j)

      Ga  = XSZ * Za**2
      Gb  = XSZ * Zb**2
      P   = 0.0d0
      dP  = 0.0d0

      r2  = dot_product(Rab,Rab)

      e1 = 0.d0
      if (r2 < 1.E-25) then

        do k = 1,3
          do l =1,3
             b1 = sqrt( Ga(k) * Gb(l) / (Ga(k) + Gb(l)) )
             c1 = XSC(k) * XSC(l)
             e1 = e1 + c1 * 4.d0 * b1**3 / (3.d0 * SQRT_PI)
          end do
        end do

        do k = 1,3
          P(k,k) = e1
        end do
      else
        r  = sqrt(r2)
        r3 = r * r2
        Y  = Rab / r

        do k = 1,3
          do l = 1,3
            dYda(k,l) = -Y(k) * Y(l) / r
          end do
          dYda(k,k) = dYda(k,k) + 1.d0 / r
        end do

        dYda2 = 0.0d0

        do k = 1,3
          do l = 1,3
            do m = 1,3
              dYda2(k,l,m) = -dYda(k,m) * Y(l)/r  &
                             -Y(k) * dYda(l,m)/r + (Y(k) * Y(l)/r2) * Y(m)
            end do
          end do
          do l = 1,3
            dYda2(k,k,l) = dYda2(k,k,l) - Y(l)/r2
          end do
        end do

        j10    = 0.0d0
        j11    = 0.0d0
        j21    = 0.0d0
        d00dr  = 0.0d0
        d00dr2 = 0.0d0
        d00dr3 = 0.0d0

        do k = 1,3
          do l = 1,3
            b1  = sqrt( Ga(k) * Gb(l) / (Ga(k) + Gb(l)) )
            br1 = b1 * r
            c1  = XSC(k) * XSC(l)
            j00 = erf(br1)/r
            e1  = (2.0d0/SQRT_PI) * b1  * exp(-br1*br1)
            dedr = -2*b1*br1 * e1
            dedr2 = -2*b1*b1 * e1 - 2.0d0*b1*br1 * dedr
            d00dr = e1/r - j00/r
            d00dr2 = dedr/r - e1/r2 - d00dr/r + j00/r2
            d00dr3 = dedr2/r - dedr/r2 - dedr/r2 + 2*e1/r3 & 
                   - d00dr2/r + d00dr/r2 + d00dr/r2 - 2*j00/r3
            j10    = j10 + c1 * d00dr 
            j11    = j11 + c1 * d00dr2
            j21    = j21 + c1 * d00dr3
          end do
        end do

        do k = 1,3
          i_k = mod(k+1,3) + 1
          do l = 1,3
            i_l = mod(l+1,3) + 1
            P(l,k) = j11 * Y(i_l) * (-Y(i_k)) - j10 * dYda(i_l,i_k)
          end do
        end do

        do k = 1,3
          i_k = mod(k+1,3) + 1
          do l = 1,3
            i_l = mod(l+1,3) + 1
            do m = 1,3
              dP(m,l,k) = j21 * Y(m) * Y(i_l) * (-Y(i_k)) + j11 & 
                        * dYda(i_l,m) * (-Y(i_k)) + j11 * Y(i_l) & 
                        * (-dYda(i_k,m)) - j11 * Y(m) * dYda(i_l,i_k) & 
                        - j10 * dYda2(i_l,i_k,m)
            end do
          end do
        end do
      end if

      do k = 1,3
        dEdRa(k) = 0.5d0 * dot_product( C(i1:i2), matmul(dP(k,:,:),C(j1:j2)) )
      end do

      grad_cpe(:,i) = grad_cpe(:,i) + dEdRa
      grad_cpe(:,j) = grad_cpe(:,j) - dEdRa 
    end do
  end do

  ! 1st order gradient
  do i = 1,nn
    i1 = (i-1)*3 + 1
    i2 = i1 + 2
    Za = CPEPar(i,1) * exp(CPEPar(i,2) * Q_qm(i))

    do j = 1,nn
      Zb  = Zeta_qm(j)
      Rab = x(:,i) - x(:,j)
      Ga  = XSZ * Za**2
      Gb  = XSZ * Zb**2
      P2  = 0.0d0
      dP2 = 0.0d0

      r2  = dot_product(Rab,Rab)

      if (r2 < 1.e-25) then
        P2 = 0.0d0
        dP2 = 0.0d0
      else
        r = sqrt(r2)
        Y = Rab/r

        do k = 1,3
          do l = 1,3
            dYda(k,l) = -Y(k)*Y(l)/r
          end do
            dYda(k,k) = dYda(k,k) + 1.d0/r
        end do

        j10 = 0.0d0
        j11 = 0.0d0
        j00 = 0.0d0
        d00dr = 0.0d0
        d00dr2 = 0.0d0
  
        do k = 1,3
          do l = 1,3
             b2 = sqrt(Ga(k) * Gb(l) / (Ga(k) + Gb(l)))
             br2 = b2 * r
             c2 = XSC(k) * XSC(l)
             e2 = (2.0d0/SQRT_PI) * b2 * exp(-br2*br2)
             dedr4  = -2*b2*br2 * e2
             j00    = erf(br2)/r
             d00dr  = e2/r - j00/r
             d00dr2 = dedr4/r - e2/r2 - d00dr/r + j00/r2
             j10    = j10 + c2 * d00dr 
             j11    = j11 + c2 * d00dr2
          end do
        end do

        P2(1,1) = j10 * Y(3)
        P2(2,1) = j10 * Y(1)
        P2(3,1) = j10 * Y(2)   
       
        do k = 1,3
           dP2(k,1,1) = j11 * Y(k) * Y(3) + j10 * dYda(3,k)
           dP2(k,2,1) = j11 * Y(k) * Y(1) + j10 * dYda(1,k)
           dP2(k,3,1) = j11 * Y(k) * Y(2) + j10 * dYda(2,k)
        end do
      end if

      ! switching function
      rlo = CPEPar(i,3) + CPEPar(j,3)
      rhi = CPEPar(i,4) + CPEPar(j,4)

      if (r .ge. rhi) then
        s  = 0.0d0
        ds = 0.0d0
      else if (r .le. rlo) then
        s  = 1.0d0
        ds = 0.0d0
      else
        f  = (rhi-r)/(rhi-rlo)
        ds = -1.0d0 / (rhi-rlo)
        sf = 10.0d0 * f**3 - 15.0d0 * f**4 + 6.0d0 * f**5
        dsf = 30.0d0 * f**2 - 60.0d0 * f**3 + 30.0d0 * f**4
        s = sf
        ds = dsf * ds
      end if

      s  = 1.0d0 - s
      ds = - ds
 
      do k = 1,3
        dEdRa(k) = s * dot_product( C(i1:i2),dP2(k,:,1)*Q_qm(j) )
      end do

      if (r > 0.0d0) then
        e2 = dot_product( C(i1:i2), P2(:,1)*Q_qm(j) )
        dEdRa = dEdRa + ( ds * Rab/r ) * e2
      end if

      grad_cpe(:,i) = grad_cpe(:,i) + dEdRa
      grad_cpe(:,j) = grad_cpe(:,j) - dEdRa

    end do
  end do
end subroutine cpe_grad
