subroutine invert_eta_cpe(N,A,AI)
  implicit none

  ! -> INPUT
  integer, intent(in) :: N
  double precision, dimension(N,N), intent(out) :: A

  ! OUTPUT -> 
  double precision, dimension(N,N), intent(out) :: AI

  ! LOCAL
  integer :: info,lwork,i,j
  double precision, parameter :: tol = 1.0e-10
  double precision, dimension(N,N) :: U, VT, WM
  double precision, dimension(N) :: W
  double precision, dimension(:), allocatable :: work
  double precision :: query 

  logical, parameter :: only_use_svd = .false.
  logical use_svd 

  AI = A
  info = 0
  use_svd = only_use_svd


  if (.not. use_svd) then
    ! cholesky factorization of matrix AI
    call dpotrf('U', N, AI(1,1), N, info)
    if (info < 0) then
      write (6,*) 'dpotrf: argument had an illegal value:', info
    else if (info > 0) then
      write (6,*) 'dpotrf: factorization not completed:', info
    else
      ! inversion of the cholesky factorized matrix, however
      ! why only the first element???
      call dpotri('U', N, AI(1,1), N, info)
      if (info < 0) then
        write (6,*) 'dpotri: argument had an illegal value:', info
      else if (info > 0) then
        write (6,*) 'dpotri: factorization not completed:', info
        use_svd = .true.
      else
        do j = 2,N
          do i = 1,j-1
            AI(j,i) = AI(i,j)
          end do
        end do
      end if
    end if
  end if

  if (use_svd) then
    lwork = -1
    call dgesvd("A","A",N,N,AI(1,1),N,W(1),U(1,1),N,VT(1,1),N,query,lwork,info)
    if (info /= 0) then
      write (6, *) 'dgesvd not completed' 
      stop
    end if

    lwork = nint(query)
    allocate( work(lwork) )
    W  = 0.d0
    U  = 0.d0
    VT = 0.d0

    call dgesvd("A","A",N,N, AI(1,1),N,W(1),U(1,1),N,VT(1,1), N, work(1), lwork, info)
    WM = 0.d0
    do i=1,N
      if ( abs(w(i)) > tol) then
        wm(i,i) = 1.0d0 / w(i)
      end if
    end do

    AI = matmul( transpose(VT), matmul(WM,transpose(U)) )

  end if
end subroutine invert_eta_cpe
