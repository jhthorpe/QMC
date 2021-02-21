!-------------------------------------------------------
! linal.f90
!	- module containing linear algebra subroutines
!-------------------------------------------------------
module linal

contains

!-------------------------------------------------------
! ddet
!       - computes the determinant of a double
!         precision matrix via LU decomposition 
!
!       - note that matrix A is destroyed in the process
!-------------------------------------------------------
! M             : int, number of elements
! A             : 2D real*8, matrix to get det of
! IPIV		: 1D int, list of pivots

real(kind=8) function ddet(M,A,IPIV)
  implicit none
  integer, intent(in) :: M
  real(kind=8), intent(inout) :: A(M,M)
  integer, intent(inout) :: IPIV(M)

  real(kind=8) :: dtmp
  integer :: i,info

  if (M .eq. 1) then
    ddet = A(1,1)
    return
  end if

  call dgetrf(M,M,A,M,IPIV,info)
  if (info .ne. 0) then
    write(*,*) "ERROR - bad exit from dgetrf"
    write(*,*) "info is",info
    stop
  end if

  !I am not sure if it is faster to copy over A(i,i) to 
  !  a working vector and then use SIMD to accelerate the
  !  product accumulation...
  !Accumulate the diagonals of U
  dtmp = A(1,1)
  i = 2
  do while (.true.)
    dtmp = dtmp * A(i,i)
    if (i .ge. M) exit
    i = i + 1
  end do
  !Account for the sign of P
  if (IPIV(1) .ne. 1) dtmp = -1*dtmp
  i = 2
  do while (.true.)
    if (IPIV(i) .ne. i) dtmp = -1*dtmp
    if (i .ge. M) exit
    i = i + 1
  end do

  ddet = dtmp

end function ddet

!-------------------------------------------------------

end module linal
