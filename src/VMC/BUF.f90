!-------------------------------------------------------
! BUF.f90
!	- module containing subroutines for dealing with
!	  the core buffers
!-------------------------------------------------------
module buf

!-------------------------------------------------------
! DBUF
!	- double prescision buffer
!-------------------------------------------------------
type DBUF
  real(kind=8), dimension(1000000) :: buf
  integer :: nelm,nfree,next
end type DBUF

!-------------------------------------------------------
! IBUF
!	-integer buffer
!-------------------------------------------------------
type IBUF
  integer, dimension(1000000) :: buf
  integer :: nelm,nfree,next
end type IBUF
!-------------------------------------------------------

contains

!-------------------------------------------------------
! DBUF_INIT
!	- initializes the double buffer
!-------------------------------------------------------
subroutine DBUF_INIT(DCORE)
  implicit none
  type(DBUF), intent(inout) :: DCORE
  DCORE%next = 1
  DCORE%nelm = 1000000
  DCORE%nfree = DCORE%nelm
  call DBUF_ZERO(DCORE)
end subroutine DBUF_INIT

!-------------------------------------------------------
! IBUF_INIT
!	- initializes the int buffer
!-------------------------------------------------------
subroutine IBUF_INIT(ICORE)
  implicit none
  type(IBUF), intent(inout) :: ICORE
  ICORE%next = 1
  ICORE%nelm = 1000000
  ICORE%nfree = ICORE%nelm
  call IBUF_ZERO(ICORE)
end subroutine IBUF_INIT

!-------------------------------------------------------
! DBUF_ZERO
!	- zeros the double buffer 
!-------------------------------------------------------
subroutine DBUF_ZERO(DCORE)
  implicit none
  type(DBUF), intent(inout) :: DCORE
  DCORE%buf(1:DCORE%nelm) = 0.d0 
end subroutine DBUF_ZERO

!-------------------------------------------------------
! IBUF_ZERO
!	- zeros the double buffer 
!-------------------------------------------------------
subroutine IBUF_ZERO(ICORE)
  implicit none
  type(IBUF), intent(inout) :: ICORE
  ICORE%buf(1:ICORE%nelm) = 0
end subroutine IBUF_ZERO

!-------------------------------------------------------
! DBUF_RESERVE
!	- checks out a section of DBUF 
!-------------------------------------------------------
! i		: int, starting index
! N		: int, length of 
subroutine DBUF_RESERVE(i,N,DCORE)
  implicit none
  integer, intent(in) :: N
  integer, intent(inout) :: i
  type(DBUF), intent(inout) :: DCORE
  if (N .le. DCORE%nfree) then
    i = DCORE%next
    DCORE%next = i + N 
    DCORE%nfree = DCORE%nfree - N
  else
    write(*,*) "ERROR - out of memory in DBUF"
    stop
  end if
end subroutine DBUF_RESERVE

!-------------------------------------------------------
! IBUF_RESERVE
!	- reserves a section of IBUF
!-------------------------------------------------------
! i		: int, starting index
! N		: int, length of 
subroutine IBUF_RESERVE(i,N,ICORE)
  implicit none
  integer, intent(in) :: N
  integer, intent(inout) :: i
  type(IBUF), intent(inout) :: ICORE
  if (N .le. ICORE%nfree) then
    i = ICORE%next
    ICORE%next = i + N 
    ICORE%nfree = ICORE%nfree - N
  else
    write(*,*) "ERROR - out of memory in IBUF"
    stop
  end if
end subroutine IBUF_RESERVE

!-------------------------------------------------------
! DBUF_POP
!	- removes the top N elements of DBUF 
!-------------------------------------------------------
! N		: int, length of 
subroutine DBUF_POP(N,DCORE)
  implicit none
  integer, intent(in) :: N
  type(DBUF), intent(inout) :: DCORE
  integer :: M
  M = min(N,DCORE%nelm - DCORE%nfree)
  DCORE%nfree = DCORE%nfree + M
  DCORE%next = DCORE%next - M 
end subroutine DBUF_POP

!-------------------------------------------------------
! IBUF_POP
!	- removes the top N elements of IBUF 
!-------------------------------------------------------
! N		: int, length of 
subroutine IBUF_POP(N,ICORE)
  implicit none
  integer, intent(in) :: N
  type(IBUF), intent(inout) :: ICORE
  integer :: M
  M = min(N,ICORE%nelm - ICORE%nfree)
  ICORE%nfree = ICORE%nfree + M
  ICORE%next = ICORE%next - M 
end subroutine IBUF_POP

!-------------------------------------------------------
! DBUF_PRINT
!	- prints information about DBUF
!-------------------------------------------------------
subroutine DBUF_PRINT(DCORE)
  implicit none
  type(DBUF), intent(inout) :: DCORE
  write(*,*) "DBUF : nelm"   ,DCORE%nelm
  write(*,*) "DBUF : nfree"  ,DCORE%nfree
  write(*,*) "DBUF : next"   ,DCORE%next
end subroutine DBUF_PRINT

!-------------------------------------------------------
! IBUF_PRINT
!	- prints information about DBUF
!-------------------------------------------------------
subroutine IBUF_PRINT(ICORE)
  implicit none
  type(IBUF), intent(inout) :: ICORE
  write(*,*) "IBUF : nelm"   ,ICORE%nelm
  write(*,*) "IBUF : nfree"  ,ICORE%nfree
  write(*,*) "IBUF : next"   ,ICORE%next
end subroutine IBUF_PRINT

end module buf
