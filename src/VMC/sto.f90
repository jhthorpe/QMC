!-------------------------------------------------------
! sto_mod
!	- contains subroutines for STO specific orbitals
!
!  ORBITAL TYPE LIST
!  1. 1_S_0  <- only this is implemented right now 
!  2. 2_S_0
!  3. 2_P_x
!  4. 2_P_y
!  5. 2_P_z
!  ibgn stores starting point of each type  
!  ilen stores length of each type 
!
!  TODO
!  1. This could probably be modifed to use ICORE
!     and DCORE to store the information more efficiently
!
!-------------------------------------------------------
module sto_mod
  use constants

!-------------------------------------------------------
! sto_dat
!	- type that stores information about the STO
!	  orbitals 
! prm		: 1D real*8, list of STO prm
! ibgn          : 1D int, starting point of each orb. type
! ilen          : 1D int, length of each orb. type
! sto_n         : 1D int, list of n (principle) quantum numbers
! sto_l         : 1D int, list of l (angular) quantum numbers
! sto_ml        : 1D int, list of ml quantum numbers
!-------------------------------------------------------
type sto_dat
  real(kind=8), allocatable :: xyz(:,:)
  real(kind=8), allocatable :: prm(:,:)
  integer, allocatable :: ibgn(:),ilen(:) 
  integer, allocatable :: n(:),l(:),ml(:)
  integer :: Nsto,Nprm
end type sto_dat
!-------------------------------------------------------

contains

!-------------------------------------------------------
! sto_init
!	- initializes the sto data type
!-------------------------------------------------------
subroutine sto_init(N,NP,sto)
  implicit none
  integer, intent(in) :: N,NP
  type(sto_dat), intent(inout) :: sto
  allocate(sto%xyz(1:3,1:N))  
  allocate(sto%prm(1:NP,1:N))
  allocate(sto%ibgn(1:N))
  allocate(sto%ilen(1:N))
  allocate(sto%n(1:N))
  allocate(sto%l(1:N))
  allocate(sto%ml(1:N))
  sto%Nsto = N
  sto%Nprm = NP
end subroutine sto_init


!-------------------------------------------------------
! sto_free
!	- frees the sto data 
!-------------------------------------------------------
subroutine sto_free(sto)
  implicit none
  type(sto_dat), intent(inout) :: sto
  if (allocated(sto%ml)) deallocate(sto%ml)
  if (allocated(sto%l)) deallocate(sto%l)
  if (allocated(sto%n)) deallocate(sto%n)
  if (allocated(sto%ilen)) deallocate(sto%ilen)
  if (allocated(sto%ibgn)) deallocate(sto%ibgn)
  if (allocated(sto%prm)) deallocate(sto%prm)
  if (allocated(sto%xyz)) deallocate(sto%xyz)
end subroutine sto_free


!-------------------------------------------------------
! sto_print
!	- prints the sto data 
!-------------------------------------------------------
subroutine sto_print(sto)
  implicit none
  type(sto_dat), intent(in) :: sto
  integer :: i,f

99 format((2x,A))
98 format((2x,A,1x,I3))
97 format((2x,A,3(1x,F10.7)))

  write(*,99) "Slater-Type Orbitals"
  write(*,98) "Number of STO        :",sto%Nsto
  write(*,98) "Number of Parameters :",sto%Nprm

  !Type 1 integrals...
  write(*,99) '1S orbitals' 
  do f=sto%ibgn(1),sto%ibgn(1)+sto%ilen(1)-1
    write(*,98) "Orbital #",f
    write(*,97) "Orbital Center @ ",sto%xyz(1,f),&
                                    sto%xyz(2,f),&
                                    sto%xyz(3,f)
    write(*,'(2x,A)',advance='no') "Orbital Parameters"
    do i=1,sto%Nprm-1
      write(*,'(1x,F10.7)',advance='no') sto%prm(i,f) 
    end do  
    write(*,'(1x,F10.7)') sto%prm(i,f) 
      
  end do
  write(*,*) 
  
end subroutine sto_print
!-------------------------------------------------------

!-------------------------------------------------------
! sto_sphr_eval_all
!	- evaluates all the STO type orbitals
!	  at the given electron positions
!-------------------------------------------------------
! Nelc		: int, number of electrons
! Nsto		: int, number of STOs
! Rfi           : 2D real*8, list of f-e radii
! Tfi		: 2D real*8, list of f-e theta
! Pfi		: 2D real*8, list of f-e phi 
! fnc		: 2D real*8, set of functions (fnc,elc) 

subroutine sto_sphr_eval_all(Nelc,Nsto,sto,Rfi,Tfi,Pfi,fnc)
  implicit none
  integer, intent(in) :: Nelc,Nsto
  type(sto_dat), intent(in) :: sto
  real(kind=8), dimension(Nsto,Nelc), intent(in) :: Rfi,Tfi,Pfi
  real(kind=8), intent(inout) :: fnc(Nsto,Nelc)

  integer :: i,f

  !-----------------------------------------
  ! TYPE 1 : 1_S_0
  ! R(r)     = 2.d0*exp(-prm*r)
  ! P(theta) = 1/sqrt(2)
  ! F(phi)   = 1/sqrt(pi*2)
  do i=1,Nelc
    do f=sto%ibgn(1),sto%ibgn(1)+sto%ilen(1)-1
     fnc(f,i) = exp(-1.d0*sto%prm(1,f)*Rfi(f,i))*isqrtpi 
    end do 
  end do 

end subroutine sto_sphr_eval_all

!-------------------------------------------------------
! sto_sphr_eval_1e
!	- evaluates all the STO type orbitals
!	  for a SINGLE electron
!-------------------------------------------------------
! Nsto		: int, number of STOs
! Rfi           : 1D real*8, list of f-e radii
! Tfi		: 1D real*8, list of f-e theta
! Pfi		: 1D real*8, list of f-e phi 
! fnc		: 1D real*8, set of functions (fnc,elc) 

subroutine sto_sphr_eval_1e(Nsto,sto,Rfi,Tfi,Pfi,fnc)
  implicit none
  integer, intent(in) :: Nsto
  type(sto_dat), intent(in) :: sto
  real(kind=8), dimension(Nsto), intent(in) :: Rfi,Tfi,Pfi
  real(kind=8), intent(inout) :: fnc(Nsto)

  integer :: f

  !-----------------------------------------
  ! TYPE 1 : 1_S_0
  ! R(r)     = 2.d0*exp(-prm*r)
  ! P(theta) = 1/sqrt(2)
  ! F(phi)   = 1/sqrt(pi*2)
  do f=sto%ibgn(1),sto%ibgn(1)+sto%ilen(1)-1
    fnc(f) = exp(-1.d0*sto%prm(1,f)*Rfi(f))*isqrtpi 
  end do 

end subroutine sto_sphr_eval_1e
!-------------------------------------------------------

end module sto_mod
!-------------------------------------------------------
