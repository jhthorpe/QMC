!-------------------------------------------------------
! orb
!	- module containing functions and subrotuines
!	  dealing with the construction and evalutation
!	  of general orbtials
!
! ORBITAL TYPES
! 1. STO  <- only this is implemented, currently
! 2. GTO
!-------------------------------------------------------
module orb_mod
  use sto_mod
  use BUF

!-------------------------------------------------------
! orb_dat
!	- type that stores information about the 
!	  orbitals
!-------------------------------------------------------
type orb_dat
  integer :: Nbas,Ntype
end type orb_dat

contains

!-------------------------------------------------------
! orb_init
!	- initializes the orbitals lists, reading
!	  from the orb.dat file
!
! FORMAT OF ORB.DAT
! line1         : number of total basis functions
! line2         : number of basis function types
! line3         : blank line
! line4         : First typeID
! line5         : number of functions of this type, number of parameters
! line6->6+n-1  : n,l,ml,x,y,z,parameters
! line6+n       : blank line
! repeat until done
!-------------------------------------------------------
! orb		: orb_dat, 
! ICORE		: IBUF, buffer
! DCORE		: DBUF, buffer

subroutine orb_init(orb,sto,ICORE,DCORE)
  implicit none
  type(orb_dat), intent(inout) :: orb
  type(sto_dat), intent(inout) :: sto
  type(IBUF), intent(inout) :: ICORE
  type(DBUF), intent(inout) :: DCORE

  character(len=6) :: otype
  integer :: fid,i,f,N,NP
  logical :: err

  N = 0
  err = .false.
 
  !read in data
  fid = 100
  open(file='orb.dat',unit=fid,status='old')
  read(fid,*) orb%Nbas
  read(fid,*) orb%Ntype
  read(fid,*)

  do i=1,orb%Ntype
    read(fid,*) otype

    !reads for different types
    if (trim(otype).eq.'STO') then
      read(fid,*) N,NP
      write(*,*) "N,NP",N,NP
      call sto_init(N,NP,sto)
      do f=1,N
        read(fid,*) sto%n(f),sto%l(f),sto%ml(f),&
                    sto%xyz(1:3,f),sto%prm(1:NP,f)
      end do        
      read(fid,*) !blank line at end

    end if

  end do 
  close(unit=fid)

  !sort the lists 
  call sto_sort(sto,ICORE,DCORE)
  call sto_print(sto)

end subroutine orb_init
!-------------------------------------------------------

end module orb_mod
!-------------------------------------------------------
