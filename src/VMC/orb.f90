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
  use linal

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
! orb_build_SD
!	- given a set of orbitals evaluated at all
!	  the electron's positions, constructs a SD
!	  matrix from a given set of MOCs (Molecular
!	  orbital coefficients)  
!
!	- for efficient use, MOs,bas,M,CF, and SD should
!	  be stored in DCORE in the calling 
!	  routine
!-------------------------------------------------------
! Nelc		: int, number of electrons
! Nmos		: int, number of orbitals
! Nbas		: 2D real*8, basis functions @ all R
! MOs		: 2D real*8, MO coefs
! bas		: 2D real*8, basis functions
! Phi		: 2D real*8, matrix of orbitals 
! CF		: 2D real*8, cofactor matrix of Phi, occ only
! SD		: real*8, Slater-Det of Phi
! ICORE		: IBUF 
! DCORE		: DBUF

subroutine orb_build_SD(Nelc,Nmos,Nbas,MOs,bas,Phi,CF,SD,ICORE,DCORE)
  implicit none
  integer, intent(in) :: Nelc,Nmos,Nbas
  real(kind=8), intent(in) :: MOs(Nbas,Nmos),bas(Nbas,Nelc)
  real(kind=8), intent(inout) :: Phi(Nmos,Nelc),CF(Nelc,Nelc)
  real(kind=8), intent(inout) :: SD
  type(IBUF), intent(inout) :: ICORE
  type(DBUF), intent(inout) :: DCORE

  integer :: i0,i1,LWORK,info

  !Construct Phi
  call dgemm('T','N',Nmos,Nbas,Nelc,1.d0,MOs,Nbas,bas,Nbas,&
              0.d0,Phi,Nmos)

  !LU factorize Phi to get IPIV and determinant
  CF(1:Nelc,1:Nelc) = transpose(Phi(1:Nelc,1:Nelc))
  call IBUF_RESERVE(i0,Nelc,ICORE)
  SD = ddet(Nelc,CF,ICORE%buf(i0:i0+Nelc-1))

  !Invert first NxN part of Phi to get cofactors
  !Get LWORK
  call dgetri(Nelc,CF,Nelc,ICORE%buf(i0:i0+Nelc-1),&
              DCORE%buf(DCORE%next),-1,info)
  LWORK = CEILING(DCORE%buf(DCORE%next))
  call DBUF_RESERVE(i1,LWORK,DCORE)

  !actually get inverse of transpose at last
  call dgetri(Nelc,CF,Nelc,ICORE%buf(i0:i0+Nelc-1),&
              DCORE%buf(i1:i1+LWORK-1),info)
  if (info .ne. 0) then
    write(*,*) "ERROR - problem out of dgetri in orb_build_SD"
    stop
  end if

  call DBUF_POP(LWORK,DCORE)
  call IBUF_POP(Nelc,ICORE)

end subroutine orb_build_SD

!-------------------------------------------------------

end module orb_mod
!-------------------------------------------------------
