program xvmc
  use orb_mod
  use sto_mod
  use BUF
  implicit none
 
  type(IBUF) :: ICORE
  type(DBUF) :: DCORE
  type(orb_dat) :: orb
  type(sto_dat) :: sto

  !initialize ICORE and DCORE
  call IBUF_INIT(ICORE)
  call DBUF_INIT(DCORE)

  !Read in the orbitals 
  call orb_init(orb,sto,ICORE,DCORE)
  
  call sto_free(sto)
end program xvmc
