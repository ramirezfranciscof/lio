!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine RMMcalc2_FockMao(DensMao,FockMao,Energy)
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm
  use garcha_mod, only:M,Md,RMM,kkind,kkinds,cool,cools,igrid2
  implicit none
  complex*16,intent(in) :: DensMao(M,M)
  real*8,intent(out)    :: FockMao(M,M)
  real*8,intent(out)    :: Energy

  real*8   :: Energy_1e
  real*8   :: Energy_Coulomb
  real*8   :: Energy_SolvT,Energy_SolvF
  real*8   :: Energy_Exchange
  integer  :: kk,idx0
  integer  :: MM,MMd,igpu
  logical  :: MEMO

!
!
! Calculate fixed-parts of fock
!--------------------------------------------------------------------!
  if (allocated(kkind))  deallocate(kkind)
  if (allocated(kkinds)) deallocate(kkinds)
  if (allocated(cool))   deallocate(cool)
  if (allocated(cools))  deallocate(cools)

  call g2g_reload_atom_positions(igrid2)
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()

  if (igpu.le.1) then
    call intsol(Energy_SolvF,Energy_SolvT,.true.)
  else
    call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
  endif

  call int2()

  if (igpu.gt.2) call aint_coulomb_init()
  if (igpu.eq.5) MEMO = .false.
  if (MEMO) then
    call g2g_timer_start('int3mem')
    call int3mem()
    call g2g_timer_stop('int3mem')
  endif

!
! Calculate unfixed Fock in RMM
!--------------------------------------------------------------------!
  call rmmput_dens(DensMao)
  call g2g_timer_start('g2g-solve + int3lu')
  call int3lu(Energy_Coulomb)
  call g2g_solve_groups(0,Energy_Exchange,0)
  call g2g_timer_stop('g2g-solve + int3lu')

  MM=M*(M+1)/2
  MMd=Md*(Md+1)/2
  idx0=3*MM+2*MMd
  Energy_1e=0.0d0
  do kk=1,MM
    Energy_1e=Energy_1e+RMM(kk)*RMM(idx0+kk)
  enddo

!  Energy=0.0d0
  Energy=Energy+Energy_1e
  Energy=Energy+Energy_Coulomb
  Energy=Energy+Energy_SolvT
  Energy=Energy+Energy_Exchange

!
! Extract FockMao from RMM
!--------------------------------------------------------------------!
  call rmmget_fock(FockMao)


  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
