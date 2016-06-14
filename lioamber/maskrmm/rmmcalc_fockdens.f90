!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmcalc_fockdens( dens_mao, efield, fock_mao, dipole, energy )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, only: M, Md, RMM

  implicit none
  complex*16, intent(in)    :: dens_mao(:,:)
  real*8,     intent(in)    :: efield(3)
  real*8,     intent(inout) :: fock_mao(:,:)
  real*8,     intent(inout) :: dipole(3)
  real*8,     intent(inout) :: energy

  real*8  :: energy_1e
  real*8  :: energy_coulomb
  real*8  :: energy_xc
  real*8  :: energy_efld

  real*8  :: efield_mod
  real*8  :: gnum, factor
  integer :: MM, MMd, idx0, kk

  call g2g_timer_start('rmmcalc_fockdens')
!
!
! Initializations
!------------------------------------------------------------------------------!
  energy_1e = 0.0d0
  energy_coulomb = 0.0d0
  energy_xc = 0.0d0
  energy_efld = 0.0d0

!  gnum=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
!  gnum=0.0d0
!  factor=(2.54d0*2.00d0)
!  efield_mod=efield(1)**2+efield(2)**2+efield(3)**2
!
!
! Fock calculation
!------------------------------------------------------------------------------!
  call rmmput_dens( dens_mao )
  call int3lu( energy_coulomb )
  call g2g_solve_groups( 0, energy_xc, 0)
!  if ( efield_mod > epsilon(efield_mod) ) 
!    call intfld( gnum, efield(1), efield(2), efield(3) )
!  end if
  call rmmget_fock( fock_mao )
!
!
! Energy calculation
!------------------------------------------------------------------------------!
  MM = M*(M+1)/2
  MMd = Md*(Md+1)/2
  idx0 = 3*MM + 2*MMd
  do kk=1,MM
    energy_1e = energy_1e + RMM(kk)*RMM(idx0+kk)
  enddo

  dipole = 0.0d0
!  call int1
!  energy_efld = energy_efld - 1.00D0 * gnum*(Fx*ux+Fy*uy+Fz*uz)/factor
!  energy_efld = energy_efld - 0.50D0 * (1.0D0-1.0D0/epsilon)*Qc2/a0

  energy = energy + energy_1e
  energy = energy + energy_Coulomb
  energy = energy + energy_xc
!
!
  call g2g_timer_stop('rmmcalc_fockdens')
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
