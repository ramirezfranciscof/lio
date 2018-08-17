!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!  This subroutine calculates the part of the fock matrix that depends on
!  the electronic density, but it must receive as input the core fock.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_fock_xx( dipole_xyz, energy_coul, energy_xc                 &
                         &, apply_efield, efield_xyz, energy_field )
!------------------------------------------------------------------------------!

   use garcha_mod  , only: M, natom, Iz, NCO, Nunp, Dbug
   use field_data  , only: a0, epsilon
   use faint_cpu77 , only: int3lu, intfld

   implicit none
   real*8    , intent(out)   :: dipole_xyz(3)
   real*8    , intent(out)   :: energy_coul
   real*8    , intent(out)   :: energy_xc

   logical   , intent(in)    :: apply_efield
   real*8    , intent(in)    :: efield_xyz(3)
   real*8    , intent(inout) :: energy_field

   integer :: MM, kk
   real*8  :: factor, g, Qc
   real*8  :: dip_times_field, strange_term
!
!
!  Calculate unfixed Fock in RMM - int3lu and solve_groups
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc3-solve3lu')
   call rmmCheckNaNs( "Start" )
   call int3lu( energy_coul )
   call rmmCheckNaNs( "Post Coulomb" )
   print*, "solve_groups i"
   call rmmgen_Write( 666, 0 )
   call g2g_solve_groups( 0, energy_xc, 0 )
   print*, "solve_groups o"
   call rmmCheckNaNs( "Post Ex-Corr" )
   call g2g_timer_stop('rmmcalc3-solve3lu')
!
!
!  Calculate unfixed Fock in RMM - electric field
!------------------------------------------------------------------------------!
   call g2g_timer_start('rmmcalc3-dipole')
   call dip( dipole_xyz(1), dipole_xyz(2), dipole_xyz(3) )
   call g2g_timer_stop('rmmcalc3-dipole')

   if (apply_efield) then
      call g2g_timer_start('rmmcalc3-field')
      g = 1.0d0
      factor = 2.54d0

      Qc = (-2.0d0) * NCO + Nunp
      do kk = 1, natom
         Qc = Qc + Iz(kk)
      end do

      call intfld( g, efield_xyz(1), efield_xyz(2), efield_xyz(3) )

      dip_times_field = 0.0d0
      dip_times_field = dip_times_field + efield_xyz(1) * dipole_xyz(1)
      dip_times_field = dip_times_field + efield_xyz(2) * dipole_xyz(2)
      dip_times_field = dip_times_field + efield_xyz(3) * dipole_xyz(3)
      strange_term = (0.5d0) * (1.0d0 - 1.0d0/epsilon) * Qc**2 / a0

      energy_field = 0.0d0
      energy_field = energy_field - g * dip_times_field / factor
      energy_field = energy_field - strange_term

      call g2g_timer_stop('rmmcalc3-field')
   endif

end subroutine rmmCalc_fock_xx
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_fock_cs_d( dens_mao, core_mao                               &
                           &, fock_mao, dipole_xyz                             &
                           &, energy_coul, energy_xc                           &
                           &, apply_efield, efield_xyz, energy_field )
!------------------------------------------------------------------------------!
   use garcha_mod  , only: M
   implicit none
   real*8    , intent(in)    :: dens_mao(M,M)
   real*8    , intent(in)    :: core_mao(M,M)
   real*8    , intent(out)   :: fock_mao(M,M)
   real*8    , intent(out)   :: dipole_xyz(3)
   real*8    , intent(out)   :: energy_coul
   real*8    , intent(out)   :: energy_xc

   logical   , intent(in)    :: apply_efield
   real*8    , intent(in)    :: efield_xyz(3)
   real*8    , intent(inout) :: energy_field

   call rmmput_core( core_mao )
   call rmmput_dens( dens_mao )
   call rmmCalc_fock_xx( dipole_xyz, energy_coul, energy_xc     &
                      &, apply_efield, efield_xyz, energy_field )
   call rmmget_fock( fock_mao )

end subroutine rmmCalc_fock_cs_d
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_fock_cs_z( dens_mao, core_mao                               &
                           &, fock_mao, dipole_xyz                             &
                           &, energy_coul, energy_xc                           &
                           &, apply_efield, efield_xyz, energy_field )
!------------------------------------------------------------------------------!
   use garcha_mod  , only: M
   implicit none
   complex*16, intent(in)    :: dens_mao(M,M)
   real*8    , intent(in)    :: core_mao(M,M)
   real*8    , intent(out)   :: fock_mao(M,M)
   real*8    , intent(out)   :: dipole_xyz(3)
   real*8    , intent(out)   :: energy_coul
   real*8    , intent(out)   :: energy_xc

   logical   , intent(in)    :: apply_efield
   real*8    , intent(in)    :: efield_xyz(3)
   real*8    , intent(inout) :: energy_field

   call rmmput_core( core_mao )
   call rmmput_dens( dens_mao )
   call rmmCalc_fock_xx( dipole_xyz, energy_coul, energy_xc     &
                      &, apply_efield, efield_xyz, energy_field )
   call rmmget_fock( fock_mao )

end subroutine rmmCalc_fock_cs_z
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_fock_os_d( densa_mao, densb_mao, core_mao                   &
                           &, focka_mao, fockb_mao, dipole_xyz                 &
                           &, energy_coul, energy_xc                           &
                           &, apply_efield, efield_xyz, energy_field )
!------------------------------------------------------------------------------!
   use garcha_mod  , only: M
   implicit none
   real*8    , intent(in)    :: densa_mao(M,M)
   real*8    , intent(in)    :: densb_mao(M,M)
   real*8    , intent(in)    :: core_mao(M,M)
   real*8    , intent(out)   :: focka_mao(M,M)
   real*8    , intent(out)   :: fockb_mao(M,M)
   real*8    , intent(out)   :: dipole_xyz(3)
   real*8    , intent(out)   :: energy_coul
   real*8    , intent(out)   :: energy_xc

   logical   , intent(in)    :: apply_efield
   real*8    , intent(in)    :: efield_xyz(3)
   real*8    , intent(inout) :: energy_field

   print*, "inside rmmCalc_fock_os_d"
   call rmmput_core( core_mao )
   call rmmput_dens( densa_mao, densb_mao )
   call rmmCalc_fock_xx( dipole_xyz, energy_coul, energy_xc     &
                      &, apply_efield, efield_xyz, energy_field )
   call rmmget_fock( focka_mao, fockb_mao )
   print*, "ending rmmCalc_fock_os_d"

end subroutine rmmCalc_fock_os_d
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCalc_fock_os_z( densa_mao, densb_mao, core_mao                   &
                           &, focka_mao, fockb_mao, dipole_xyz                 &
                           &, energy_coul, energy_xc                           &
                           &, apply_efield, efield_xyz, energy_field )
!------------------------------------------------------------------------------!
   use garcha_mod  , only: M
   implicit none
   complex*16, intent(in)    :: densa_mao(M,M)
   complex*16, intent(in)    :: densb_mao(M,M)
   real*8    , intent(in)    :: core_mao(M,M)
   real*8    , intent(out)   :: focka_mao(M,M)
   real*8    , intent(out)   :: fockb_mao(M,M)
   real*8    , intent(out)   :: dipole_xyz(3)
   real*8    , intent(out)   :: energy_coul
   real*8    , intent(out)   :: energy_xc

   logical   , intent(in)    :: apply_efield
   real*8    , intent(in)    :: efield_xyz(3)
   real*8    , intent(inout) :: energy_field

   print*, "in rmmCalc_fock_os_z"
   call rmmput_core( core_mao )
   call rmmput_dens( densa_mao, densb_mao )
   call rmmCalc_fock_xx( dipole_xyz, energy_coul, energy_xc     &
                      &, apply_efield, efield_xyz, energy_field )
   call rmmget_fock( focka_mao, fockb_mao )

end subroutine rmmCalc_fock_os_z
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
