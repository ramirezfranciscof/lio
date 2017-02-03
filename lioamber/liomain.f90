!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIOMAIN.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the liomain subroutine, which performs several common     !
! procedures before and after calling SCF either by LIO alone or when          !
! performing in tantem with AMBER/GROMACS. Currently the routine:              !
! * Allocates matrices needed for SCF calculations.                            !
! * Calls SCF or SCFOP for closed/open shell calculations.                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 
subroutine liomain(E, dipxyz)
    use garcha_mod, only : M, Smat, RealRho, OPEN, writeforces, energy_freq,   &
                           restart_freq, npas, sqsm, mulliken, lowdin, dipole, &
                           do_ehrenfest, first_step, nucpos, &
                           Eorbs, fukui, print_coeffs
    use ecp_mod   , only : ecpmode, IzECP
    use ehrenfest
 
    implicit none
    REAL*8, intent(inout) :: dipxyz(3), E

    if (.not.allocated(Smat))    allocate(Smat(M,M))
    if (.not.allocated(RealRho)) allocate(RealRho(M,M))
    if (.not.allocated(sqsm))    allocate(sqsm(M,M))
    if (.not.allocated(Eorbs))   allocate(Eorbs(M))

!------------------------------------------------------------------------------!
! FFR - Ehrenfest
    if ( do_ehrenfest ) then
        if ( first_step ) call SCF( E, dipxyz )
        call ehrendyn( E, dipxyz )
    else
!------------------------------------------------------------------------------!
! FFR - Normal Code
    if(OPEN) then
        if (ecpmode) stop "ECP is unavailable for Open Shell systems."
        call SCFOP(E, dipxyz)
    else
        call SCF(E)
    endif
!------------------------------------------------------------------------------! 
! FFR - Ehrenfest
    endif
!------------------------------------------------------------------------------! 

    ! Perform Mulliken and Lowdin analysis, get fukui functions and dipole.
    if (MOD(npas, energy_freq).eq.0) then

        if (mulliken.or.lowdin) call do_population_analysis()

        if (dipole) call do_dipole(dipxyz, 69)
  
        if (fukui) call do_fukui()

        if (writeforces) then
            if (ecpmode) stop "ECP does not feature forces calculation."
            call do_forces(123)
        endif

        if (print_coeffs) call write_orbitals(29)
    endif

    if ((restart_freq.gt.0).and.(MOD(npas, restart_freq).eq.0)) call write_restart(88)

    return
end subroutine liomain
