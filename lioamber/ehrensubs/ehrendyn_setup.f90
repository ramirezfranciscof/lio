!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn_setup( natoms, nbasis, realrho )
   use ehrendata,  only: stored_densM1, stored_densM2
   use garcha_mod, only: qm_forces_ds, qm_forces_total

   implicit none
   integer, intent(in) :: natoms
   integer, intent(in) :: nbasis
   real*8,  intent(in) :: realrho( nbasis, Nbasis )

   natoms = natoms
   nbasis = nbasis

   if (allocated(stored_densM1)) deallocate(stored_densM1)
   allocate(stored_densM1( nbasis, nbasis ))
   stored_densM1 = DCMPLX( 0.0d0 )

   if (allocated(stored_densM2)) deallocate(stored_densM2)
   allocate(stored_densM2( nbasis, nbasis ))
   stored_densM2 = DCMPLX( realrho )

   if (allocated(qm_forces_total)) deallocate(qm_forces_total)
   allocate( qm_forces_total(3, natoms) )
   qm_forces_total = 0.0d0

   if (allocated(qm_forces_ds)) deallocate(qm_forces_ds)
   allocate( qm_forces_ds(3, natoms) )
   qm_forces_ds = 0.0d0

end subroutine ehrendyn_setup
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
