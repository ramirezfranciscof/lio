!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrentest( energy, dipmom )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use ehrenfest
  use garcha_mod, &
  only: M, rmm, natom, Iz, delta_time_au, total_time_au, atom_mass, &
      & nucpos, nucvel, is_first_step, qm_forces_ds, qm_forces_total, &
      & rhosave_old, rhosave_now, rhosave_new

  implicit none
  real*8, intent(inout)               :: energy
  real*8, intent(inout), dimension(3) :: dipmom

  if ( is_first_step ) then
    if ( allocated(rhosave_now) ) deallocate(rhosave_now)
    allocate( rhosave_now(M,M) )

    if ( allocated(qm_forces_ds) ) deallocate(qm_forces_ds)
    allocate( qm_forces_ds(3,natom) )
    qm_forces_ds(:,:)=0.0d0

    if ( allocated(qm_forces_total) ) deallocate(qm_forces_total)
    allocate( qm_forces_total(3,natom) )
    qm_forces_total(:,:)=0.0d0
  end if
  call SCF( energy, dipmom )

  if ( is_first_step ) then
!    call int1(En)
!    call spunpack('L',M,RMM(M5),Smat)
!    call spunpack('L',M,RMM(M1),RealRho)
!    call fixrho(M,RealRho)
!    rhosave_now=DCMPLX(RealRho)
  end if


end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
