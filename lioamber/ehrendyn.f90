!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn( energy, dipmom )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use ehrenfest
  use garcha_mod, &
  only: M, rmm, natom, Iz, delta_time_au, total_time_au, atom_mass, &
      & nucpos, nucvel, is_first_step, qm_forces_ds, qm_forces_total, &
      & rhosave_old, rhosave_now, rhosave_new

  implicit none
  real*8, intent(inout)               :: energy
  real*8, intent(inout), dimension(3) :: dipmom

  integer :: ii, jj

  if ( is_first_step ) then
    total_time_au = 0.0d0

    if ( allocated(rhosave_old) ) deallocate(rhosave_old)
    allocate( rhosave_old(M,M) )

    if ( allocated(rhosave_now) ) deallocate(rhosave_now)
    allocate( rhosave_now(M,M) )

    if ( allocated(rhosave_new) ) deallocate(rhosave_new)
    allocate( rhosave_new(M,M) )

    if ( allocated(qm_forces_ds) ) deallocate(qm_forces_ds)
    allocate( qm_forces_ds(3,natom) )
    qm_forces_ds(:,:)=0.0d0

    if ( allocated(qm_forces_total) ) deallocate(qm_forces_total)
    allocate( qm_forces_total(3,natom) )
    qm_forces_total(:,:)=0.0d0

    if (allocated(atom_mass)) deallocate(atom_mass)
    allocate(atom_mass(natom))
    call ehren_masses(natom,Iz,atom_mass)

    call SCF( energy, dipmom )
  end if

  print*, ' IN!'
  print*, ' IN!'
  print*, ' IN!'
  print*, ' IN!'

  write(998,*) ''
  do ii=1,M
  do jj=1,M
    write(998,*) rhosave_old(ii,jj), rhosave_now(ii,jj)
  enddo
  enddo

  call ehrenstep_verlet( rhosave_old, rhosave_now, rhosave_new, qm_forces_ds, &
                       & dipmom, energy, atom_mass, nucpos, nucvel, &
                       & qm_forces_total, delta_time_au, &
                       & is_first_step=is_first_step )

  write(999,*) ''
  do ii=1,M
  do jj=1,M
    write(999,*) rhosave_new(ii,jj)
  enddo
  enddo

  print*, ' OUT!'
  print*, ' OUT!'
  print*, ' OUT!'
  print*, ' OUT!'

  rhosave_old = rhosave_now
  rhosave_now = rhosave_new
  total_time_au = total_time_au + delta_time_au
  is_first_step = .false.


end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
