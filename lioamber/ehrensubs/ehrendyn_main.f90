!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn_main( energy_o, dipmom_o )
!------------------------------------------------------------------------------!
!
!  stored_densM1 and stored_densM2 are stored in ON basis, except for the
!  first step.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, &
   &  only: M, natom, time_step, atom_mass, nucpos, nucvel                     &
   &      , qm_forces_ds, qm_forces_total

   use lionml_data, &
   &  only: ndyn_steps, edyn_steps, edyn_prop                                  &
   &      , rsti_loads, rsti_fname, rsto_saves, rsto_nfreq, rsto_fname

   use ehrendata, &
   &  only: stored_time, stored_energy, stored_dipmom, stored_forces           &
   &      , stored_densM1, stored_densM2                                       &
   &      , rsti_funit, rsto_funit, nustep_count

   implicit none
   real*8,intent(inout) :: dipmom_o(3), energy_o
   real*8               :: dipmom(3)  , energy  , energy0

   real*8  :: time, dtn, dte, dtx
   integer :: elstep_local, elstep_total, elstep_keeps
   integer :: substep, substeps
   integer :: nn, kk

   logical :: first_nustep
   logical :: load_restart
   logical :: rhomid_in_ao
   logical :: missing_last

   real*8, allocatable, dimension(:,:) :: nucfor_ds
   real*8, allocatable, dimension(:,:) :: Smat, Sinv
   real*8, allocatable, dimension(:,:) :: Lmat, Umat, Linv, Uinv
   real*8, allocatable, dimension(:,:) :: Fock, Fock0
   real*8, allocatable, dimension(:,:) :: Bmat, Dmat

   complex*16, allocatable, dimension(:,:) :: RhoOld, RhoMid, RhoNew
   complex*16, allocatable, dimension(:,:) :: RhoMidF
   complex*16, allocatable, dimension(:,:) :: Tmat
!
!
!
!  INITIALIZATIONS AND LOGICALS ASSIGNMENT
!------------------------------------------------------------------------------!
   allocate( nucfor_ds(3,natom) )
   allocate( Smat(M,M), Sinv(M,M) )
   allocate( Lmat(M,M), Umat(M,M), Linv(M,M), Uinv(M,M) )
   allocate( Fock(M,M), Fock0(M,M) )
   allocate( Bmat(M,M), Dmat(M,M), Tmat(M,M) )
   allocate( RhoOld(M,M), RhoMid(M,M), RhoNew(M,M), RhoMidF(M,M) )

   dtn = time_step
   dte = (time_step / edyn_steps)

   if ( (nustep_count == 0).and.(rsti_load) ) then
      load_restar    = .true.
      backwards_step = .false.
      substep_dt     = dte
      substep_totals = 1

   else if ( (nustep_count == 0).and.(.not.rsti_load) ) then
      load_restar    = .false.
      backwards_step = .true.
      substep_dt     = -(dte / 40)
      substep_totals = 1

   else if ( (nustep_count == 1).and.(.not.rsti_load) ) then
      load_restar    = .false.
      backwards_step = .false.
      substep_dt     = (dte / 20)
      substep_totals = 20

   else
      load_restar    = .false.
      backwards_step = .false.
      substep_dt     = dte
      substep_totals = 1

   endif
!
!
!
!  LOAD RESTART OR CALC SCF
!------------------------------------------------------------------------------!
   if (load_restar) then
      print*,'RESTARTS NOT AVAILABLE YET'; stop
!      call ehrenrsti_load( rsti_fname, rsti_funit, natom, qm_forces_total,  &
!                         & nucvel, M, stored_densM1, stored_densM2 )
   else
      call SCF( energy_o, dipmom_o )
      dens_in_ao = .true.

      stored_energy = energy_o
      stored_dipmom = dipmom_o
      stored_time   = 0.0d0
   endif
!
!
!
!  PREPARATIONS
!------------------------------------------------------------------------------!
   energy0 = 0.0d0
   call RMMcalc0_Init()
   call RMMcalc1_Overlap( Smat,  energy0 )
   call RMMcalc2_FockMao( Fock0, energy0 )

   call ehren_cholesky( M, Smat, Lmat, Umat, Linv, Uinv, Sinv )

   RhoOld = stored_densM1
   RhoMid = stored_densM2
   if (dens_in_ao) then
      RhoMid = matmul(RhoMid, Lmat)
      RhoMid = matmul(Umat, RhoMid)
      RhoOld = RhoMid
      stored_densM2 = RhoMid
   endif

   dtx = (0.5d0) * (dtn)
   call ehrenaux_updatevels( dtx, natom, atom_mass, qm_forces_total, nucvel )
!
!
!
!  INTERNAL ELECTRONIC PROPAGATION CYCLE
!------------------------------------------------------------------------------!
   elstep_keeps = ceiling( real(elstep_total) / 2.0d0 )

   do elstep_local = 1, elstep_total

      dipmom(:) = 0.0d0
      energy = energy0
      Fock   = Fock0

      call ehrendyn_prep &
      &  ( )

      call ehrendyn_step &
      &  ( elec_prop, time, dtaux, M, natom, nucpos, nucvel, nucfor_ds,        &
      &    Sinv, Uinv, Linv, RhoOld, RhoMid, RhoNew, Fock, energy, dipmom )

      if ( elstep_local == elstep_keeps ) qm_forces_ds = nucfor_ds
      RhoOld = RhoMid
      RhoMid = RhoNew

      dtx = (0.5d0) * (elecstep_dt)
      call ehrenaux_updatevels( dtx, natom, atom_mass, qm_forces_total, nucvel )

   enddo

   stored_densM1 = RhoOld
   stored_densM2 = RhoMid
!
!
!
!  Finalizations
!------------------------------------------------------------------------------!
   dipmom_o = stored_dipmom
   energy_o = stored_energy

   stored_dipmom = dipmom
   stored_energy = energy
   stored_time   = time

   ehrenout_dipole( dipmom, 134, write_header )

   if (rsto_saves) then
      print*,'RESTARTS NOT AVAILABLE YET'; stop
!      call ehrenrsto_save( rsto_fname, rsto_funit, rsto_nfreq, ndyn_steps,     &
!         & nustep_count, Natom, qm_forces_total, nucvel,                       &
!         & M, stored_densM1, stored_densM2)
   endif

   time = time + dtn * 0.0241888d0
   nustep_count = nustep_count + 1
   deallocate( nucfor_ds )
   deallocate( Smat, Sinv )
   deallocate( Lmat, Umat, Linv, Uinv )
   deallocate( Fock, Fock0 )
   deallocate( Bmat, Dmat, Tmat )
   deallocate( RhoOld, RhoMid, RhoNew, RhoMidF )
end subroutine ehrendyn_main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
