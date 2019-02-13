!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenfest_main( energy_o )

   implicit none
   real(dpk), intent(inout) :: energy_io
!
!
!
!  Preliminaries to reduce
!------------------------------------------------------------------------------!
   call g2g_timer_start('ehrenfest main')

   timeau   = stored_timeau
   energy   = stored_energy
   dipmom   = stored_dipmom
   forces   = stored_forces
   dens_now = stored_denson
!
!
!
!  Update forces and velocities
!------------------------------------------------------------------------------!
   call recalc_Core( Smat, Hmao, )
   if (forces_recalc) then
      call rmmput_denson
      call recalc_Fock( Hmao, Fmao )
      call calc_forces
   end if

   call force_to_accel( )
   call update_velvec( dt_n )
   call recalc_Dmat( Dmat )
   call recalc_Fock( Hmao, Fmao )
!
!
!
!  Electron dynamics cycle
!------------------------------------------------------------------------------!
   do elecdyn_step = 1,
      call evolve_density( dt_e )
   end do
!
!
!
!  Prepare outputs and save data for next step
!------------------------------------------------------------------------------!
   energy_o = stored_energy
   dipmom_o = stored_dipmom

   stored_timeau =
   stored_energy =
   stored_dipmom =
   stored_forces =
   stored_denson =
!
!
!
!  Deinitialization of variables
!------------------------------------------------------------------------------!

   call g2g_timer_stop('ehrenfest - lio step')

end subroutine ehrendyn_main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!



!  BEWARE OF COUNTING => EXTRA STEP WHEN NOT DOING RESTART...
!   if (first_step) return
   nustep_count = nustep_count + 1
   time = stored_time

   allocate( nucfor_ds(3,natom) )
   allocate( Smat(M,M), Sinv(M,M) )
   allocate( Lmat(M,M), Umat(M,M), Linv(M,M), Uinv(M,M) )
   allocate( Fock(M,M), Fock0(M,M) )
   allocate( RhoOld(M,M), RhoMid(M,M), RhoNew(M,M), RhoMidF(M,M) )
   allocate( Bmat(M,M), Dmat(M,M), Tmat(M,M) )

   dtn = tdstep
   dte = ( tdstep / edyn_steps )

   first_nustep = (nustep_count == 1)
   load_restart = (first_nustep).and.(rsti_loads)
   rhomid_in_ao = (first_nustep).and.(.not.rsti_loads)
   missing_last = (first_nustep).and.(.not.rsti_loads)

   if (first_nustep) stored_energy = energy_o
   if (load_restart) call ehrenaux_rsti( rsti_fname, &
   &  natom, qm_forces_total, nucvel, M, stored_densM1, stored_densM2 )
!
!
!  Update velocities, calculate fixed fock, load last step dens matrices
!------------------------------------------------------------------------------!
   if (velocity_recalc) then
      dtaux = dtn/2.0d0 - dte/2.0d0
      call ehrenaux_updatevel( natom, atom_mass, qm_forces_total, nucvel, dtaux )
   else
      call ehrenaux_updatevel( natom, atom_mass, qm_forces_total, nucvel, dtn )
   endif

   energy0 = 0.0d0
   call RMMcalc0_Init()
   call RMMcalc1_Overlap( Smat, energy0 )
   call ehrenaux_cholesky( M, Smat, Lmat, Umat, Linv, Uinv, Sinv )
   call RMMcalc2_FockMao( Fock0, energy0 )

   RhoOld = stored_densM1
   RhoMid = stored_densM2
   if (rhomid_in_ao) then
      RhoMid = matmul(RhoMid, Lmat)
      RhoMid = matmul(Umat, RhoMid)
      stored_densM2 = RhoMid
   endif
!
!
!
!  ELECTRONIC STEP CYCLE
!------------------------------------------------------------------------------!
   elstep_keeps = ceiling( real(edyn_steps) / 2.0 )

   do elstep_local = 1, edyn_steps
      call g2g_timer_start('ehrendyn - electronic step')
      elstep_count = elstep_count + 1
      dipmom(:) = 0.0d0
      energy = energy0
      Fock = Fock0

      if (velocity_recalc) call ehrenaux_updatevel &
      &  ( natom, atom_mass, qm_forces_total, nucvel, dte )

      call ehrendyn_step( missing_last, propagator, time, dte, M, natom,       &
                        & nucpos, nucvel, nucfor_ds, Sinv, Uinv, Linv,         &
                        & RhoOld, RhoMid, RhoNew, Fock, dipmom, energy )

      RhoOld = RhoMid
      RhoMid = RhoNew

      if ( elstep_local == elstep_keeps ) qm_forces_ds = nucfor_ds
      time = time + dte * 0.0241888d0

      call ehrenaux_writedip(elstep_count, wdip_nfreq, stored_time, dipmom,    &
      &    wdip_fname)

      if (rsto_saves) call ehrenaux_rsto( &
      &  rsto_fname, rsto_nfreq, ndyn_steps*edyn_steps, elstep_count,          &
      & natom, qm_forces_total, nucvel, M, RhoOld, RhoMid )

      call g2g_timer_stop('ehrendyn - electronic step')

   enddo
!
!
!
!  Finalizations
!------------------------------------------------------------------------------!
   stored_densM1 = RhoOld
   stored_densM2 = RhoMid

   dipmom_o = stored_dipmom
   energy_o = stored_energy
   stored_dipmom = dipmom
   stored_energy = energy
   stored_time = time

   deallocate( Smat, Sinv )
   deallocate( Lmat, Umat, Linv, Uinv )
   deallocate( Fock, Fock0 )
   deallocate( RhoOld, RhoMid, RhoNew, RhoMidF )
   deallocate( Bmat, Dmat, Tmat )
   call g2g_timer_stop('ehrenfest main')

end subroutine ehrendyn_main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
