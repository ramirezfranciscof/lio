!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn( Energy_o, DipMom_o )
!------------------------------------------------------------------------------!
!
! RhoSaveA and RhoSaveB are stored in ON basis, except for the first step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, &
   &   only: M, natom, total_time, first_step, atom_mass               &
   &       , nucpos, nucvel, qm_forces_ds, qm_forces_total

   use liokeys, &
   &   only: ndyn_steps, verbose_levels                                        &
   &       , rstinp, rstinp_fname, rstout, rstout_nfreq, rstout_fname

   use ehrendata, &
   &   only: stored_energy, step_number, rstinp_funit, rstout_funit            &
   &       , dt_elec, dt_nucl                                                  &
   &       , RhoSaveA, RhoSaveB

   implicit none
   real*8,intent(inout) :: Energy_o, DipMom_o(3)

   real*8 :: Energy, DipMom(3)

   real*8,allocatable,dimension(:,:)     :: Smat,Sinv,Lmat,Umat,Linv,Uinv
   real*8,allocatable,dimension(:,:)     :: Fock,FockInt
   complex*16,allocatable,dimension(:,:) :: RhoOld,RhoMid,RhoNew

   real*8,allocatable,dimension(:,:)     :: Bmat,Dmat
   complex*16,allocatable,dimension(:,:) :: Tmat
   real*8                                :: dipxyz(3), dipole_norm
   integer                               :: ii,jj,kk,idx
   logical                               :: save_this_step

! Preliminaries
!------------------------------------------------------------------------------!
   call g2g_timer_start('ehrendyn call')
   print*,'holi'
   if ( verbose_levels > 0 ) then
      print*,'Doing ehrenfest!'
   endif

   step_number = step_number + 1

   allocate( Smat(M,M), Sinv(M,M) )
   allocate( Lmat(M,M), Umat(M,M), Linv(M,M), Uinv(M,M) )
   allocate( Fock(M,M), FockInt(M,M) )
   allocate( RhoOld(M,M), RhoMid(M,M), RhoNew(M,M) )
   allocate( Bmat(M,M), Dmat(M,M), Tmat(M,M) )

   if ( first_step ) then
   if ( rstinp ) then
      open( unit=rstinp_funit, file=rstinp_fname )
      print*,'Using restart'
      call rstload( rstinp_funit, Natom, qm_forces_total, M, RhoSaveA, RhoSaveB )
      close( unit=rstinp_funit )
   endif
   endif


! Update velocities
!------------------------------------------------------------------------------!
   do ii=1,natom
   do kk=1,3
      nucvel(kk,ii)=nucvel(kk,ii)+(1.5d0)*dt_elec*qm_forces_total(kk,ii)/atom_mass(ii)
   enddo
   enddo

! Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
   Energy=0.0d0
   call RMMcalc0_Init()
   call RMMcalc1_Overlap(Smat,Energy)
   call ehren_cholesky(M,Smat,Lmat,Umat,Linv,Uinv,Sinv)
   print*,'holio'

! Esto deja la Rho correcta en RMM, pero habria que ordenarlo mejor
   RhoMid=RhoSaveB
   if (.not.first_step) then
      RhoMid=matmul(RhoMid,Linv)
      RhoMid=matmul(Uinv,RhoMid)
   endif
   print*,'holix'
   call RMMcalc2_FockMao(RhoMid,Fock,DipMom,Energy)
   print*,'holiy'
   call calc_forceDS(natom,M,nucpos,nucvel,RhoMid,Fock,Sinv,Bmat,qm_forces_ds)
   print*,'holiz'


! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
   Fock=matmul(Fock,Uinv)
   Fock=matmul(Linv,Fock)
   Dmat=calc_Dmat(M,Linv,Uinv,Bmat)
   Tmat=DCMPLX(Fock)+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
   print*,'holib'

   RhoOld=RhoSaveA
   RhoMid=RhoSaveB
   if (first_step) then
      RhoMid=matmul(RhoMid,Lmat)
      RhoMid=matmul(Umat,RhoMid)
      call ehren_verlet_e(M,-(dt_elec/2.0d0),Tmat,RhoMid,RhoMid,RhoOld)
   endif
   call ehren_verlet_e(M,dt_elec,Tmat,RhoOld,RhoMid,RhoNew)
   RhoSaveA=RhoMid
   RhoSaveB=RhoNew


! Saving restart
!------------------------------------------------------------------------------!
!   call ehrensave( step_number, 
   print*,'holia'
   save_this_step = .false.
   if ( rstout ) then
      if ( rstout_nfreq > 0 ) then
         if ( modulo(step_number,rstout_nfreq) == 1 ) save_this_step = .true.
      endif
      if ( step_number == (ndyn_steps+1) ) save_this_step = .true.
   endif

   if ( save_this_step ) then
      open( unit=rstout_funit, file=rstout_fname )
      call rstsave( rstout_funit, Natom, qm_forces_total, M, RhoSaveA, RhoSaveB )
      close( unit=rstout_funit )
   endif

! Calculation of the dipole moment
!------------------------------------------------------------------------------!
   if (first_step) then
      call write_dipole(dipxyz, 0, 134, .true.)
      total_time=0.0d0
   else
      call dip(dipxyz)
      dipole_norm = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)
      call write_dipole(dipxyz, dipole_norm, 134, .false.)  

      print*,''
      print*,' Timer: ',total_time
      print*,''
      total_time=total_time+dt_elec*0.0241888d0
   endif

   DipMom_o = DipMom
   Energy_o = stored_energy
   stored_energy = Energy

   deallocate( Smat, Sinv )
   deallocate( Lmat, Umat, Linv, Uinv )
   deallocate( Fock, FockInt )
   deallocate( RhoOld, RhoMid, RhoNew )
   deallocate( Bmat, Dmat, Tmat )
   call g2g_timer_stop('ehrendyn call')

901 format(F15.9,2x,F15.9)
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
