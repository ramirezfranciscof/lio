!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehren_dipole &
!------------------------------------------------------------------------------!
  ( dens_old, dens_now, dens_new, force_term, system_dipole, system_energy &
  & nucpos, nucvel, nucfor, time_step, is_first_step )
!------------------------------------------------------------------------------!
!
! RhoSaveA and RhoSaveB are stored in ON basis, except for the first step
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm,    only: rmmcalc_focknucl, rmmcalc_fockdens

  implicit none
  real*8,  intent(inout) :: dens_old(:,:)
  real*8,  intent(inout) :: dens_now(:,:)
  real*8,  intent(inout) :: dens_new(:,:)
  real*8,  intent(inout) :: force_term(3,:)
  real*8,  intent(inout) :: system_dipole(3)
  real*8,  intent(inout) :: system_energy
  real*8,  intent(in)    :: nucpos(3,:)
  real*8,  intent(in)    :: nucvel(3,:)
  real*8,  intent(in)    :: nucfor(3,:)
  real*8,  intent(in)    :: time_step
  logical, intent(in)    :: is_first_step

  real*8,allocatable,dimension(:,:)     :: Smat,Sinv,Lmat,Umat,Linv,Uinv
  real*8,allocatable,dimension(:,:)     :: Fock,FockInt
  complex*16,allocatable,dimension(:,:) :: RhoOld,RhoMid,RhoNew

  real*8,allocatable,dimension(:,:)     :: Bmat,Dmat
  complex*16,allocatable,dimension(:,:) :: Tmat
  real*8                                :: dt
  real*8                                :: ux,uy,uz
  integer                               :: nn,kk,ii,jj,idx

  call g2g_timer_start('ehrendyn step')
!
!
! Preliminaries
!------------------------------------------------------------------------------!
  dt=time_step
  nn=numof_basis
  allocate( Smat(nn,nn), Sinv(nn,nn) )
  allocate( Lmat(nn,nn), Umat(nn,nn), Linv(nn,nn), Uinv(nn,nn) )
  allocate( Fock(nn,nn), FockInt(nn,nn) )
  allocate( RhoOld(nn,nn), RhoMid(nn,nn), RhoNew(nn,nn) )
  allocate( Bmat(nn,nn), Dmat(nn,nn), Tmat(nn,nn) )

  if ( .not.allocated(qm_forces_total) ) then
    allocate( qm_forces_total(3,numof_atoms) )
    qm_forces_total=0.0d0
  endif

  if ( .not.allocated(qm_forces_ds) ) then
    allocate( qm_forces_ds(3,numof_atoms) )
    qm_forces_ds=0.0d0
  endif

! Update velocities
!------------------------------------------------------------------------------!
  do ii=1,numof_atoms
  do kk=1,3
    nucvel(kk,ii)=nucvel(kk,ii)+(1.5)*dt*qm_forces_total(kk,ii)/atom_mass(ii)
  enddo
  enddo

! Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
  Energy=0.0d0
  call rmmcalc_focknucl( Smat, Energy )

  call RMMcalc0_Init()
  call RMMcalc1_Overlap(Smat,Energy)
  call ehren_cholesky(numof_basis,Smat,Lmat,Umat,Linv,Uinv,Sinv)

! Esto deja la Rho correcta en RMM, pero habria que ordenarlo mejor
  RhoMid=RhoSaveB
  if (.not.is_first_step) then
    RhoMid=matmul(RhoMid,Linv)
    RhoMid=matmul(Uinv,RhoMid)
  endif
  call rmmcalc_fockdens( RhoMid, efield, Fock, Dipole, Energy )
  call calc_forceDS(numof_atoms,numof_basis,nucpos,nucvel,RhoMid,Fock,Sinv,Bmat,qm_forces_ds)



! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
  Fock=matmul(Fock,Uinv)
  Fock=matmul(Linv,Fock)
  Dmat=calc_Dmat(numof_basis,Linv,Uinv,Bmat)
  Tmat=DCMPLX(Fock)+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)

  RhoOld=RhoSaveA
  RhoMid=RhoSaveB
  if (is_first_step) then
    RhoMid=matmul(RhoMid,Lmat)
    RhoMid=matmul(Umat,RhoMid)
    call ehren_verlet_e(numof_basis,-(dt/2.0d0),Tmat,RhoMid,RhoMid,RhoOld)
  endif
  call ehren_verlet_e(numof_basis,dt,Tmat,RhoOld,RhoMid,RhoNew)
  RhoSaveA=RhoMid
  RhoSaveB=RhoNew


! Calculation of the dipole moment
!------------------------------------------------------------------------------!
  if (is_first_step) then
    open(unit=134,file='x.dip')
    open(unit=135,file='y.dip')
    open(unit=136,file='z.dip')
    write(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
    write(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
    write(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
    total_time=0.0d0
  endif

  if (.not.is_first_step) then
   call dip(ux,uy,uz)
   write(134,901) total_time,ux
   write(135,901) total_time,uy
   write(136,901) total_time,uz
   print*,''
   print*,' Timer: ',total_time
   print*,''
   total_time=total_time+dt*0.0241888
  endif

  call g2g_timer_stop('ehrendyn step')

 901 format(F15.9,2x,F15.9)

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
