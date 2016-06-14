!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenstep_verlet &
!------------------------------------------------------------------------------!
  ( dens_old, dens_now, dens_new, force_term, system_dipole, system_energy, &
  & nucmass, nucpos, nucvel, nucfor, time_step, &
  & is_first_step, info )
!
! dens_xxx should be in the ON base
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm,    only: rmmcalc_focknucl, rmmcalc_fockdens

  implicit none
  complex*16, intent(inout), dimension(:,:) :: dens_old
  complex*16, intent(inout), dimension(:,:) :: dens_now
  complex*16, intent(inout), dimension(:,:) :: dens_new
  real*8,     intent(inout), dimension(:,:) :: force_term
  real*8,     intent(inout), dimension(3)   :: system_dipole
  real*8,     intent(inout)                 :: system_energy
  real*8,     intent(in),    dimension(:)   :: nucmass
  real*8,     intent(in),    dimension(:,:) :: nucpos
  real*8,     intent(in),    dimension(:,:) :: nucvel
  real*8,     intent(in),    dimension(:,:) :: nucfor
  real*8,     intent(in)                    :: time_step

  logical, intent(in),  optional :: is_first_step
  integer, intent(out), optional :: info

  integer :: local_stat
  integer :: numof_atoms
  integer :: numof_basis
  real*8  :: dt
  real*8  :: newvel
  integer :: nn, kk, ii, jj, idx

  real*8,     allocatable, dimension(:,:) :: Smat, Sinv
  real*8,     allocatable, dimension(:,:) :: Lmat, Umat, Linv, Uinv
  real*8,     allocatable, dimension(:,:) :: Fock
  complex*16, allocatable, dimension(:,:) :: Dens
  real*8,     allocatable, dimension(:,:) :: Bmat, Dmat
  complex*16, allocatable, dimension(:,:) :: Tmat
  real*8,     allocatable, dimension(:,:) :: nucvel_now

  call g2g_timer_start('ehrenstep_verlet')
!
!
!
! Preliminary sets and checks
!------------------------------------------------------------------------------!
  numof_atoms = size( force_term, 2 )
  numof_basis = size( dens_old, 1)

  local_stat = 0

  if ( numof_basis /= size(dens_old,2) )  local_stat=1
  if ( numof_basis /= size(dens_now,1) )  local_stat=1
  if ( numof_basis /= size(dens_now,2) )  local_stat=1
  if ( numof_basis /= size(dens_new,1) )  local_stat=1
  if ( numof_basis /= size(dens_new,2) )  local_stat=1

  if ( numof_atoms /= size(nucmass,1) )   local_stat=2
  if ( numof_atoms /= size(nucpos,2) )    local_stat=2
  if ( numof_atoms /= size(nucvel,2) )    local_stat=2
  if ( numof_atoms /= size(nucfor,2) )    local_stat=2

  if ( 3 /= size(force_term,1) )  local_stat=3
  if ( 3 /= size(nucpos,1) )      local_stat=3
  if ( 3 /= size(nucvel,1) )      local_stat=3
  if ( 3 /= size(nucfor,1) )      local_stat=3


  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 1
      return
    else
      print*,'ehrenstep_verlet : incompatible size between arguments'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
!
! Allocation of memory
!------------------------------------------------------------------------------!
  nn=numof_basis
  allocate( Smat(nn,nn), Sinv(nn,nn), Fock(nn,nn), Dens(nn,nn), &
          & Lmat(nn,nn), Umat(nn,nn), Linv(nn,nn), Uinv(nn,nn), &
          & Bmat(nn,nn), Dmat(nn,nn), Tmat(nn,nn), &
          & nucvel_now(3,numof_atoms), stat=local_stat )

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 2
      return
    else
      print*,'ehrenstep_verlet : incompatible size between arguments'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
!
! Nuclear Force Calculations (matrices in AO)
!------------------------------------------------------------------------------!
  dt = time_step
  do ii=1,numof_atoms
  do kk=1,3
    newvel = (1.5d0) * time_step * nucfor(kk,ii) / nucmass(ii)
    nucvel_now(kk,ii) = nucvel(kk,ii) + newvel
    ! 1.5 because nucvel corresponds to one step before
  end do
  end do

  system_energy=0.0d0
  call RMMcalc0_Init()
  call RMMcalc1_Overlap( Smat, system_energy )
!  call rmmcalc_focknucl( Smat, system_energy )
  call ehren_cholesky( Smat, Sinv, Lmat, Umat, Linv, Uinv )

  Dens = dens_now
  Dens = matmul( Dens, Linv ) ! to atomic base
  Dens = matmul( Uinv, Dens ) ! to atomic base
  call RMMcalc2_FockMao( Dens, Fock, system_energy)
!  call rmmcalc_fockdens( Dens, efield, Fock, dipole, energy )
  call calc_forceDS( numof_atoms, numof_basis, nucpos, nucvel_now, &
                   & Dens, Fock, Sinv, Bmat, force_term )
  Fock = matmul( Fock, Uinv )
  Fock = matmul( Linv, Fock )
!
!
!
! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
  Dmat = calc_Dmat( numof_basis, Linv, Uinv, Bmat )
  Tmat = DCMPLX(Fock) + DCMPLX(0.0d0,1.0d0) * DCMPLX(Dmat)

  if ( present(is_first_step) ) then
    if ( is_first_step ) then
      dt = -( time_step / 2.0d0 )
      call ehren_verlet_e( numof_basis, dt, Tmat, dens_now, dens_now, dens_old )
    endif
  endif
  dt = time_step

  write(997,*) 'time step: ', dt, numof_basis
  do ii=1,size(dens_old,1)
  do jj=1,size(dens_old,2)
!    write(997,*) Tmat(ii,jj)
  enddo
  enddo
  call ehren_verlet_e( numof_basis, dt, Tmat, dens_old, dens_now, dens_new )

  do ii=1,size(dens_old,1)
  do jj=1,size(dens_old,2)
    write(997,*) dens_new(ii,jj)
  enddo
  enddo


  call g2g_timer_stop('ehrenstep_verlet')
!
!
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
