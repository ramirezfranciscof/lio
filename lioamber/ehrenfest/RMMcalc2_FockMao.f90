!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine RMMcalc2_FockMao( DensMao, FockMao, DipMom, Energy )
!
! Time is in ps
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use maskrmm
  use garcha_mod, only:M,Md,RMM,kkind,kkinds,cool,cools,igrid2, &
                     & total_time,field,Fx,Fy,Fz,epsilon,a0,    &
                     & natom,tdstep,Iz,NCO,Nunp
  implicit none
  complex*16,intent(in) :: DensMao(M,M)
  real*8,intent(out)    :: FockMao(M,M)
  real*8,intent(out)    :: DipMom(3)
  real*8,intent(out)    :: Energy

  real*8   :: Energy_1e
  real*8   :: Energy_Coulomb
  real*8   :: Energy_SolvT,Energy_SolvF
  real*8   :: Energy_Exchange
  real*8   :: Energy_Efield

  real*8   :: FieldNow(3)
  real*8   :: factor, g, Qc
  real*8   :: alpha, timeof_pert, timeof_peak
  real*8   :: time_gaus, dip_times_field, strange_term
  integer  :: kk,idx0
  integer  :: MM,MMd,igpu
  logical  :: MEMO

  logical  :: laser_is_on
  real*8   :: efield_shape, laser_freq, laser_long
  
!
!
! Calculate fixed-parts of fock
!--------------------------------------------------------------------!
  if (allocated(kkind))  deallocate(kkind)
  if (allocated(kkinds)) deallocate(kkinds)
  if (allocated(cool))   deallocate(cool)
  if (allocated(cools))  deallocate(cools)

  call g2g_reload_atom_positions(igrid2)
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()

  if (igpu.le.1) then
    call intsol(Energy_SolvF,Energy_SolvT,.true.)
  else
    call aint_qmmm_fock(Energy_SolvF,Energy_SolvT)
  endif

  call int2()

  if (igpu.gt.2) call aint_coulomb_init()
  if (igpu.eq.5) MEMO = .false.
  if (MEMO) then
    call g2g_timer_start('int3mem')
    call int3mem()
    call g2g_timer_stop('int3mem')
  endif

!
! Calculate unfixed Fock in RMM
!--------------------------------------------------------------------!
  call rmmput_dens(DensMao)
  call g2g_timer_start('g2g-solve + int3lu')
  call int3lu(Energy_Coulomb)
  call g2g_solve_groups(0,Energy_Exchange,0)
  call g2g_timer_stop('g2g-solve + int3lu')

  Energy_Efield = 0.0d0
  alpha       = 0.2 / ( tdstep * 0.0241888 )**2
  timeof_peak =  50 * ( tdstep * 0.0241888 )
  timeof_pert = 100 * ( tdstep * 0.0241888 )
  if ( field .and. ( total_time <= timeof_pert ) ) then
    Qc=-2*NCO+Nunp
    do kk=1,natom
      Qc=Qc+Iz(kk)
    end do

    g = 1.0d0
    factor = 2.54d0
    time_gaus = exp( (-alpha) * ( total_time - timeof_peak )**2 )
    FieldNow(1) = Fx * time_gaus
    FieldNow(2) = Fy * time_gaus
    FieldNow(3) = Fz * time_gaus
    write(999,*) ' ==> ', FieldNow

    call dip( DipMom(1), DipMom(2), DipMom(3) )
    call intfld( g, FieldNow(1), FieldNow(2), FieldNow(3) )
    dip_times_field = 0.0d0
    dip_times_field = dip_times_field + FieldNow(1) * DipMom(1)
    dip_times_field = dip_times_field + FieldNow(2) * DipMom(2)
    dip_times_field = dip_times_field + FieldNow(3) * DipMom(3)
    strange_term = (0.5d0) * (1.0d0 - 1.0d0/epsilon) * Qc**2 / a0

    Energy_Efield = Energy_Efield - g * dip_times_field / factor
    Energy_Efield = Energy_Efield - strange_term
  end if

  laser_is_on = .true.
  laser_long = 139.3 !nm
  laser_freq = (6.28318530718) * (299.792458) / (laser_long)
  ! [freq fs-1] = [2pi] * [c in nm/fs] / [long]
  if ( laser_is_on ) then
    Qc=-2*NCO+Nunp
    do kk=1,natom
      Qc=Qc+Iz(kk)
    end do

    g = 1.0d0
    factor = 2.54d0
    efield_shape = sin( laser_freq * total_time )
    FieldNow(1) = Fx * efield_shape
    FieldNow(2) = Fy * efield_shape
    FieldNow(3) = Fz * efield_shape
    print*,'Doing laser: ',FieldNow

    call dip( DipMom(1), DipMom(2), DipMom(3) )
    call intfld( g, FieldNow(1), FieldNow(2), FieldNow(3) )
    dip_times_field = 0.0d0
    dip_times_field = dip_times_field + FieldNow(1) * DipMom(1)
    dip_times_field = dip_times_field + FieldNow(2) * DipMom(2)
    dip_times_field = dip_times_field + FieldNow(3) * DipMom(3)
    strange_term = (0.5d0) * (1.0d0 - 1.0d0/epsilon) * Qc**2 / a0

    Energy_Efield = Energy_Efield - g * dip_times_field / factor
    Energy_Efield = Energy_Efield - strange_term
  endif
!
! Calculate Energy
!--------------------------------------------------------------------!
  MM=M*(M+1)/2
  MMd=Md*(Md+1)/2
  idx0=3*MM+2*MMd
  Energy_1e=0.0d0
  do kk=1,MM
    Energy_1e=Energy_1e+RMM(kk)*RMM(idx0+kk)
  enddo

!  Energy=0.0d0
  Energy=Energy+Energy_1e
  Energy=Energy+Energy_Coulomb
  Energy=Energy+Energy_SolvT
  Energy=Energy+Energy_Exchange
  Energy=Energy+Energy_Efield

!
! Extract FockMao from RMM
!--------------------------------------------------------------------!
  call rmmget_fock(FockMao)


  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
