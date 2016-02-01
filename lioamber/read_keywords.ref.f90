!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine read_keywords(file_name)
!
! THIS SUBROUTINE SHOULD BETTER BE ADDED AS PAR OF GARCHA_MOD
! ALSO: GARCHA_MOD SHOULD BE SEPARATED IN DIFFERENT MODULES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, only:                                                        &
  natom, nsol, OPEN, NMAX, Nunp, VCINP, frestartin, GOLD, told, rmax, rmaxs,   &
  predcoef, idip, writexyz, intsoldouble, DIIS, ndiis, dgtrig, Iexch, integ,   &
  dens, igrid, igrid2, timedep, tdstep, ntdstep, propagator, NBCH, field,      &
  a0, epsilon, exter, Fx, Fy, Fz, tdrestart, writedens, basis_set,             &
  fitting_set, int_basis, cubegen_only, cube_res, cube_dens, cube_dens_file,   &
  cube_orb, cube_sel, cube_orb_file, cube_elec, cube_elec_file

  implicit none
  character(len=*),intent(in) :: file_name
  logical                     :: file_exists
  integer                     :: io_stat

  integer :: charge
  logical :: writeforces


  namelist /lio/                                                               &
  natom, nsol, OPEN, NMAX, Nunp, VCINP, frestartin, GOLD, told, rmax, rmaxs,   &
  predcoef, idip, writexyz, intsoldouble, DIIS, ndiis, dgtrig, Iexch, integ,   &
  dens, igrid, igrid2, timedep, tdstep, ntdstep, propagator, NBCH, field,      &
  a0, epsilon, exter, Fx, Fy, Fz, tdrestart, writedens, basis_set,             &
  fitting_set, int_basis, cubegen_only, cube_res, cube_dens, cube_dens_file,   &
  cube_orb, cube_sel, cube_orb_file, cube_elec, cube_elec_file,                &
  charge, writeforces

!------------------------------------------------------------------------------!

  inquire(file=file_name,exist=file_exists)
  if (.not.file_exists) then
     write(*,*) 'ERROR: Problem finding Input file "',adjustl(file_name),'".'
     stop
  endif

  open(unit=100,file=file_name,iostat=io_stat)
  if (io_stat.ne.0) then
     write(*,*) 'ERROR: Problem opening input file "',adjustl(file_name),'".'
     stop
  endif

  read(unit=100,nml=lio,iostat=io_stat)
  if (io_stat.ne.0) then
     write(*,*) 'ERROR: Problem reading input file "',adjustl(file_name),'".'
     stop
  endif

  close(unit=100,iostat=io_stat)
  if (io_stat.ne.0) then
     write(*,*) 'ERROR: Problem closing input file "',adjustl(file_name),'".'
     stop
  endif

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
