!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine ehrendyn(Energy,Dipmom)
!------------------------------------------------------------------------------!
!
! Description Pending
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use ehrensubs
  use garcha_mod, only: M, natom, nucpos, nucvel, RhoSave, qm_forces_ds &
                      , tdstep, first_step, rmm, atom_mass, total_time  &
                      , first_step, RhoCero !, Fock_a, Fock_b
  implicit none
  real*8,intent(inout) :: Energy,Dipmom(3)

  complex*16,allocatable,dimension(:,:) :: DensOld,DensInt,DensNew
  real*8,allocatable,dimension(:,:)     :: Fock,FockInt
  real*8,allocatable,dimension(:,:)     :: Bmat,Dmat,Tmat
  real*8,allocatable,dimension(:,:)     :: Smat,Sinv,Lmat,Umat,Linv,Uinv
  real*8                                :: dt
  integer                               :: ii,jj,kk,idx

  real*8                                :: RhoReal,RhoImag,ux,uy,uz
  logical                               :: do_magnus

! Preliminaries
!------------------------------------------------------------------------------!
  allocate(Fock(M,M),FockInt(M,M))
  allocate(DensOld(M,M),DensInt(M,M),DensNew(M,M))
  allocate(Smat(M,M),Sinv(M,M))
  allocate(Lmat(M,M),Umat(M,M),Linv(M,M),Uinv(M,M))
  allocate(Bmat(M,M),Dmat(M,M),Tmat(M,M))

  if (.not.allocated(qm_forces_ds)) then
    allocate(qm_forces_ds(3,natom))
    qm_forces_ds=0.0d0
  endif
  dt=tdstep


! Update velocities
!------------------------------------------------------------------------------!
  do ii=1,natom
  do kk=1,3
    nucvel(kk,ii)=nucvel(kk,ii)+dt*qm_forces_ds(kk,ii)/atom_mass(ii)
  enddo
  enddo


! Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
  call Calculate_Overlap(Smat)
  call ehren_cholesky(M,Smat,Lmat,Umat,Linv,Uinv,Sinv)
!  if (.not.first_step) then
!    DensOld=matmul(DensOld,Linv)
!    DensOld=matmul(Uinv,DensOld)
!  endif


  call Calculate_Fock(RhoSave,Fock,Energy)
  call calc_forceDS(natom,M,nucpos,nucvel,RhoSave,Fock,Sinv,Bmat,qm_forces_ds)
  write(888,*) '## Energy = ',Energy


! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
  DensOld=RhoSave
  DensOld=matmul(DensOld,Lmat)
  DensOld=matmul(Umat,DensOld)
  Fock=matmul(Fock,Uinv)
  Fock=matmul(Linv,Fock)
  Dmat=calc_Dmat(M,Linv,Uinv,Bmat)


  do_magnus=.false.
  if (do_magnus) then
    Tmat=DCMPLX(Fock)
!    Tmat=DCMPLX((7.0d0/4.0d0)*FockB-(3.0d0/4.0d0)*FockA)
!    Tmat=Tmat+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
    call ehren_magnus(M,10,dt/2,Tmat,DensOld,DensInt)

    DensInt=matmul(DensInt,Linv)
    DensInt=matmul(Uinv,DensInt)
    call Calculate_Fock(DensInt,FockInt,Energy)
    FockInt=matmul(FockInt,Uinv)
    FockInt=matmul(Linv,FockInt)

    Tmat=DCMPLX(FockInt)
!    Tmat=Tmat+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
    call ehren_magnus(M,10,dt,Tmat,DensOld,DensNew)
  else
    if (first_step) call ehren_verlet_e(M,-(dt/2.0d0),Fock,DensOld,DensOld,RhoCero)
    call ehren_verlet_e(M,dt,Fock,RhoCero,DensOld,DensNew)
    RhoCero=DensOld
  endif


  RhoSave=DensNew
  RhoSave=matmul(DensNew,Linv)
  RhoSave=matmul(Uinv,RhoSave)


! Store Dens in RMM
!------------------------------------------------------------------------------!
  DensOld=matmul(DensOld,Linv)
  DensOld=matmul(Uinv,DensOld)
  do jj=1,M
    ii=jj
    idx=ii+(2*M-jj)*(jj-1)/2
    RMM(idx)=(REAL(DensOld(ii,jj)))
    do ii=jj+1,M
       idx=ii+(2*M-jj)*(jj-1)/2
       RMM(idx)=REAL(DensOld(ii,jj))+REAL(DensOld(jj,ii))
    enddo
  enddo


! Calculation of the dipole moment
!------------------------------------------------------------------------------!
  if (first_step) then
    open(unit=134,file='x.dip')
    open(unit=135,file='y.dip')
    open(unit=136,file='z.dip')
    write(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
    write(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
    write(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
    total_time=0.0d0
  endif

  call dip(ux,uy,uz)
  write(134,901) total_time,ux
  write(135,901) total_time,uy
  write(136,901) total_time,uz
  print*,''
  print*,' Timer: ',total_time
  print*,''
  total_time=total_time+dt*0.0241888

 901 format(F15.9,2x,F15.9)

  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
