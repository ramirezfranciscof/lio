!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine ehrendyn(Energy,Dipmom)
!------------------------------------------------------------------------------!
!
! RhoSave should always remain in OA basis
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use ehrensubs
  use garcha_mod, only: M, natom, nucpos, nucvel, RhoSaveA, RhoSaveB & 
                      , qm_forces_ds, tdstep, first_step, rmm, atom_mass &
                      , total_time, first_step, qm_forces_total !, Fock_a, Fock_b
  implicit none
  real*8,intent(inout) :: Energy,Dipmom(3)

  real*8,allocatable,dimension(:,:)     :: Smat,Sinv,Lmat,Umat,Linv,Uinv
  real*8,allocatable,dimension(:,:)     :: Fock,FockInt
  real*8,allocatable,dimension(:,:)     :: DensAO
  complex*16,allocatable,dimension(:,:) :: RhoOld,RhoMid,RhoNew

  real*8,allocatable,dimension(:,:)     :: Bmat,Dmat
  complex*16,allocatable,dimension(:,:) :: Tmat
  real*8                                :: dt
  integer                               :: ii,jj,kk,idx
  real*8                                :: ux,uy,uz

! Preliminaries
!------------------------------------------------------------------------------!
  allocate(Smat(M,M),Sinv(M,M),Lmat(M,M),Umat(M,M),Linv(M,M),Uinv(M,M))
  allocate(Fock(M,M),FockInt(M,M))
  allocate(DensAO(M,M))
  allocate(RhoOld(M,M),RhoMid(M,M),RhoNew(M,M))
  allocate(Bmat(M,M),Dmat(M,M),Tmat(M,M))

  if (.not.allocated(qm_forces_total)) then
    allocate(qm_forces_total(3,natom))
    qm_forces_total=0.0d0
  endif

  if (.not.allocated(qm_forces_ds)) then
    allocate(qm_forces_ds(3,natom))
    qm_forces_ds=0.0d0
  endif

  dt=tdstep


! Update velocities
!------------------------------------------------------------------------------!
  do ii=1,natom
  do kk=1,3
    nucvel(kk,ii)=nucvel(kk,ii)+(1.5)*dt*qm_forces_total(kk,ii)/atom_mass(ii)
  enddo
  enddo


! Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
  call Calculate_Overlap(Smat)
  call ehren_cholesky(M,Smat,Lmat,Umat,Linv,Uinv,Sinv)

  DensAO=real(RhoSaveB)
!  DensAO=RealPart(RhoSaveB)
!  DensAO=matmul(DensAO,Linv)
!  DensAO=matmul(Uinv,DensAO)

! Esto deja la Rho correcta en RMM, pero habria que ordenarlo mejor
  call Calculate_Fock(DensAO,Fock,Energy)
  call calc_forceDS(natom,M,nucpos,nucvel,RhoSaveB,Fock,Sinv,Bmat,qm_forces_ds)


! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
  Fock=matmul(Fock,Uinv)
  Fock=matmul(Linv,Fock)
  Dmat=calc_Dmat(M,Linv,Uinv,Bmat)
  Tmat=DCMPLX(Fock)+CMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)

  RhoOld=RhoSaveA
  RhoMid=matmul(RhoMid,Lmat)
  RhoMid=matmul(Umat,RhoMid)
  RhoMid=RhoSaveB
  RhoMid=matmul(RhoMid,Lmat)
  RhoMid=matmul(Umat,RhoMid)

  if (first_step) then
    call ehren_verlet_e(M,-(dt/2.0d0),Tmat,RhoMid,RhoMid,RhoOld)
  endif
  call ehren_verlet_e(M,dt,Tmat,RhoOld,RhoMid,RhoNew)

  RhoMid=matmul(RhoMid,Linv)
  RhoMid=matmul(Uinv,RhoMid)
  RhoSaveA=RhoMid
  RhoNew=matmul(RhoNew,Linv)
  RhoNew=matmul(Uinv,RhoNew)
  RhoSaveB=RhoNew

  do jj=1,M
    ii=jj
    idx=ii+(2*M-jj)*(jj-1)/2
    RMM(idx)=(REAL(RhoMid(ii,jj)))
    do ii=jj+1,M
       idx=ii+(2*M-jj)*(jj-1)/2
       RMM(idx)=REAL(RhoMid(ii,jj))+REAL(RhoMid(jj,ii))
    enddo
  enddo


!!  do_magnus=.false.
!!  if (do_magnus) then
!!    Tmat=DCMPLX(Fock)
!    Tmat=DCMPLX((7.0d0/4.0d0)*FockB-(3.0d0/4.0d0)*FockA)
!    Tmat=Tmat+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
!!    call ehren_magnus(M,10,dt/2,Tmat,DensOld,DensInt)

!!    DensInt=matmul(DensInt,Linv)
!!    DensInt=matmul(Uinv,DensInt)
!!    call Calculate_Fock(DensInt,FockInt,Energy)
!!    FockInt=matmul(FockInt,Uinv)
!!    FockInt=matmul(Linv,FockInt)

!!    Tmat=DCMPLX(FockInt)
!    Tmat=Tmat+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
!!    call ehren_magnus(M,10,dt,Tmat,DensOld,DensNew)
!!  else
!!    Tmat=DCMPLX(Fock)+CMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
!!    if (first_step) call ehren_verlet_e(M,-(dt/2.0d0),Tmat,DensOld,DensOld,RhoCero)
!!    call ehren_verlet_e(M,dt,Tmat,RhoCero,DensOld,DensNew)
!!    RhoCero=DensOld
!!  endif



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
