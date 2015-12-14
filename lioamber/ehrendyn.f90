!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine ehrendyn(Energy,Dipmom)
!------------------------------------------------------------------------------!
!
! Description Pending
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use ehrensubs
  use garcha_mod, only: M,natom,nucpos,nucvel,RhoSave,qm_forces_ds,tdstep
  implicit none
  real*8,intent(inout) :: Energy,Dipmom(3)

  complex*16,allocatable,dimension(:,:) :: Dens,DensInt
  real*8,allocatable,dimension(:,:)     :: Fock,FockInt
  real*8,allocatable,dimension(:,:)     :: Bmat,Dmat,Tmat
  real*8,allocatable,dimension(:,:)     :: Smat,Sinv,Lmat,Umat,Linv,Uinv
  real*8                                :: dt
  integer :: ii,kk


! Preliminaries
!------------------------------------------------------------------------------!
  allocate(Fock(M,M),FockInt(M,M),Dens(M,M),DensInt(M,M))
  allocate(Smat(M,M),Sinv(M,M))
  allocate(Lmat(M,M),Umat(M,M),Linv(M,M),Uinv(M,M))
  allocate(Bmat(M,M),Dmat(M,M),Tmat(M,M))

  if (allocated(qm_forces_ds)) deallocate(qm_forces_ds)
  allocate(qm_forces_ds(3,natom))

  dt=tdstep

! Nuclear Force Calculation (works in AO)
!------------------------------------------------------------------------------!
  Dens=RhoSave
  call Calculate_Overlap(Smat)
  call ehren_cholesky(M,Smat,Lmat,Umat,Linv,Uinv,Sinv)
  call Calculate_Fock(Dens,Fock)
  call calc_forceDS(natom,M,nucpos,nucvel,Dens,Fock,Sinv,Bmat,qm_forces_ds)
  print*,''
  print*,''
  print*,'     natom    dir'
  do ii=1,natom
  do kk=1,3
     print*,ii,kk,qm_forces_ds(kk,ii)
  enddo
  enddo
  print*,''
  print*,''


! Density Propagation (works in ON)
!------------------------------------------------------------------------------!
!  Dens=basechange(M,Umat,Dens,Lmat)
!  Fock=basechange(M,Linv,Fock,Uinv)
  Dens=matmul(Dens,Lmat)
  Dens=matmul(Umat,Dens)
  Fock=matmul(Fock,Uinv)
  Fock=matmul(Linv,Fock)
  Dmat=calc_Dmat(M,Lmat,Umat,Bmat)
!  Tmat=DCMPLX((7.0d0/4.0d0)*FockB-(3.0d0/4.0d0)*FockA)
  Tmat=DCMPLX(Fock)
  Tmat=Tmat+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
  DensInt=ehren_magnus(M,10,Tmat,Dens,dt/2.0)

  DensInt=matmul(DensInt,Linv)
  DensInt=matmul(Uinv,DensInt)
  call Calculate_Fock(DensInt,FockInt)
  FockInt=matmul(FockInt,Uinv)
  FockInt=matmul(Linv,FockInt)

  Tmat=DCMPLX(FockInt)
  Tmat=Tmat+DCMPLX(0.0d0,1.0d0)*DCMPLX(Dmat)
  Dens=ehren_magnus(M,10,Tmat,Dens,dt)
!  RhoSave=basechange(M,Uinv,Dens,Linv)
  RhoSave=matmul(Dens,Linv)
  RhoSave=matmul(Uinv,RhoSave)

  deallocate(Fock,FockInt,Dens,DensInt)
  deallocate(Smat,Sinv)
  deallocate(Lmat,Umat,Linv,Uinv)
  deallocate(Bmat,Dmat,Tmat)
  return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
