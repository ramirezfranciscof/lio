!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  function forceDS(natoms,nbasis,nucpos,nucvel,Pao,Fao,Sinv)
!--------------------------------------------------------------------!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use testmod
  implicit none
  real*8                :: forceDS(3,Natoms)
  integer,intent(in)    :: natoms
  integer,intent(in)    :: nbasis
  real*8,intent(in)     :: nucpos(3,natoms)
  real*8,intent(in)     :: nucvel(3,natoms)
  complex*16,intent(in) :: Pao(Nbasis,nbasis)
  real*8,intent(in)     :: Fao(Nbasis,nbasis)
  real*8,intent(in)     :: Sinv(Nbasis,nbasis)

  complex*16,allocatable :: Matin(:,:),Mat1(:,:),Mat2(:,:)
  real*8,allocatable     :: MatB(:,:),TrpB(:,:)
  complex*16,allocatable :: fterm1(:,:),fterm2(:,:),fterm3(:,:)
!--------------------------------------------------------------------!
  allocate(Matin(nbasis,nbasis))
  allocate(Mat1(nbasis,nbasis),Mat2(nbasis,nbasis))
  allocate(MatB(nbasis,nbasis),TrpB(nbasis,nbasis))
  allocate(fterm1(3,natoms),fterm2(3,natoms),fterm3(3,natoms))

  fterm1=dcmplx(0.0d0,0.0d0)
  fterm2=dcmplx(0.0d0,0.0d0)
  fterm3=dcmplx(0.0d0,0.0d0)

  Mat1=matmul(Pao,Fao)
  Mat1=matmul(Mat1,Sinv)
  Mat2=matmul(Sinv,Fao)
  Mat2=matmul(Mat1,Pao)
  Matin=transpose(Mat1)+Mat2
  call forceDS_dss(natoms,nbasis,nucpos,nucvel,Matin,MatB,fterm1)
  TrpB=transpose(MatB)

!  Mat1=matmul(Pao,TrpB)
!  Mat1=matmul(Mat1,Sinv)
!  Mat1=Mat1*dcmplx(0.0d0, 1.0d0)
!  Mat2=matmul(Sinv,MatB)
!  Mat2=matmul(Mat2,Pao)
!  Mat2=Mat2*dcmplx(0.0d0,-1.0d0)
!  Matin=transpose(Mat1)+Mat2
!  call forceDS_dss(natoms,nbasis,nucpos,nucvel,Matin,MatB,fterm2)

!  Mat1=Pao*dcmplx(0.0d0,-1.0d0)
!  Mat2=Pao*dcmplx(0.0d0, 1.0d0)
!  Matin=transpose(Mat1)+Mat2
!  call forceDS_dds(natoms,nbasis,nucpos,nucvel,Matin,fterm3)

  forceDS=dble(real(fterm1+fterm2+fterm3))

  return;end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
