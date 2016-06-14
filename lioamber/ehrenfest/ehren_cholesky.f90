!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehren_cholesky( Smat, Sinv, Lmat, Umat, Linv, Uinv, info)
!------------------------------------------------------------------------------!
!
! TO DO: abstract the decomposition and the inversion
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  implicit none
  real*8, intent(in), dimension(:,:)   :: Smat
  real*8, intent(out), dimension(:,:)  :: Sinv
  real*8, intent(out), dimension(:,:)  :: Lmat, Umat
  real*8, intent(out), dimension(:,:)  :: Linv, Uinv
  integer, intent(out), optional       :: info

  integer :: local_stat
  integer :: Msize
  integer :: ii,jj,iost
!
!
! Checks
!------------------------------------------------------------------------------!
  local_stat = 0
  Msize = size( Smat, 1 )  
  if ( Msize /= size(Smat,2) )  local_stat=1
  if ( Msize /= size(Sinv,1) )  local_stat=1
  if ( Msize /= size(Sinv,2) )  local_stat=1
  if ( Msize /= size(Lmat,1) )  local_stat=1
  if ( Msize /= size(Lmat,2) )  local_stat=1
  if ( Msize /= size(Umat,1) )  local_stat=1
  if ( Msize /= size(Umat,2) )  local_stat=1
  if ( Msize /= size(Linv,1) )  local_stat=1
  if ( Msize /= size(Linv,2) )  local_stat=1
  if ( Msize /= size(Uinv,1) )  local_stat=1
  if ( Msize /= size(Uinv,2) )  local_stat=1

  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 1
      return
    else
      print*,'ehren_cholesky : incompatible size between arguments'
      print*,'local info: ', local_stat
      stop
    end if
  end if
!
!
! Decomposition
!------------------------------------------------------------------------------!
  Lmat=Smat
  call dpotrf('L', Msize, Lmat, Msize, local_stat)
  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 2
      return
    else
      print*,'ehren_cholesky : problem using dpotrf'
      print*,'local info: ', local_stat
      stop
    end if
  end if

  do ii=1,Msize-1
  do jj=ii+1,Msize
    Lmat(ii,jj)=0.0d0
  enddo
  enddo
!
!
! Inversion
!------------------------------------------------------------------------------!
  Linv=Lmat
  call dtrtri('L', 'N', Msize, Linv, Msize, local_stat)
  if ( local_stat /= 0 ) then
    if ( present(info) ) then
      info = 3
      return
    else
      print*,'ehren_cholesky : problem using dtrtri'
      print*,'local info: ', local_stat
      stop
    end if
  end if

  Umat=transpose(Lmat)
  Uinv=transpose(Linv)
  Sinv=matmul(Uinv,Linv)
!
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
