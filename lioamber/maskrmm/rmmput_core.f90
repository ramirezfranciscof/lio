!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine rmmput_core_r(CoreMao)
  use garcha_mod, only:M,Md,RMM
  implicit none
  real*4,intent(in) :: CoreMao(M,M)
  integer           :: MM,MMd,ii,jj,idx,idx0

  MM = M*(M+1)/2
  MMd = Md*(Md+1)/2
  idx0=M*(M+1)
  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     RMM(idx)=CoreMao(jj,ii)
  enddo
  enddo

  return; end subroutine
!--------------------------------------------------------------------!
  subroutine rmmput_core_d(CoreMao)
  use garcha_mod, only:M,Md,RMM
  implicit none
  real*8,intent(in) :: CoreMao(M,M)
  integer           :: MM,MMd,ii,jj,idx,idx0

  MM = M*(M+1)/2
  MMd = Md*(Md+1)/2
  idx0=M*(M+1)
  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     RMM(idx)=CoreMao(jj,ii)
  enddo
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
