!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine Calculate_Overlap(Smat)
  use garcha_mod, only:M,RMM

  implicit none
  real*8,intent(out)    :: Smat(M,M)
  real*8                :: Energy
  integer               :: ii,jj,idx,idx0


! Calculate new Overlap in RMM
!------------------------------------------------------------------------------!
  call int1(Energy)


! Extract FockMao from RMM
!------------------------------------------------------------------------------!
  idx0=M*(M+1)
  do ii=1,M
  do jj=1,ii
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     Smat(ii,jj)=RMM(idx)
     Smat(jj,ii)=RMM(idx)
  enddo
  enddo


  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
