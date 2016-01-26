!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine Calculate_Overlap(Smat)
  use garcha_mod, only:M,RMM,natom,d,r

  implicit none
  real*8,intent(out)    :: Smat(M,M)
  real*8                :: Energy
  integer               :: ii,jj,idx,idx0
  integer               :: igpu


! Calculate new Overlap in RMM
!------------------------------------------------------------------------------!
  do ii=1,natom
  do jj=1,natom
    d(ii,jj)=0.0d0
    d(ii,jj)=d(ii,jj)+(r(ii,1)-r(jj,1))**2
    d(ii,jj)=d(ii,jj)+(r(ii,2)-r(jj,2))**2
    d(ii,jj)=d(ii,jj)+(r(ii,3)-r(jj,3))**2
  enddo
  enddo

  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
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
