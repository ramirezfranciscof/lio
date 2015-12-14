!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine Calculate_Fock(DensMao,FockMao)
  use garcha_mod, only:M,RMM,kkind,kkinds,cool,cools

  implicit none
  complex*16,intent(in) :: DensMao(M,M)
  real*8,intent(out)    :: FockMao(M,M)
  real*8                :: Energy1,Energy2
  integer               :: ii,jj,idx,idx0
  integer               :: igpu


! Store DensMao in RMM
!------------------------------------------------------------------------------!
  do jj=1,M
    do ii=jj,M
       idx=ii+(2*M-jj)*(jj-1)/2
       RMM(idx)=(REAL(DensMao(ii,jj)))*2
    enddo
    idx=jj+(2*M-jj)*(jj-1)/2
    RMM(idx)=RMM(idx)/2
  enddo


! Calculate new Fock in RMM
!------------------------------------------------------------------------------!
  if (allocated(kkind))  deallocate(kkind)
  if (allocated(kkinds)) deallocate(kkinds)
  if (allocated(cool))   deallocate(cool)
  if (allocated(cools))  deallocate(cools)

! Cuando se mueven los nucleos
  call int1(Energy1)
  call aint_query_gpu_level(igpu)
  if (igpu.le.1) then
    call intsol(Energy1,Energy2,.true.)
  else
    call aint_qmmm_fock(Energy1,Energy2)
  endif
  call int2()
  call int3mem()

! Si los nucleaos estan quietos
  call int3lu(Energy1)
  call g2g_solve_groups(0,Energy1,0)


! Extract FockMao from RMM
!------------------------------------------------------------------------------!
  idx=0
  idx0=M*(M+1)
  do jj=1,M
  do ii=jj,M
     idx=ii+(2*M-jj)*(jj-1)/2+idx0
     FockMao(ii,jj)=RMM(idx)
     FockMao(jj,ii)=RMM(idx)
     write(603,*) idx,jj,ii,RMM(idx)
  enddo
  enddo

  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
