!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmget_fock_r( FockMao )
   use garcha_mod, only: M, RMM
   implicit none
   real*4, intent(out) :: FockMao(M,M)
   integer             :: ii, jj, idx, idx0

   idx0 = M * (M+1)
   do jj = 1, M
   do ii = jj, M
      idx = ii + (2*M-jj) * (jj-1) / 2
      FockMao(ii,jj) = RMM(idx+idx0)
      FockMao(jj,ii) = RMM(idx+idx0)
   enddo
   enddo

end subroutine rmmget_fock_r
!
!
!--------------------------------------------------------------------!
subroutine rmmget_fock_d( FockMao )
   use garcha_mod, only: M, RMM
   implicit none
   real*8, intent(out) :: FockMao(M,M)
   integer             :: ii, jj, idx, idx0

   idx0 = M * (M+1)
   do jj = 1, M
   do ii = jj, M
      idx = ii + (2*M-jj) * (jj-1) / 2
      FockMao(ii,jj) = RMM(idx+idx0)
      FockMao(jj,ii) = RMM(idx+idx0)
   enddo
   enddo

end subroutine rmmget_fock_d
!
!
!--------------------------------------------------------------------!
subroutine rmmget_fockos_r( FockMao_a, FockMao_b )
   use garcha_mod, only: M, RMM
   implicit none
   real*4, intent(out) :: FockMao_a(M,M)
   real*4, intent(out) :: FockMao_b(M,M)
   integer             :: ii, jj, idx, idx0a, idx0b

   idx0a = M * (M+1)
   idx0b = idx0a / 2
   do jj = 1, M
   do ii = jj, M
      idx = ii + (2*M-jj) * (jj-1) / 2 
      FockMao_a(ii,jj) = RMM(idx+idx0a)
      FockMao_a(jj,ii) = RMM(idx+idx0a)
      FockMao_b(ii,jj) = RMM(idx+idx0b)
      FockMao_b(jj,ii) = RMM(idx+idx0b)
   enddo
   enddo

end subroutine rmmget_fockos_r
!
!
!--------------------------------------------------------------------!
subroutine rmmget_fockos_d( FockMao_a, FockMao_b )
   use garcha_mod, only: M, RMM
   implicit none
   real*8, intent(out) :: FockMao_a(M,M)
   real*8, intent(out) :: FockMao_b(M,M)
   integer             :: ii, jj, idx, idx0a, idx0b

   idx0a = M * (M+1)
   idx0b = idx0a / 2
   do jj = 1, M
   do ii = jj, M
      idx = ii + (2*M-jj) * (jj-1) / 2
      FockMao_a(ii,jj) = RMM(idx+idx0a)
      FockMao_a(jj,ii) = RMM(idx+idx0a)
      FockMao_b(ii,jj) = RMM(idx+idx0b)
      FockMao_b(jj,ii) = RMM(idx+idx0b)
   enddo
   enddo

end subroutine rmmget_fockos_d
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
