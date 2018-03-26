!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmput_densos_r( DensMao_a, DensMao_b )
   use garcha_mod, only: M, RMM, rhoalpha, rhobeta
   implicit none
   real*4, intent(in)     :: DensMao_a(M,M)
   real*4, intent(in)     :: DensMao_b(M,M)
   integer                :: ii, jj, idx

   do jj = 1, M
      idx = jj + (2*M-jj) * (jj-1) / 2
      rhoalpha(idx) = DensMao_a(jj,jj)
      rhobeta(idx)  = DensMao_b(jj,jj)
      RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      do ii = jj+1, M
         idx = ii + (2*M-jj) * (jj-1) / 2
         rhoalpha(idx) = DensMao_a(jj,jj) * 2
         rhobeta(idx)  = DensMao_b(jj,jj) * 2
         RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      enddo
   enddo

end subroutine rmmput_densos_r
!
!
!--------------------------------------------------------------------!
subroutine rmmput_densos_d( DensMao_a, DensMao_b )
   use garcha_mod, only: M, RMM, rhoalpha, rhobeta
   implicit none
   real*8, intent(in)     :: DensMao_a(M,M)
   real*8, intent(in)     :: DensMao_b(M,M)
   integer                :: ii, jj, idx

   do jj = 1, M
      idx = jj + (2*M-jj) * (jj-1) / 2
      rhoalpha(idx) = DensMao_a(jj,jj)
      rhobeta(idx)  = DensMao_b(jj,jj)
      RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      do ii = jj+1, M
         idx = ii + (2*M-jj) * (jj-1) / 2
         rhoalpha(idx) = DensMao_a(jj,jj) * 2
         rhobeta(idx)  = DensMao_b(jj,jj) * 2
         RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      enddo
   enddo

end subroutine rmmput_densos_d
!
!
!--------------------------------------------------------------------!
subroutine rmmput_densos_c( DensMao_a, DensMao_b )
   use garcha_mod, only: M, RMM, rhoalpha, rhobeta
   implicit none
   complex*8, intent(in)  :: DensMao_a(M,M)
   complex*8, intent(in)  :: DensMao_b(M,M)
   integer                :: ii, jj, idx

   do jj = 1, M
      idx = jj + (2*M-jj) * (jj-1) / 2
      rhoalpha(idx) = REAL(DensMao_a(jj,jj))
      rhobeta(idx)  = REAL(DensMao_b(jj,jj))
      RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      do ii = jj+1, M
         idx = ii + (2*M-jj) * (jj-1) / 2
         rhoalpha(idx) = REAL(DensMao_a(jj,jj)) * 2
         rhobeta(idx)  = REAL(DensMao_b(jj,jj)) * 2
         RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      enddo
   enddo

end subroutine rmmput_densos_c
!
!
!--------------------------------------------------------------------!
subroutine rmmput_densos_z( DensMao_a, DensMao_b )
   use garcha_mod, only: M, RMM, rhoalpha, rhobeta
   implicit none
   complex*16, intent(in) :: DensMao_a(M,M)
   complex*16, intent(in) :: DensMao_b(M,M)
   integer                :: ii, jj, idx

   do jj = 1, M
      idx = jj + (2*M-jj) * (jj-1) / 2
      rhoalpha(idx) = REAL(DensMao_a(jj,jj))
      rhobeta(idx)  = REAL(DensMao_b(jj,jj))
      RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      do ii = jj+1, M
         idx = ii + (2*M-jj) * (jj-1) / 2
         rhoalpha(idx) = REAL(DensMao_a(jj,jj)) * 2
         rhobeta(idx)  = REAL(DensMao_b(jj,jj)) * 2
         RMM(idx) = rhoalpha(idx) + rhobeta(idx)
      enddo
   enddo

end subroutine rmmput_densos_z
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
