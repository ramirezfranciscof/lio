!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmget_core_r( CoreMao )

   use garcha_mod, only: M, Md, RMM
   implicit none
   real*4,intent(out) :: CoreMao(M,M)
   integer            :: MM, MMd, ii, jj, idx, idx0

   MM  = M  * (M+1)  / 2
   MMd = Md * (Md+1) / 2
   idx0 = 3*MM + 2*MMd
   do jj = 1, M
   do ii = jj, M
      idx = ii + (2*M-jj) * (jj-1) / 2 + idx0
      CoreMao(ii,jj) = RMM(idx)
      CoreMao(jj,ii) = RMM(idx)
   enddo
   enddo

end subroutine rmmget_core_r
!--------------------------------------------------------------------!
subroutine rmmget_core_d( CoreMao )

   use garcha_mod, only: M, Md, RMM
   implicit none
   real*8,intent(out) :: CoreMao(M,M)
   integer            :: MM, MMd, ii, jj, idx, idx0

   MM  = M  * (M+1)  / 2
   MMd = Md * (Md+1) / 2
   idx0 = 3*MM + 2*MMd
   do jj = 1, M
   do ii = jj, M
      idx = ii + (2*M-jj) * (jj-1) / 2 + idx0
      CoreMao(ii,jj) = RMM(idx)
      CoreMao(jj,ii) = RMM(idx)
   enddo
   enddo

end subroutine rmmget_core_d
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
