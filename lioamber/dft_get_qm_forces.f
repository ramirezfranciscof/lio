!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine dft_get_qm_forces(dxyzqm)
!------------------------------------------------------------------------------!
!
! This subroutine calculates the derivatives of E with
! respect to the atomic coordinates of the quantum part
! of the system.
!
! Note that the forces will be minus this derivative;
! that needs to be considered by the calling routine.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:natom
c       use qmmm_module, only : qmmm_struct
       IMPLICIT NONE

       REAL*8,INTENT(INOUT)              :: dxyzqm(3,natom)
       REAL*8,DIMENSION(:,:),ALLOCATABLE :: ff
       REAL*8                            :: ftot(3),factor
       INTEGER                           :: ii,jj
!
!------------------------------------------------------------------------------!
!
       allocate(ff(natom,3))
       ff=0
!
       call g2g_timer_start('int1G')
       call int1G(ff)
       call g2g_timer_stop('int1G')
!
       call g2g_timer_start('intSG')
       call intSG(ff)
       call g2g_timer_stop('intSG')
c       write(77,*) ff
!
       call g2g_timer_start('int3G')
       call int3G(ff,.true.)
       call g2g_timer_stop('int3G')
!
c       factor=627.509391D0/0.5291772108D0
       factor=1.D0
       do ii=1,natom 
       do jj=1,3
         dxyzqm(jj,ii)=ff(ii,jj)*factor
       enddo
       enddo
!
!------------------------------------------------------------------------------!
       deallocate (ff)
       end subroutine dft_get_qm_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
