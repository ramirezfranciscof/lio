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
!       WRITE(6,*) "-----------------------------------"
!       WRITE(6,*) "FORCES:"
!
       call g2g_timer_start('int1G')
       call int1G(ff)
       call g2g_timer_stop('int1G')
!       WRITE(6,666) 'hola int1g',ff(1,:)
!
       call g2g_timer_start('intSG')
       call intSG(ff)
       call g2g_timer_stop('intSG')
!       WRITE(6,666) 'hola intSG',ff(1,:)
!
       call g2g_timer_start('int3G')
       call int3G(ff,.true.)
       call g2g_timer_stop('int3G')
!       WRITE(6,666) 'hola int3G',ff(1,:)
!       WRITE(6,*) "-----------------------------------"
       ff=0 ! SOLO PARA DEBUG

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
 666   FORMAT(A,3(2X,E15.5))
       end subroutine dft_get_qm_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
