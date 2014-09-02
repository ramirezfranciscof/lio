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
! (a more suitable name for this subroutine would have
! been "dft_get_qm_Ederivatives").
!
!------------------------------------------------------------------------------!
       USE garcha_mod,ONLY:natom
       IMPLICIT NONE
       REAL*8,INTENT(INOUT)              :: dxyzqm(3,natom)
       REAL*8,DIMENSION(:,:),ALLOCATABLE :: force
       REAL*8,DIMENSION(:,:),ALLOCATABLE :: force1G,forceSG,force3G
       REAL*8                            :: factor
       INTEGER                           :: ii,jj
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       allocate(force(natom,3))
       allocate(force1G(natom,3),forceSG(natom,3),force3G(natom,3))
       open(unit=200,file='forces.dat',position="append")
       force   = 0.0d0
       force1G = 0.0d0
       forceSG = 0.0d0
       force3G = 0.0d0

       call g2g_timer_start('int1G')
       call int1G(force1G)
       call g2g_timer_stop('int1G')
       call g2g_timer_start('intSG')
       call intSG(forceSG)
       call g2g_timer_stop('intSG')
       call g2g_timer_start('int3G')
       call int3G(force3G,.true.)
       call g2g_timer_stop('int3G')
       force=force1G+forceSG+force3G

!       call print_forceset(200,natom,force1G,'F1G:')
       call print_forceset(200,natom,forceSG,'FSG:')
!       call print_forceset(200,natom,force3G,'F3G:')

       factor=1.0d0
c       factor=627.509391D0/0.5291772108D0
       do ii=1,natom;do jj=1,3
         dxyzqm(jj,ii)=force(ii,jj)*factor
       enddo;enddo

!------------------------------------------------------------------------------!
       close(200)
       deallocate(force,force1G,forceSG,force3G)
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine print_forceset(nunit,natoms,force,descript)
       integer,intent(in) :: nunit,natoms
       real*8,intent(in)  :: force(natoms,3)
       character(len=4)   :: descript
       integer            :: ii,jj
!------------------------------------------------------------------------------!
       do ii=1,natoms
       write(nunit,200) descript,ii,(force(ii,jj),jj=1,3)
       enddo
       write(nunit,'(A)')
     > '------------------------------------------------------------'
!------------------------------------------------------------------------------!
 200   format(1X,A4,1X,I2,3(2X,E15.7),1X)
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
