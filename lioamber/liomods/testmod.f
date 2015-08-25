!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module testmod
       implicit none
       real*8,allocatable :: DSX(:,:,:),DSY(:,:,:),DSZ(:,:,:)
       real*8,allocatable :: DSXc(:,:,:),DSYc(:,:,:),DSZc(:,:,:)
       real*8,allocatable :: DSXt(:,:,:),DSYt(:,:,:),DSZt(:,:,:)

       contains
!------------------------------------------------------------------------------!
       subroutine testinit()
       use garcha_mod
       implicit none
       allocate(DSX(M,M,natom),DSY(M,M,natom),DSZ(M,M,natom))
       allocate(DSXc(M,M,natom),DSYc(M,M,natom),DSZc(M,M,natom))
       allocate(DSXt(M,M,natom),DSYt(M,M,natom),DSZt(M,M,natom))

       return;end subroutine
       end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

