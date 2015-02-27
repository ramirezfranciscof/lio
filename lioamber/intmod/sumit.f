!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine sumit(orbint,TermT,TermD,nvel,force,bij)
!--------------------------------------------------------------------!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       implicit none
       complex*16,intent(in)    :: orbint(3)
       complex*16,intent(in)    :: TermT,TermD
       real*8,intent(in)        :: nvel(3)

       complex*16,intent(inout) :: force(3)
       complex*16,intent(inout) :: bij

       integer :: kk

       do kk=1,3
         force(kk)=force(kk)+TermT*orbint(kk)
         force(kk)=force(kk)+TermD*conjg(orbint(kk))
         bij=bij+nvel(kk)*orbint(kk)
         print*,kk,force(kk)
       enddo

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
