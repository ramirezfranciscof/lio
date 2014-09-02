!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE mulliken(OutID)
!------------------------------------------------------------------------------!
!
! DESCRIPTION PENDING
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:M,Md,RMM,natom,Iz,Nuc
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: OutID

       INTEGER            :: ii,jj,idx,MM,M1,M2,M5
       REAL*8             :: Energy,t0
       REAL*8,ALLOCATABLE :: q(:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
! SETUP
!----------------------------------------------------------!
       PRINT*,'AVER'
       call g2g_timer_start('mulliken')
       call int1(Energy)
       allocate(q(natom))
       do ii=1,natom
         q(ii)=Iz(ii)
       enddo

       MM=M*(M+1)/2
       M1=1
       M2=2*M
       M5=M1+MM*2
!
!
! CALC MULLIKEN
!----------------------------------------------------------!
       do ii=1,M
         do jj=1,ii-1
           idx=ii+(M2-jj)*(jj-1)/2
           t0=RMM(idx)*RMM(M5+idx-1)/2.D0
           q(Nuc(ii))=q(Nuc(ii))-t0
         enddo

         idx=ii+(M2-ii)*(ii-1)/2
         t0=RMM(idx)*RMM(M5+idx-1)
         q(Nuc(ii))=q(Nuc(ii))-t0

         do jj=ii+1,M
           idx=jj+(M2-ii)*(ii-1)/2
           t0=RMM(idx)*RMM(M5+idx-1)/2.D0
           q(Nuc(ii))=q(Nuc(ii))-t0
         enddo
       enddo
!
!
! WRITE OUTPUT
!----------------------------------------------------------!
       write(OutID,300)
       write(OutID,200) 'MULLIKEN POPULATION ANALYSIS'
       write(OutID,200)
       write(OutID,201) 'ATOM #','ATOM TYPE','POPULATION'
       do ii=1,natom
         write(OutID,202) ii,Iz(ii),q(ii)
       enddo
       write(OutID,200)
!
! ROUTINE EXIT
!----------------------------------------------------------!
       deallocate(q)
       call g2g_timer_stop('mulliken')
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 200   FORMAT(A)
 201   FORMAT(A,4x,A,4x,A)
 202   FORMAT(I3,9X,I3,6X,F14.8)
 300   FORMAT('##################################################
     >##########')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
