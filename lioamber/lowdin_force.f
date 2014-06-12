!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE lowdin_force(ff,M,Fock,DSM,UInv,Rho,eigvec,eigval)
!------------------------------------------------------------------------------!
!
!
! Fuentes:
! (X) JChemPhys 130 224106 (2009)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       IMPLICIT NONE
       REAL*8,INTENT(OUT)               :: ff
       INTEGER,INTENT(IN)               :: M
       REAL*8,INTENT(IN),DIMENSION(M,M) :: Fock,DSM,UInv,Rho,eigvec
       REAL*8,INTENT(IN),DIMENSION(M)   :: eigval


       REAL*8,ALLOCATABLE :: csi(:,:),csj(:,:),fsi(:,:),fsj(:,:)
       REAL*8,ALLOCATABLE :: col(:,:),fil(:,:)
       REAL*8,ALLOCATABLE :: MAT(:,:),TSM(:,:),DUM(:,:),TRA(:,:)
       REAL*8             :: denom,num,TSN(1,1)
       INTEGER            :: ii,jj,kk
       INTEGER            :: TestID
!
!------------------------------------------------------------------------------!
!
       ff=0.0d0
       ALLOCATE(csi(M,1),csj(M,1),fsi(1,M),fsj(1,M))
       ALLOCATE(col(M,1),fil(1,M))
       ALLOCATE(MAT(M,M),TSM(M,M),DUM(M,M))
!
!----------------------------------------------------------!
!
       DO ii=1,M
       DO jj=1,M
         DUM(ii,jj)=0.0d0
       ENDDO
       ENDDO

       DO ii=1,M
       DO jj=1,M

         DO kk=1,M
           csi(kk,1)=eigvec(kk,ii)
           fsi(1,kk)=eigvec(kk,ii)
           csj(kk,1)=eigvec(kk,jj)
           fsj(1,kk)=eigvec(kk,jj)
         ENDDO

         col=MATMUL(DSM,csj)
         TSN=MATMUL(fsi,col)

         fil=MATMUL(TSN,fsj)
         denom=SQRT(eigval(ii))+SQRT(eigval(jj))
         fil=fil/denom
         DUM=DUM+MATMUL(csi,fil)

       ENDDO
       ENDDO

!       WRITE(6,*) 'DUM:'
!       WRITE(6,20) DUM(1,:)
!       WRITE(6,20) DUM(2,:)
!       WRITE(6,20) DUM(3,:)
!       WRITE(6,20) DUM(4,:)
!       WRITE(6,*) '---------------------'


       MAT=MATMUL(UInv,Rho)
       TSM=MATMUL(DUM,MAT)
       MAT=MATMUL(Fock,TSM)

       DO kk=1,M
         ff=ff+MAT(kk,kk)*2
       ENDDO
!
!----------------------------------------------------------!
!
       DEALLOCATE(MAT,TSM,DUM)
       DEALLOCATE(col,fil)
       DEALLOCATE(csi,csj,fsi,fsj)

20     FORMAT(4(E15.8,2x))

       RETURN;END SUBROUTINE
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
