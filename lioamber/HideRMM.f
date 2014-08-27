!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  FFR_HIDERMM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! For easy access to the matrixes in RMM vector
! rmm_get_rho(M,RhoMat)
! rmm_put_rho(M,RhoMat)
!
! rmm_get_sfm(M,SFMat)
! rmm_get_sfv(MM,SFVec)
! rmm_get_energy(Energy)
!
! rmm_put_qmm(M,QMat)
! rmm_exc_qmm(M,Uinv,Rho,Fock,Linv)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! When taken out of RMM packed storage (SP), rho matrix has doubled
! elements in non diagonal positions. This subroutine fixes that.
! When taken out of RMM packed storage (SP), rho matrix has doubled
! elements in non diagonal positions. So before packing rho, one has
! to mess it up.
!
!------------------------------------------------------------------------------!
       SUBROUTINE rmm_get_rho(M,RhoMat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in)         :: M
       real*8,intent(out)         :: RhoMat(M,M)
       integer                    :: ptr
!
       ptr=1
       call spunpack('L',M,RMM(ptr),RhoMat)
       call fixrho(M,RhoMat)
!
       RETURN;END SUBROUTINE

!------------------------------------------------------------------------------!
       SUBROUTINE rmm_put_rho(M,RhoMat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in)         :: M
       real*8,intent(in )         :: RhoMat(M,M)
       integer                    :: ptr
!
       ptr=1
       call messrho(M,RhoMat)
       call sprepack('L',M,RMM(ptr),RhoMat)
!
       RETURN;END SUBROUTINE

!------------------------------------------------------------------------------!
       SUBROUTINE fixrho(NM,Matrix)
!------------------------------------------------------------------------------!
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(OUT) :: Matrix(NM,NM)
       INTEGER            :: ii,jj
!
       DO ii=1,NM
       DO jj=1,NM
         IF (ii.ne.jj) THEN
           Matrix(ii,jj)=Matrix(ii,jj)/2
         ENDIF
       ENDDO
       ENDDO
!
       RETURN;END SUBROUTINE

!------------------------------------------------------------------------------!
       SUBROUTINE messrho(NM,Matrix)
!------------------------------------------------------------------------------!
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: NM
       REAL*8,INTENT(OUT) :: Matrix(NM,NM)
       INTEGER            :: ii,jj
!
       DO ii=1,NM
       DO jj=1,NM
         IF (ii.ne.jj) THEN
           Matrix(ii,jj)=2*Matrix(ii,jj)
         ENDIF
       ENDDO
       ENDDO
!
       RETURN;END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE rmm_get_sfm(M,SFMat)
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in)         :: M
       real*8,intent(out)         :: SFMat(M,M)
       integer                    :: ptr
!------------------------------------------------------------------------------!
       ptr=1+M*(M+1)
       call spunpack('L',M,RMM(ptr),SFMat)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE rmm_get_sfv(MM,SFVec)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in)         :: MM
       real*8,intent(out)         :: SFVec(MM)
       integer                    :: kk,ptr

       ptr=1+2*MM
       do kk=1,MM
         SFVec(kk)=RMM(ptr-1+kk)
       enddo
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE rmm_get_energy(Energy)
       use garcha_mod, only:M,Md,RMM
       implicit none
       real*8,intent(out)         :: Energy
       integer                    :: MM,MMd,kk,ptr1,ptr2
!------------------------------------------------------------------------------!
       MM=M*(M+1)/2
       MMd=Md*(Md+1)/2
       ptr1=1               ! M1 rho
       ptr2=ptr1+3*MM+2*MMd ! M11

       Energy=0.0d0
       do kk=1,MM
         Energy=Energy+RMM(ptr1-1+kk)*RMM(ptr2-1+kk)
       enddo
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE rmm_put_qmm(M,QMat)
       use garcha_mod, only:Md,RMM
       implicit none
       integer,intent(in)         :: M
       real*8,intent(in )         :: QMat(M,M)
       integer                    :: ptr,idx,ii,jj,MM,MMd
!------------------------------------------------------------------------------!
       MM=M*(M+1)/2
       MMd=Md*(Md+1)/2
       ptr=1+M+4*MM+2*MMd
       do ii=1,M
       do jj=1,ii
           idx=(ptr-1)+ii+(2*M-jj)*(jj-1)/2
           RMM(idx)=QMat(ii,jj)+QMat(jj,ii)
           IF (ii.EQ.jj) RMM(idx)=RMM(idx)/2
       enddo
       enddo
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE rmm_exc_qmm(M,Uinv,Rho,Fock,Linv)
       implicit none
       integer,intent(in) :: M
       real*8,intent(in)  :: Uinv(M,M),Rho(M,M),Fock(M,M),Linv(M,M)
       real*8,allocatable :: Qaux(:,:)
!------------------------------------------------------------------------------!
       allocate(Qaux(M,M))
       call prepQM(M,Qaux,UInv,Rho,Fock,LInv)
       call rmm_put_qmm(M,Qaux)
       deallocate(Qaux)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
