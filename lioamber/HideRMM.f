!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  HideRMM
!------------------------------------------------------------------------------!
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
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Rho - The Density Matrix
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! In the RMM vector, Rho is stored in packed storage but
! each non diagonal position has dobule the real value.
! So before putting a density matrix in RMM, non diagonal
! positions need to be x2 (messrho) and when taking it
! out of RMM, the non diagonal positions need to be
! divided by two (fixrho).
!
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
       return;end subroutine

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
       return;end subroutine

!------------------------------------------------------------------------------!
       SUBROUTINE fixrho(M,Matrix)
!------------------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: M
       real*8,intent(out) :: Matrix(M,M)
       integer            :: ii,jj
!
       do ii=1,M
       do jj=1,M
         if (ii.ne.jj) then
           Matrix(ii,jj)=Matrix(ii,jj)/2
         endif
       enddo
       enddo
!
       return;end subroutine

!------------------------------------------------------------------------------!
       SUBROUTINE messrho(M,Matrix)
!------------------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: M
       real*8,intent(out) :: Matrix(M,M)
       integer            :: ii,jj
!
       do ii=1,M
       do jj=1,M
         if (ii.ne.jj) then
           Matrix(ii,jj)=2*Matrix(ii,jj)
         endif
       enddo
       enddo
!
       return;end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Fock, S and Energy
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! These subroutines take and put matrixes in the spaces
! of RMM that are reserved for S and Fock matrixes.
! Additionally, it includes one subroutine for the
! re-calculation of energy.
!
!
!------------------------------------------------------------------------------!
       SUBROUTINE rmm_get_sfm(M,SFMat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:RMM
       implicit none
       integer,intent(in)         :: M
       real*8,intent(out)         :: SFMat(M,M)
       integer                    :: ptr

       ptr=1+M*(M+1)
       call spunpack('L',M,RMM(ptr),SFMat)

       return;end subroutine

!------------------------------------------------------------------------------!
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

       return;end subroutine

!------------------------------------------------------------------------------!
       SUBROUTINE rmm_get_energy(Energy)
!------------------------------------------------------------------------------!
       use garcha_mod, only:M,Md,RMM
       implicit none
       real*8,intent(out)         :: Energy
       integer                    :: MM,MMd,kk,ptr1,ptr2

       MM=M*(M+1)/2
       MMd=Md*(Md+1)/2
       ptr1=1               ! M1 rho
       ptr2=ptr1+3*MM+2*MMd ! M11

       Energy=0.0d0
       do kk=1,MM
         Energy=Energy+RMM(ptr1-1+kk)*RMM(ptr2-1+kk)
       enddo

       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! QMM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 
! Subroutine to configure the matrix that is used to
! calculate forces in intSG. SCF sets the matrix Q that
! is referenced in Szabo, page 441 (or an analagous one)
! This subrutines substitutes it for one apropriate for
! the calculation of the forces in the case of excited
! electronic states, according to what is indicated in
!
! JCP 114, 9758 (2001); doi: 10.1063/1.1372182
!
! for the case of cholesky diagonalization of S.
!
!------------------------------------------------------------------------------!
       SUBROUTINE rmm_put_qmm(M,QMat)
!------------------------------------------------------------------------------!
       use garcha_mod, only:Md,RMM
       implicit none
       integer,intent(in)         :: M
       real*8,intent(in )         :: QMat(M,M)
       integer                    :: ptr,idx,ii,jj,MM,MMd

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

       return;end subroutine

!------------------------------------------------------------------------------!
       SUBROUTINE rmm_exc_qmm(M,Uinv,Rho,Fock,Linv)
!------------------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: M
       real*8,intent(in)  :: Uinv(M,M),Rho(M,M),Fock(M,M),Linv(M,M)
       real*8,allocatable :: Qaux(:,:)

       allocate(Qaux(M,M))
       call prepQM(M,Qaux,UInv,Rho,Fock,LInv)
       call rmm_put_qmm(M,Qaux)
       deallocate(Qaux)

       return;end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
