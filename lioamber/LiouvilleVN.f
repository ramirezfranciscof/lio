!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE LiouvilleVN(Energy,DipMom)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:M,FockMat,RhoMat,dt,natom,istep
       IMPLICIT NONE
       REAL*8,INTENT(INOUT) :: Energy,DipMom(3)

       INTEGER :: MM,ii,jj,kk,numit,maxit

       REAL*8  :: Enn,Esvn,Esve,Eelec,Ecoul,Exc,Efld
     > ,          crit,tolf,DipMod,trazarho

       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: FockNew,FockGuess,RhoReal
     > ,                                    Lmat,Linv,Umat,Uinv

       COMPLEX*8,ALLOCATABLE :: RhoNew(:,:)
       REAL*8,ALLOCATABLE    :: SMV(:),force(:,:)
       REAL*8,ALLOCATABLE    :: Matrix(:,:),Vector(:)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       call g2g_timer_start('LiouvilleVN')
       istep=istep+1
       write(6,301)
       write(6,200) 'LvN - Starting step: ',istep
       write(6,300)

!
!
!
! SETUP
!----------------------------------------------------------!
       MM=M*(M+1)/2
       allocate(Matrix(M,M),Vector(MM))
       allocate(SMV(MM),force(natom,3))
       allocate(RhoNew(M,M),RhoReal(M,M))
       allocate(FockNew(M,M),FockGuess(M,M))
       allocate(Lmat(M,M),Linv(M,M),Umat(M,M),Uinv(M,M))
!
!
!
! FIRST TIME? CALL SCF AND SETUP FOCKMAT AND RHOMAT
!----------------------------------------------------------!
       IF (istep.EQ.1) THEN
         call SCF(Energy,DipMom)
         call rmm_get_rho(M,RhoReal)

         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)

         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockMat)

         call MatrixTransform(M,Linv,FockMat,Uinv)
         call MatrixTransform(M,Umat,RhoReal,Lmat)
         RhoMat=CMPLX(RhoReal)
       ENDIF
!
!
!
! NOT FIRST TIME? DO LIOUVILLE (SC LOOP)
!----------------------------------------------------------!
       IF (istep.GT.1) THEN
         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV) ! ¿Efld?
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)

         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockGuess)
         call MatrixTransform(M,Linv,FockGuess,Uinv)

         call magnus_sc(M,RhoMat,FockMat,FockGuess,Uinv,Linv,dt,
     >        RhoNew,FockNew,Energy)
         RhoMat=RhoNew
         RhoReal=real(RhoNew)
         FockMat=FockNew
       ENDIF
!
!
! ENERGY CALCULATIONS & OTHERS
!----------------------------------------------------------!
       call rmm_get_energy(Eelec)
       Energy=Eelec+Ecoul+Exc+Enn+Esvn ! ¿Efld?
       write(6,201) Eelec,Ecoul,Exc
       write(6,202) Enn,Efld
       write(6,203) Esvn,Esve
       write(6,200)
       write(6,200) '  Total Energy: ',Energy
       write(6,300)

       call dip(DipMom(1),DipMom(2),DipMom(3))
       DipMod=sqrt(DipMom(1)**2+DipMom(2)**2+DipMom(3)**2)
       write(6,205) DipMod,DipMom
       call MullikenPop(85)

       force=0.0d0
       call intSG(force)
       WRITE(777,*) '-------------------'
       WRITE(777,*) force
       call rmm_exc_qmm(M,Uinv,RhoReal,FockMat,Linv)
       force=0.0d0
       call intSG(force)
       WRITE(777,*) force
!
!
!
! ROUTINE EXIT
!----------------------------------------------------------!
       deallocate(SMV,force)
       deallocate(RhoNew,RhoReal)
       deallocate(FockNew,FockGuess)
       deallocate(Lmat,Linv,Umat,Uinv)

       write(6,300)
       write(6,200) 'LvN - Done. '
       write(6,200)
       call g2g_timer_stop('LiouvilleVN')
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 200   FORMAT(A,F20.12)
 201   FORMAT('  Eelec: ',F12.5,3X,'Ecoul: ',F12.5,3X,'Exc: ',F12.5)
 202   FORMAT('  Enn:   ',F12.5,3X,'Efld:  ',F12.5)
 203   FORMAT('  Esvn:  ',F12.5,3X,'Esve:  ',F12.5)
 205   FORMAT('  Mod: ',F12.5,' Vec:',3(1X,F12.5))
 300   FORMAT('-----------------------------------------------------
     >-----------------')
 301   FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     >%%%%%%%%%%%%%%%%% ')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

