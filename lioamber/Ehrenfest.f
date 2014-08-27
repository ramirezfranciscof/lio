!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE Ehrenfest(Energy,DipMom)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY:M,FockMat,RhoMat,dt,natom,nstep
       IMPLICIT NONE
       REAL*8,INTENT(INOUT) :: Energy,DipMom(3)

       INTEGER :: MM,ii,jj,kk,numit,maxit,ne_steps,jstep

       REAL*8  :: Enn,Esvn,Esve,Eelec,Ecoul,Exc,Efld
     > ,          crit,tolf,DipMod,trazarho,dt2

       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: FockNew,FockGuess,RhoReal
     > ,                                    Lmat,Linv,Umat,Uinv

       COMPLEX*8,ALLOCATABLE :: RhoNew(:,:)
       REAL*8,ALLOCATABLE    :: SMV(:),force(:,:)
       REAL*8,ALLOCATABLE    :: Matrix(:,:),Vector(:)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       call g2g_timer_start('Ehrenfest3dt')
       nstep=nstep+1
       ne_steps=10
       dt2=dt/ne_steps
       write(6,303)
       write(6,201) ' EHR - Starting nuclear step number  ',nstep
       write(6,302)
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
       IF (nstep.EQ.1) THEN
         call SCF(Energy,DipMom)
         call rmm_get_rho(M,RhoReal)

         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)

         call MatrixTransform(M,Umat,RhoReal,Lmat) ! RhoOA  => RhoON
         RhoMat=CMPLX(RhoReal)
         call MatrixTransform(M,Uinv,RhoReal,Linv) ! RhoON  => RhoOA 
       ELSE
         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)
       ENDIF
!
!
!
! NOT FIRST TIME? PROPAGA
!----------------------------------------------------------!
       DO jstep=1,ne_steps
         call rmm_put_rho(M,RhoReal)
         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockGuess)
         call MatrixTransform(M,Linv,FockMat,Uinv) ! FockOA => FockON

         call MagnusSC(M,RhoMat,FockMat,FockGuess,Uinv,Linv
     >   ,dt2,RhoNew,FockNew,Energy)
         RhoReal=real(RhoNew)
         RhoMat=RhoNew

         call rmm_get_energy(Eelec)
         Energy=Eelec+Ecoul+Exc+Enn+Esvn+Efld
         write(6,211) Eelec,Ecoul,Exc
         write(6,212) Enn,Efld
         write(6,213) Esvn,Esve
         write(6,202) '  Total Energy: ',Energy
         write(6,301)
       ENDDO
!
!
! ENERGY CALCULATIONS & OTHERS
!----------------------------------------------------------!
       call rmm_get_energy(Eelec)
       Energy=Eelec+Ecoul+Exc+Enn+Esvn+Efld
       write(6,201)
       write(6,211) Eelec,Ecoul,Exc
       write(6,212) Enn,Efld
       write(6,213) Esvn,Esve
       write(6,201)
       write(6,202) '  Total Energy: ',Energy
       write(6,302)

       call dip(DipMom(1),DipMom(2),DipMom(3))
       DipMod=sqrt(DipMom(1)**2+DipMom(2)**2+DipMom(3)**2)
       write(6,214) DipMod,DipMom

       call MatrixTransform(M,Umat,RhoReal,Lmat) ! RhoOA  => RhoON
       call rmm_exc_qmm(M,Uinv,RhoReal,FockMat,Linv)
       call MatrixTransform(M,Uinv,RhoReal,Linv) ! RhoON  => RhoOA
       force=0.0d0
       call intSG(force)
       write(6,215) force
!
!
!
! ROUTINE EXIT
!----------------------------------------------------------!
       deallocate(SMV,force)
       deallocate(RhoNew,RhoReal)
       deallocate(FockNew,FockGuess)
       deallocate(Lmat,Linv,Umat,Uinv)
       write(6,201) ' EHR - Done'
       call g2g_timer_stop('Ehrenfest3dt')
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 201   FORMAT(A,I6)
 202   FORMAT(A,F20.12)
 211   FORMAT('  Eelec: ',F12.5,3X,'Ecoul: ',F12.5,3X,'Exc: ',F12.5)
 212   FORMAT('  Enn:   ',F12.5,3X,'Efld:  ',F12.5)
 213   FORMAT('  Esvn:  ',F12.5,3X,'Esve:  ',F12.5)
 214   FORMAT('  Mod:   ',F12.5,' Vec:',3(1X,F12.5))
 215   FORMAT('  Force: ',3(1X,F12.5))
 301   FORMAT('----------------------------------------')
 302   FORMAT('-----------------------------------------------------
     >-----------------')
 303   FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     >%%%%%%%%%%%%%%%%% ')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
