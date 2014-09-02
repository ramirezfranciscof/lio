!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE LiouvilleVN(Energy,DipMom)
!------------------------------------------------------------------------------!
!
! DESCRIPTION PENDING
!
!------------------------------------------------------------------------------!
       USE garcha_mod, ONLY:M,FockMat,RhoMat,dt,natom,nstep
       IMPLICIT NONE
       REAL*8,INTENT(INOUT) :: Energy,DipMom(3)

       INTEGER :: MM,ii,jj,kk,numit,maxit
       REAL*8  :: 
     > Enn,Esvn,Esve,Eelec,Ecoul,Exc,Efld,crit,tolf,DipMod,trazarho

       COMPLEX*16,ALLOCATABLE :: RhoNew(:,:)
       REAL*8,ALLOCATABLE     :: SMV(:)
       REAL*8,ALLOCATABLE,DIMENSION(:,:) ::
     > FockNew,FockGuess,RhoReal,Lmat,Linv,Umat,Uinv

       REAL*8,ALLOCATABLE     :: Matrix(:,:),Vector(:)
       REAL*8,ALLOCATABLE     :: force(:,:),forceb(:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
! SETUP
!----------------------------------------------------------!
       call g2g_timer_start('LiouvilleVN')
       nstep=nstep+1
       write(6,416)
       write(6,200) 'LvN - Starting step: ',nstep
       write(6,416)

       MM=M*(M+1)/2
       allocate(SMV(MM))
       allocate(RhoNew(M,M),RhoReal(M,M))
       allocate(FockNew(M,M),FockGuess(M,M))
       allocate(Lmat(M,M),Linv(M,M),Umat(M,M),Uinv(M,M))

       allocate(Matrix(M,M),Vector(MM))
       allocate(force(natom,3),forceb(natom,3))
!
!
!
! FIRST TIME: CALL SCF AND SETUP FOCKMAT AND RHOMAT
!----------------------------------------------------------!
       IF (nstep.EQ.1) THEN
         call SCF(Energy,DipMom)
         call rmm_get_rho(M,RhoReal)
         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockMat)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)
         call MatrixTransform(M,Linv,FockMat,Uinv) ! FockOA => FockOM
         call MatrixTransform(M,Umat,RhoReal,Lmat) ! RhoOA  => RhoOM

         IF (.FALSE.) THEN
         call SCF(Energy,DipMom)
         call rmm_get_rho(M,RhoReal)

         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)
         call rmm_get_sfm(M,Matrix)

         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockMat)
         call MatrixTransform(M,Linv,FockMat,Uinv) ! FockOA => FockOM
         call MatrixTransform(M,Umat,RhoReal,Lmat) ! RhoOA  => RhoOM
         RhoMat=cmplx(RhoReal)
         CALL ffr_print_conmut(666,M,RhoReal,FockMat)
         write(666,*) '--------'
         ENDIF
       ENDIF
!
!
!
! NOT FIRST TIME: DO LIOUVILLE (SC LOOP)
!----------------------------------------------------------!
       IF (nstep.GT.1) THEN
         call ffr_reset_mem()
         call SCF(Energy,DipMom)
         call rmm_get_rho(M,RhoReal)
         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockMat)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)
         call MatrixTransform(M,Linv,FockMat,Uinv) ! FockOA => FockOM
         call MatrixTransform(M,Umat,RhoReal,Lmat) ! RhoOA  => RhoOM

         IF (.FALSE.) THEN
         call ffr_reset_mem()
         call ffr_setup_S(Enn,Esve,Esvn)
         call rmm_get_sfv(MM,SMV)
         call cholesky_sdiag(M,SMV,Lmat,Linv,Umat,Uinv)
         call rmm_get_sfm(M,Matrix)

         call ffr_setup_F(Ecoul,Exc,Eelec,Efld)
         call rmm_get_sfm(M,FockGuess)
         call MatrixTransform(M,Linv,FockGuess,Uinv) ! FockOA => FockOM
         CALL ffr_print_conmut(666,M,RhoReal,FockMat)
         write(666,*) '--------'
         call MagnusSC
     >   (M,RhoMat,FockMat,FockGuess,Uinv,Linv,dt,RhoNew,FockNew,Energy)

         RhoMat=RhoNew
         RhoReal=real(RhoNew)
         FockMat=FockNew
         ENDIF
       ENDIF
!
!
!
! PREPARATION FOR FORCE CALCULATIONS
!----------------------------------------------------------!
!       call rmm_exc_qmm(M,Uinv,RhoReal,FockMat,Linv) ! NECESITA OM
       IF (.FALSE.) THEN
       call MatrixTransform(M,Uinv,RhoReal,Linv) ! RhoON  => RhoOA
       call rmm_put_rho(M,RhoReal)
       call MatrixTransform(M,Umat,RhoReal,Lmat) ! RhoOA  => RhoOM
       ENDIF
!
!
!
! ENERGY CALCULATIONS & OTHERS
!----------------------------------------------------------!
       IF (.FALSE.) THEN
       call rmm_get_energy(Eelec)
       call LiouvilleVN_printE(6,Eelec,Ecoul,Exc,Enn,Efld,Esvn,Esve)
       Energy=Eelec+Ecoul+Exc+Enn+Esvn+Efld

       call dip(DipMom(1),DipMom(2),DipMom(3))
       call LiouvilleVN_printD(6,DipMom)
!       call mulliken(85) ! ESTO NO FUNCA
       ENDIF
!
!
!
! ROUTINE EXIT
!----------------------------------------------------------!
       deallocate(Matrix,Vector,force,forceb)

       deallocate(SMV)
       deallocate(RhoNew,RhoReal)
       deallocate(FockNew,FockGuess)
       deallocate(Lmat,Linv,Umat,Uinv)

       write(6,200);write(6,406)
       write(6,200) 'LvN - Done. '
       write(6,200);write(6,200)
       call g2g_timer_stop('LiouvilleVN')
!
!
!
!------------------------------------------------------------------------------!
 200   FORMAT(A,I10)
 406   FORMAT('----------------------------------------',
     >        '--------------------')
 416   FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',
     >        '%%%%%%%%%%%%%%%%%%%%')
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE LiouvilleVN_printE
     > (nunit,Eelec,Ecoul,Exc,Enn,Efld,Esvn,Esve)
       implicit none
       integer,intent(in) :: nunit
       real*8,intent(in)  :: Eelec,Ecoul,Exc,Enn,Efld,Esvn,Esve
       real*8             :: Etot
!------------------------------------------------------------------------------!
       Etot=Eelec+Ecoul+Exc+Enn+Esvn+Efld
       write(nunit,201)
!       write(nunit,405)
       write(nunit,201) 'Total Energy:',Etot
       write(nunit,305)
       write(nunit,202) 'Eelec:',Eelec,'Ecoul: ',Ecoul
       write(nunit,202) 'Exc:  ',Exc,'Enn:   ',Enn
       write(nunit,202) 'Esvn: ',Esvn,'Esve: ',Esve
       write(nunit,202) 'Efld: ',Efld
       write(nunit,305)
!       write(nunit,405)
       write(nunit,201)
!------------------------------------------------------------------------------!
  201  FORMAT(2X,A,F15.8)
  202  FORMAT(2(2X,A6,1X,F14.8,2X))
  301  FORMAT('----------')
  305  FORMAT('--------------------------------------------------')
  405  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE LiouvilleVN_printD
     > (nunit,DipMom)
       implicit none
       integer,intent(in) :: nunit
       real*8,intent(in)  :: DipMom(3)
       real*8             :: DipTot
!------------------------------------------------------------------------------!
       DipTot=sqrt(DipMom(1)**2+DipMom(2)**2+DipMom(3)**2)
!       write(nunit,201)
       write(nunit,201) 'Dipole Moment:',DipTot
       write(nunit,202) 
     > 'DMX:',DipMom(1),'DMY:',DipMom(2),'DMZ',DipMom(3)
!       write(nunit,201)
!------------------------------------------------------------------------------!
  201  FORMAT(2X,A,F15.6)
  202  FORMAT(3(2X,A4,1X,F12.7,1X))
  301  FORMAT('----------')
  305  FORMAT('--------------------------------------------------')
  405  FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
