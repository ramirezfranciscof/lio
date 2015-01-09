!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor(F1a,F1b,FON,rho2,devPtrX,factorial,
     > Fxx,Fyy,Fzz,g)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!Predictor-Corrector Cheng, V.Vooris.PhysRevB.2006.74.155112
! Esta rutina recibe: F1a,F1b,rho2
! Tira: F5 = F(t+(deltat/2))      
       use garcha_mod
       REAL*8,intent(inout) :: F1a(M,M),F1b(M,M),FON(M,M)
       integer*8,intent(in) :: devPtrX
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       integer :: i,j,k,kk,stat
       real*8 :: E2, tdstep1
      external CUBLAS_INIT, CUBLAS_SHUTDOWN
      integer CUBLAS_INIT
      REAL*8,intent(in) :: factorial(NBCH)
      REAL*8,intent(in) :: g,Fxx,Fyy,Fzz
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: rho2(M,M)
       COMPLEX*8,allocatable :: rho4(:,:),rho2t(:,:)
#else
       COMPLEX*16, intent(in) :: rho2(M,M)
       COMPLEX*16,allocatable :: rho4(:,:),rho2t(:,:)
#endif
!-----------------------------------------------------------------------------n
      ALLOCATE(rho4(M,M),rho2t(M,M),F3(M,M),FBA(M,M))
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -predictor"
        call CUBLAS_SHUTDOWN
        stop
      endif
      M2=2*M
      MM=M*(M+1)/2
c first i
      M1=1
c now Fold
      M3=M1+MM
c now S, F also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
c Least squares
      M17=M15+MM
c vectors of MO
      M18=M17+MMd
c weights (in case of using option )
      M19=M18+M*NCO
c
* RAM storage of two-electron integrals (if MEMO=T)
      M20 = M19 + natom*50*Nang
c Initializations/Defaults
c xmm es la primer matriz de (M,M) en el 
!------------------------------------------------------------------------------!
! Codigo del predictor:
!------------------------------------------------------------------------------!
! tdstep predictor es 0.5 tdstep magnum
       tdstep1=tdstep*0.5
! Paso1: Con las matrices pasadas F1a y F1b extrapolamos a F3----> Extrapolacion
       F3=(7.D0/4.D0)*F1b-(3.D0/4.D0)*F1a
!       F3=1.750D0*F1b-0.750D0*F1a
! Paso2: Usando H3, la matriz densidad rho2 es propagada a rho4----> Prediccion
       rho2t=rho2
       call cumagnusfac(F3,rho2,rho4,M,NBCH,tdstep1,factorial)
! Paso3: Escribimos rho4 en el RMM para poder obtener F5 en el siguiente paso.
!       call g2g_timer_start('cumatmul_predictor')
!       call cumxp(rho4,devPtrX,rho2t,M)
!       call cumpxt(rho2t,devPtrX,rho2t,M)
!       call g2g_timer_stop('cumatmul_predictor')
       call complex_rho_on_to_ao(rho4,devPtrX,rho2t,M)
!       do j=1,M
!       do k=j,M
!         if(j.eq.k) then
!           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))
!         else
!           RMM(k+(M2-j)*(j-1)/2)=REAL(rho2t(j,k))*2
!         endif
!       enddo
!       enddo
      call sprepack_ctr('L',M,RMM,rho2t)
! Paso4: La matriz densidad 4 es usada para calcular F5------> Corrector
      call int3lu(E2)
      call g2g_solve_groups(0,Ex,0)
      stat=CUBLAS_INIT()
      if (stat.NE.0) then
        write(*,*) "initialization failed -predictor"
        call CUBLAS_SHUTDOWN
        stop
      endif
      if (field) then
         write(*,*) 'FIELD PREDICTOR'
         call dip2(g,Fxx,Fyy,Fzz)
      endif
!------------------------------------------------------------------------------!
! Escritura de fock cuadrada
!------------------------------------------------------------------------------!
! Parte Inferior Izquierda (con diagonal)
!       do j=1,M
!         do k=1,j
!           FBA(j,k)=RMM(M5+j+(M2-k)*(k-1)/2-1)
!         enddo
! Parte Superior Derecha (sin diagonal)
!         do k=j+1,M
!           FBA(j,k)=RMM(M5+k+(M2-j)*(j-1)/2-1)
!         enddo
!       enddo
       call spunpack('L',M,RMM(M5),FBA)
! Ahora tenemos F5 transformada en base de ON y en su forma cuadrada
!       call g2g_timer_start('cumatmul2_predictor')
!       call cumxtf(FBA,devPtrX,FON,M)
!       call cumfx(FON,DevPtrX,FON,M)
!       call g2g_timer_stop('cumatmul2_predictor')
       call fock_ao_to_on(FBA,devPtrX,FON,M)
       call CUBLAS_SHUTDOWN
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine cupredictor_op(F1a_a,F1b_a,F1a_b,F1b_b,FON_a,FON_b,
     > rho2_a,rho2_b,factorial,devPtrX,Fxx,Fyy,Fzz,g)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This routine recives: F1a,F1b,rho2
! And gives: F5 = F(t+(deltat/2))      
       use garcha_mod
       REAL*8,intent(inout) :: F1a_a(M,M),F1b_a(M,M),
     > F1a_b(M,M),F1b_b(M,M),FON_a(M,M),FON_b(M,M)
       integer*8,intent(in) :: devPtrX
       REAL*8,intent(in) :: g,Fxx,Fyy,Fzz
       REAL*8, intent(in) :: factorial(NBCH)
       REAL*8,allocatable :: F3(:,:),FBA(:,:)
       integer :: i,j,k,kk,stat
       real*8 :: E2, tdstep1
#ifdef TD_SIMPLE
       COMPLEX*8, intent(in) :: rho2_a(M,M),rho2_b(M,M)
       COMPLEX*8,allocatable :: rho4(:,:),rho2t(:,:)
#else
       COMPLEX*16, intent(in) :: rho2_a(M,M),rho2_b(M,M)
       COMPLEX*16,allocatable :: rho4(:,:),rho2t(:,:)
#endif
       integer*8 devPtrScratch1
!------------------------------------------------------------------------------!
       ALLOCATE(rho4(M,M),rho2t(M,M),F3(M,M),FBA(M,M))
c
       M2=2*M
       MM=M*(M+1)/2
c first i
       M1=1
c now Fold
       M3=M1+MM
c now S, F also uses the same position after S was used
       M5=M3+MM
c now G
! tdstep predictor is 0.5*tdstep magnus
       tdstep1=tdstep*0.50D0
!
       F3=(7.D0/4.D0)*F1b_a-(3.D0/4.D0)*F1a_a
!
!       call magnus(F3,rho2_a,rho2t,M,NBCH,tdstep1,factorial)
       call cumagnusfac(F3,rho2_a,rho2t,M,NBCH,tdstep1,factorial)
!
!       call matmulnanoc(rho2t,xtrans,rho4,M)
!       call rho_transform(rho2t,devPtrX,rho4,M)
       call complex_rho_on_to_ao(rho2t,devPtrX,rho4,M)
       call sprepack_ctr('L',M,rhoalpha,rho4)
       call sprepack_ctr('L',M,RMM,rho4)
!
       F3=(7.D0/4.D0)*F1b_b-(3.D0/4.D0)*F1a_b
!
!       call magnus(F3,rho2_b,rho2t,M,NBCH,tdstep1,factorial)
        call cumagnusfac(F3,rho2_b,rho2t,M,NBCH,tdstep1,factorial)
!
       rho4=0
!       call matmulnanoc(rho2t,xtrans,rho4,M)
!       call rho_transform(rho2t,devPtrX,rho4,M)
       call complex_rho_on_to_ao(rho2t,devPtrX,rho4,M)
       call sprepack_ctr('L',M,rhobeta,rho4)
       DO i=1,MM
          RMM(i)=RMM(i)+rhobeta(i)
       ENDDO
! Step4: Density matrix 4 is used to calculate F5
       call int3lu(E2)
       call g2g_solve_groups(0,Ex,0)
      if (field) then
         write(*,*) 'FIELD PREDICTOR'
         call dip2(g,Fxx,Fyy,Fzz)
      endif
!
       call spunpack('L',M,RMM(M5),FBA)
!       call matmulnano(FBA,X,FON_a,M)
       call fock_ao_to_on(FBA,devPtrX,FON_a,M)
       call spunpack('L',M,RMM(M3),FBA)
!       call matmulnano(FBA,X,FON_b,M)
       call fock_ao_to_on(FBA,devPtrX,FON_b,M)
       DEALLOCATE(rho4,rho2t,F3,FBA)
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

