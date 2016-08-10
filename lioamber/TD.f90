!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE TD()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! REAL TIME-TDDFT
! 
! Dario Estrin, 1992
! Nano, Dario, Uriel, Damian 2012
!
!  This subrutine takes the conve.GE. density matrix from an SCF calculation
!  and evolves it in time. In the input FILE the total number of propagation
!  steps is specified (nstep) as well as the time of each evolution step 
!  (tdstep). 
!  This implementation has two alternatives to evolve the density in time. The 
!  first one (propaga.OR.1) is the Verlet al.OR.thm that uses a convination of 
!  Liouville von Newmann expresion .OR.the time derivative of the density matrix 
!  and a first.OR.er Tay.OR.expansion of the density matrix. The second one 
!  (propaga.OR.2) is the Magnus propagation scheme that uses Backer Campbell
!  Haus.OR.f (BCH) .OR.ula. .OR.this reason when Magnus is used the number of 
!  total conmuta.OR. in the BCH espansion has to be specified (NBCH, default=10). 
!  A narrow gaussian type electric field can be introduced during the time 
!  evolution in.OR.er to excite all electronic f.EQ.encies with the same intensity.
!  Once this perturbation is turned on (Field=t, exter=t) each component of the
!  EXTERNAL electric field has to be specified in the input FILE (Fx,Fy,Fz).
!  In each step of the propagation the cartesian components of the sistems dipole
!  are s.OR.d in FILEs x.dip, y.dip, z.dip.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       USE garcha_mod, ONLY : RMM,propagator,NBCH,tdstep,idip,tdrestart,&
       exists,M,NCO,natom,Nang,Nunp,Iz,natomc,d,r,atmin,rmax,jatc,&
       nshell,nuc, nnps, nnpp, nnpd, igrid2, predcoef, npas, nsol, pc &
       ,X, Smat, Md, MEMO, ntdstep,tdstep, field, exter, Fx, Fy, Fz, a0 &
       ,WRITEdens, sol,kkind,kkinds,cool,cools, GRAD, epsilon, Iz

       USE ECP_mod, ONLY : ecpmode, term1e, VAAA, VAAB, VBAC,IzECP
       USE mathsubs

#ifdef CUBLAS
       USE cublasmath
#endif
       IMPLICIT NONE
!!!!!!!!!!!!! agregadas .OR.nick
       INTEGER :: ipop,Ndens, MM, MM2, MMd, Md2, M2, niter, Nel, igpu, initial_step
       INTEGER :: M1,M3,M5,M7,M9,M11,M13,M15,M17,M18,M19,M20, info
       REAL*8 :: E,E1,En,Ex,Es, sq2,Qc,Qc2,E1s,Ens, ux, uy, uz, g,factor &
       fxx,fyy,fzz,Enick
       INTEGER :: i,j,k,kk,n !auxiliares
       REAL*8 :: zij, ti, tj, alf,rexp,t0,ff !auxiliares
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

       INTEGER :: istep
       REAL*8 :: t,E2
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: xnano2,xmm,xtrans,ytrans,Y,&
       fock,F1a,F1b,overlap,rhoscratch
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: elmu
#ifdef TD_SIMPLE
       COMPLEX*8 :: Im,Ix
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) :: rho,rhonew,rhold,xnano,&
       rho1
#else
       COMPLEX*16 :: Im,Ix
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: rho,rhonew,rhold,xnano,&
       rho1
#endif
       INTEGER, DIMENSION (natom) :: q
       REAL*8,DIMENSION(:),ALLOCATABLE :: factorial
       INTEGER            :: LWORK,ii,jj
       REAL*8,ALLOCATABLE :: WORK(:)
!!------------------------------------!!
!! FFR ADD
       INTEGER :: pert_steps,lpfrg_steps,chkpntF1a,chkpntF1b
       REAL*8 :: dt_magnus,dt_lpfrg
       LOGICAL :: just_int3n,ematalloct
!! CUBLAS
#ifdef CUBLAS
       INTEGER :: sizeof_REAL
       PARAMETER(sizeof_REAL=8)
       INTEGER :: sizeof_complex
#ifdef TD_SIMPLE
       PARAMETER(sizeof_complex=8)
#else
       PARAMETER(sizeof_complex=16)
#endif
       INTEGER :: stat
       INTEGER*8 :: devPtrX, devPtrY,devPtrXc
       EXTERNAL :: CUBLAS_INIT, CUBLAS_SET_MATRIX
       EXTERNAL :: CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_GET_MATRIX
       EXTERNAL :: CUBLAS_FREE
       INTEGER :: CUBLAS_ALLOC, CUBLAS_SET_MATRIX,CUBLAS_GET_MATRIX
#endif
!!   GROUP OF CHARGES
       LOGICAL             :: groupcharge
       INTEGER             :: ngroup
       INTEGER,ALLOCATABLE :: group(:)
       REAL*8,ALLOCATABLE  :: qgr(:)
       REAL*8 :: tiempo1000
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       CALL g2g_timer_start('TD')
       CALL g2g_timer_start('inicio')
       just_int3n = .false.
       initial_step=0 !para restart
       ALLOCATE(factorial(NBCH))
!!------------------------------------!!
! Mulliken
       ipop=1
! Group of cha.GE.
       groupcha.GE..false.
!!------------------------------------!!
#ifdef CUBLAS
       WRITE(*,*) 'USING CUBLAS'
       CALL CUBLAS_INIT()
       IF (stat.NE.0) THEN
         WRITE(*,*) "initialization failed -TD"
         CALL CUBLAS_SHUTDOWN
         STOP
       ENDIF
#endif
#ifdef TD_SIMPLE
       WRITE(*,*) 'simple presition complex'
#else
       WRITE(*,*) 'double presition complex'
#endif
       IF(propaga.OR.EQ.2) THEN
         dt_magnus=tdstep
         dt_lpfrg=tdstep*0.10D0
         factorial(1)=1.0D0
#ifdef CUBLAS
         DO ii=1,NBCH
           factorial(ii)=1.0D0/ii
         ENDDO
#else     
         DO ii=2,NBCH
           factorial(ii)=factorial(ii-1)/ii
          ENDDO
#endif
       ENDIF

       IF(propaga.OR.EQ.1) dt_lpfrg=tdstep

!!------------------------------------!!
!! FFR ADD:
       pert_steps=100
       lpfrg_steps=200
       chkpntF1a=185
       chkpntF1b=195
!--------------------------------------------------------------------!
! Pointers -
       Ndens=1
       E=0.0D0
       E1=0.0D0
       En=0.0D0
       E2=0.0D0
       Ens=0.d0 !agregada .OR.nick
       Ex=0.d0 !agregada .OR.nick
       idip=1
       Im=(0.0D0,2.0D0)
       sq2=sqrt(2.D0)
       MM=M*(M+1)/2 
       MM2=M**2
       MMd=Md*(Md+1)/2
       Md2=2*Md
       M2=2*M
!
       ALLOCATE(xnano(M,M),xnano2(M,M),fock(M,M),rhonew(M,M), &
       rhold(M,M),rho(M,M),xmm(M,M),xtrans(M,M),Y(M,M),ytrans(M,M),&
       rho1(M,M))
!
       IF(propaga.OR.EQ.2) ALLOCATE (F1a(M,M),F1b(M,M))


!--------------------------------------------------------------------!
!%%%%%%%%%%%%%%%%%%% RESTART CASE %%%%%%%%%%%%%%%%%%%%%%%!
       IF (tdrestart) THEN  !pasar a un modulo y en FORMATo binario, Nick
! We READ the density matrix s.OR.d in RMM(1,2,3,...,MM) and it is copied in rho matrix.
         CALL READ_TD_RESTART(rho,RMM,F1a,F1b,M2,initial_step,t)
       ELSE
         CALL spunpack_rtc('L',M,RMM,rho)!agarra rho de rmm salida del SCF
       ENDIF
!%%%%%%%%%%%%%%%%%%% RESTART CASE %%%%%%%%%%%%%%%%%%%%%%%!

!------------------------------------------------------------------------------!

!!%%%%%%%%%%%%%%%%%%% POINTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%!
       M1=1 !c first i
       M3=M1+MM !c now Fold
       M5=M3+MM !c now S, F also uses the same position after S was used
       M7=M5+MM !c now G
       M9=M7+MMd !c now Gm
       M11=M9+MMd !c now H
       M13=M11+MM !c W ( e.GE.values ), also this space is used in least squares
       M15=M13+M !c aux ( vec.OR..OR.ESSl)
       M17=M15+MM !c Least squares
       M18=M17+MMd !c vec.OR. of MO
       M19=M18+M*NCO !c weights (in case of using option )
       M20 = M19 + natom*50*Nang   !c RAM s.OR.GE.of two-electron integrals (IF MEMO=T)
       Nel=2*NCO+Nunp

!c Initializations/Defaults
       WRITE(*,*) ' Starting TD CALCULATION  '
!--------------------------------------!
       niter=0
!--------------------------------------!
       Qc=0.0D0
       DO i=1,natom
         Qc=Qc+Iz(i)
       ENDDO
         Qc=Qc-Nel
         Qc2=Qc**2
!------------------------------------------------------------------------------!
! Two electron integral with neigh.OR.list., Ver si esto se usa -pasar a un modulo (y agregarlo a scf)
!
       DO i=1,natom
         natomc(i)=0
!
         DO j=1,natom
           d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+&
           (r(i,3)-r(j,3))**2

           zij=atmin(i)+atmin(j)
           ti=atmin(i)/zij
           tj=atmin(j)/zij
           alf=atmin(i)*tj
           rexp=alf*d(i,j)
           IF (rexp.LT.rmax) THEN
             natomc(i)=natomc(i)+1
             jatc(natomc(i),i)=j
           ENDIF 
         ENDDO
       ENDDO

       DO ii=nshell(0),1,-1
         nnps(nuc(ii))=ii
       ENDDO
       DO ii=nshell(0)+nshell(1),nshell(0)+1,-1
         nnpp(nuc(ii))=ii
       ENDDO
       DO ii=M,nshell(0)+nshell(1)+1,-1
         nnpd(nuc(ii))=ii
       ENDDO
!------------------------------------------------------------------!
!c
!c Create integration grid .OR.XC here
!c Assign points to groups (spheres/cubes)
!c Assign significant functions to groups
!c -Calculate point weights
!c
       CALL g2g_timer_sum_start('Excha.GE..OR.elation grid setup')
       CALL g2g_reload_atom_positions(igrid2)
       CALL g2g_timer_sum_STOP('Excha.GE..OR.elation grid setup')

       CALL aint_query_gpu_level(igpu)
       IF (igpu.GT.1) CALL aint_new_step()

!------------------------------------------------------------------------------!
! H H .OR., 1 electron matrix elements
       CALL g2g_timer_sum_start('1-e Fock')
       CALL g2g_timer_sum_start('Nuclear attraction')
       CALL int1(En)

       IF (ecpmode) THEN
         DO k=1,MM
           term1e(k)=RMM(M11+k-1) !copia los terminos de 1e-
           RMM(M11+k-1)=RMM(M11+k-1)+VAAA(k)+VAAB(k)+VBAC(k) !agrega el ECP a los terminos de 1 e
         ENDDO
       ENDIF

       CALL g2g_timer_sum_STOP('Nuclear attraction')
       IF(nsol.GT.0.OR.igpu.GE.4) THEN
         CALL g2g_timer_sum_start('QM/MM')
         IF (igpu.LE.1) THEN
           CALL g2g_timer_start('intsol')
           CALL intsol(E1s,Ens,.true.)
           CALL g2g_timer_STOP('intsol')
         ELSE
           CALL aint_qmmm_init(nsol,r,pc)
           CALL g2g_timer_start('aint_qmmm_fock')
           CALL aint_qmmm_fock(E1s,Ens)
           CALL g2g_timer_STOP('aint_qmmm_fock')
         ENDIF
           CALL g2g_timer_sum_STOP('QM/MM')
       ENDIF
!--------------------------------------!

       E1=0.D0
       DO k=1,MM
         E1=E1+RMM(k)*RMM(M11+k-1)
       ENDDO
       CALL g2g_timer_sum_STOP('1-e Fock')
       IF(ipop.EQ.1)THEN
         ALLOCATE(overlap(M,M),rhoscratch(M,M))
         CALL spunpack('L',M,RMM(M5),overlap) !copy overlap matrix calculated in SCF from RMM 
       ENDIF

!--------------------------------------!
!c Diagonalization of S matrix, after this is.NOT.needed any.OR.
!c s is in RMM(M13,M13+1,M13+2,...,M13+MM)
!--------------------------------------!
! ESSL OPTION
#ifdef essl
       CALL DSPEV(1,RMM(M5),RMM(M13),X,M,M,RMM(M15),M2) !diagonaliza overlap en RMM5, en RMM13 pone autovect y en RMM15 autova.OR.s
#endif
!--------------------------------------!
! LAPACK OPTION
#ifdef pack
       DO ii=1,M; DO jj=1,M
         X(ii,jj)=Smat(ii,jj)
       ENDDO; ENDDO

       IF (ALLOCATEd(WORK)) DEALLOCATE(WORK); ALLOCATE(WORK(1))
       CALL dsyev('V','L',M,X,M,RMM(M13),WORK,-1,info)
       LWORK=int(WORK(1));  DEALLOCATE(WORK); ALLOCATE(WORK(LWORK))
       CALL dsyev('V','L',M,X,M,RMM(M13),WORK,LWORK,info)
#endif
!--------------------------------------!


! Here, we obtain the transFORMATion matrices X and Y .OR.converting 
! from the atomic.OR.ital to a molecular.OR.ital basis (truncated
! during linear dependency elimination). 
! S is the overlap matrix
! s is the diagonal e.GE.value matrix of S
! U is the e.GE.vec.OR.matrix of S
! X=U s^(-1/2)
! matrix X's DIMENSION is M*3M. In the first M*M terms it contains
! the transFORMATion matrices and in the other M*2M terms it contains
! auxiliar matrices.
       CALL g2g_timer_start('inicio1')
         DO j=1,M
           IF (RMM(M13+j-1).LT.1.0D-06) THEN
             WRITE(*,*) 'LINEAR DEPENDENCY DETECTED'
             DO i=1,M
               X(i,j)=0.0D0
               Y(i,j)=0.0D0
             ENDDO
           ELSE
             DO i=1,M
               Y(i,j)=X(i,j)*sqrt(RMM(M13+j-1))
               X(i,j)=X(i,j)/sqrt(RMM(M13+j-1))       
             ENDDO
           ENDIF
         ENDDO
!------------------------------------------------------------------------------!
#ifdef CUBLAS
         DO i=1,M
           DO j=1,M
             rho1(i,j)=cmplx(X(i,j),0.0D0)
           ENDDO
         ENDDO
         stat = CUBLAS_ALLOC(M*M, sizeof_REAL, devPtrX)
         stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
         stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
         
         IF (stat.NE.0) THEN
           WRITE(*,*) "X and.OR.Y me.OR. allocation failed"
           CALL CUBLAS_SHUTDOWN
           STOP
         ENDIF
         stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrXc,M)
         stat=CUBLAS_SET_MATRIX(M,M,sizeof_REAL,x,M,devPtrX,M)
         DO i=1,M
           DO j=1,M
             rho1(i,j)=cmplx(Y(i,j),0.0D0)
           ENDDO
         ENDDO
         stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrY,M)
         IF (stat.NE.0) THEN
           WRITE(*,*) "X and.OR.Y setting failed"
           CALL CUBLAS_SHUTDOWN
           STOP
         ENDIF
         rho1=0
#endif
!------------------------------------------------------------------------------!
! the transFORMATion matrices is copied in xmm
!
         DO i=1,M
           DO j=1,M
             xmm(i,j)=X(i,j)
           ENDDO
         ENDDO
! the tranposed matrixes are calculated
         DO i=1,M
           DO j=1,M
             xtrans(j,i)=X(i,j)
             ytrans(j,i)=Y(i,j)
           ENDDO
         ENDDO
!------------------------------------------------------------------------------!
! External Electric Field components
!
!       WRITE(*,*) 'fx =', fx
!       WRITE(*,*) 'fy =', fy
!       WRITE(*,*) 'fz =', fz
!------------------------------------------------------------------------------!
! Rho is trans.OR.ed to the.OR.ho.OR.al basis
! with matmul:
#ifdef CUBLAS

         CALL g2g_timer_start('complex_rho_ao_to_on-cu')
         rho1=basecha.GE.cublas(M,rho,devPtrY,'dir')
         rho=rho1
         CALL g2g_timer_STOP('complex_rho_ao_to_on-cu')
#else
         rho=matmul(ytrans,rho)
         rho=matmul(rho,Y)
#endif
! with matmulnanoc
! (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
!            CALL matmulnanoc(rho,Y,rho1,M)
!            rho=rho1
!            rho=basecha.GE.M,Ytrans,rho,Y)
!--------------------------------------!
!c Precalculate three-index (two in MO basis, one in density basis) matrix
!c used in density fitting / Coulomb F element calculation here
!c (t_i in Dunlap)
!c
         CALL aint_query_gpu_level(igpu)
         IF (igpu.GT.2) THEN
           CALL aint_coulomb_init()
         ENDIF


         IF (igpu.EQ.5) MEMO = .false.
         !MEMO=.true.
         IF (MEMO) THEN
           CALL g2g_timer_start('int3mem')
           CALL g2g_timer_sum_start('Coulomb precalc')
!c La.GE.elements of t_i put into double-precision cool here
!c Size criteria based on size of pre-fac.OR.in Gaussian Product Th.OR.m
!c (applied to MO basis indices)
           CALL int3mem()
!c Small elements of t_i put into single-precision cools here
!c         CALL int3mems()
           CALL g2g_timer_STOP('int3mem')
           CALL g2g_timer_sum_STOP('Coulomb precalc')
         ENDIF
#ifdef CUBLAS
         CALL CUBLAS_FREE(devPtrY)
#endif

         CALL g2g_timer_STOP('inicio')
!##############################################################################!
! HERE STARTS THE TIME EVOLUTION
!##############################################################################!

         WRITE(*,*) 'PROPAGATION'
         DO 999 istep=1, ntdstep
!--------------------------------------!
           CALL g2g_timer_start('TD step')
           CALL CALC_CURRENT_TIME(propaga.OR.lpfrg_steps, initial_step &
           ,istep, tdstep,t) !calcula t del paso actual
!--------------------------------------!
           CALL int3lu(E2) !calvula fock
           CALL g2g_solve_groups(0,Ex,0) !otra parte de fock, intercambio y .OR.elacion (no calcula energia)
!           CALL g2g_solve_groups(1,Ex,0) !calcula energia, creo q es mas caro en el calculo, nick

           E=E1+E2+En+Ex
           WRITE(*,*) '! step & energy -Exc', istep+initial_step-1,E !despues poner step -1  a.OR. puse +500 para comparar facil
!              WRITE(*,*) "!# E1,E2,En,Ex",E1,E2,En,Ex

           E1=0.0D0

!c ELECTRIC FIELD CASE - Type=gaussian (ON)  !rutina field o fld !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(istep+initial_step.LT.pert_steps) THEN
             CALL FIELD_INT(istep+ initial_step,Qc2,E1) !calcula los terminos de fock y la energia debidos al campo electrico aplicado
!               WRITE(*,*) "campo prendido!!!!!!!!!!!!!!!!!"
           ELSE
             field=.false.
           ENDIF
!------------------------------------------------------------------------------!
! E1 includes solvent 1 electron contributions
           DO k=1,MM
             E1=E1+RMM(k)*RMM(M11+k-1) !agrega otros terminos de 1e
           ENDDO
!        WRITE(*,*) '1 electron contribution',E1
!------------------------------------------------------------------------------!
! Here we obtain the fock matrix in the molecular.OR.ital (MO) basis.
! where U matrix with e.GE.vec.OR. of S , and s is vec.OR.with
! e.GE.values
           CALL g2g_timer_start('fock')
           CALL spunpack('L',M,RMM(M5),fock) !secopia fock de rmm
#ifdef CUBLAS
           xnano2 = basecha.GE.cublas(M,fock,devPtrX,'dir') !pasa fock a la base.OR.o.OR.al
           fock=xnano2

#else
           xnano2=matmul(xtrans,fock) !xtrans xtranx trnspuesta de X (cambio de base)
           fock=matmul(xnano2,xmm) !termian el cambio de base
!             CALL fock_ao_to_on(fock,x,fock,M)
#endif

           CALL sprepack('L',M,RMM(M5),fock) !pasa fock a vec.OR.en rmm
           CALL g2g_timer_STOP('fock')
!c Now fock is s.OR.d in molecular.OR.ital basis.
!c
!  s.OR.s F1a and F1b .OR.magnus propagation
           IF((propaga.OR.EQ.2) .AND. (.NOT.tdrestart)) THEN
             IF(istep.EQ.chkpntF1a) F1a=fock         
             IF(istep.EQ.chkpntF1b) F1b=fock         
           ENDIF
!--------------------------------------------------------------------!
           IF ((propaga.OR.EQ.1).OR. &
           (((propaga.OR.EQ.2).and.(istep.LT.lpfrg_steps)) &
           .AND. (.NOT.tdrestart))) THEN
             WRITE(*,*) 'Verlet'
!c In the first step of the propagation we extrapolate rho back in time
!c using Verlet al.OR.thm to calculate rhold.
!c using matmul 
!           IF(istep.EQ.1) THEN
!             rhold=rho+(dt_lpfrg*Im*(matmul(fock,rho)))
!             rhold=rhold-(dt_lpfrg*Im*(matmul(rho,fock)))
!           ENDIF
!c using commutator
             IF (istep.EQ.1) THEN
#ifdef CUBLAS
               CALL g2g_timer_start('cuconmut')
               rhold=commuta.OR.cublas(fock,rho)
               rhold=rho+dt_lpfrg*(Im*rhold)
               CALL g2g_timer_STOP('cuconmut')
#else
               CALL g2g_timer_start('conmutc')
               rhold=commuta.OR.fock,rho)
               rhold=rho+dt_lpfrg*(Im*rhold)
               CALL g2g_timer_STOP('conmutc')
#endif
             ENDIF
!####################################################################!
! DENSITY MATRIX PROPAGATION USING VERLET ALGORITHM
! using matmul:
!           rhonew=rhold-(dt_lpfrg*Im*(matmul(fock,rho)))
!           rhonew=rhonew+(dt_lpfrg*Im*(matmul(rho,fock)))
!c--------------------------------------c
! using commuta.OR.
             CALL g2g_timer_start('commutator')
#ifdef CUBLAS
             rhonew=commuta.OR.cublas(fock,rho) !dp/dt
             rhonew=rhold-dt_lpfrg*(Im*rhonew) !p(t+dt)
#else
             rhonew=commuta.OR.fock,rho)
             rhonew=rhold-dt_lpfrg*(Im*rhonew)
#endif
             CALL  g2g_timer_STOP('commutator')
!c Density update (rhold-->rho, rho-->rhonew)
             DO i=1,M
               DO j=1,M
                 rhold(i,j)=rho(i,j)
                 rho(i,j)=rhonew(i,j)
               ENDDO
             ENDDO
! END OF VERLET PROPAGATOR
!####################################################################!
           ELSE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! DENSITY MATRIX PROPAGATION USING MAGNUS ALGORITHM
             WRITE(*,*) 'Magnus'
#ifdef CUBLAS
             CALL g2g_timer_start('cupredictor')
             CALL cupredictor(F1a,F1b,fock,rho,devPtrX,factorial, &
      fxx,fyy,fzz,g,devPtrXc) 
             WRITE(*,*) 'Magnus1'
             CALL g2g_timer_STOP('cupredictor')
             CALL g2g_timer_start('cumagnus')
             CALL cumagnusfac(fock,rho,rhonew,M,NBCH,dt_magnus, &
      factorial)
             WRITE(*,*) 'Magnus2'
             CALL g2g_timer_STOP('cumagnus')
!                rhold=rhonew
!                CALL g2g_timer_start('MAGNUS_MODIFIED')
!                CALL magnus_cublas(fock,rho,rhonew,M,NBCH,dt_magnus,
!     >factorial) 
!                CALL g2g_timer_STOP('MAGNUS_MODIFIED')
!                rhold=rhonew-rhold
!                WRITE(22222222,*) rhold
!                STOP 'hemos escrito rhold'
#else
             CALL g2g_timer_start('predictor')
             CALL predictor(F1a,F1b,fock,rho,factorial, &
       fxx,fyy,fzz,g)
             WRITE(*,*) 'Magnus3'
             CALL g2g_timer_STOP('predictor')
             CALL g2g_timer_start('magnus')
             CALL magnus(fock,rho,rhonew,M,NBCH,dt_magnus,factorial)
             WRITE(*,*) 'Magnus4'
             CALL g2g_timer_STOP('magnus')
#endif
             F1a=F1b
             F1b=fock
             rho=rhonew
! END OF MAGNUS PROPAGATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
           ENDIF
!####################################################################!
!c Here we trans.OR. the density to the atomic.OR.ital basis and take the REAL part of it. The imaginary part of the density 
!c can be descarted since .OR.a basis set of purely REAL functions the fock matrix is REAL and symetric and depends only on 
!c the REAL part of the complex density matrix. (This wont be true in the case of hybrid functionals)
!c with matmul:
#ifdef CUBLAS
           CALL g2g_timer_start('complex_rho_on_to_ao-cu')
           rho1=basecha.GE.cublas(M,rho,devPtrXc,'inv') !pasa rho a OA
           CALL g2g_timer_STOP('complex_rho_on_to_ao-cu')
#else
           CALL g2g_timer_start('complex_rho_on_to_ao')
           rho1=matmul(x,rho)
           rho1=matmul(rho1,xtrans)
           CALL g2g_timer_STOP('complex_rho_on_to_ao')
#endif
!       rho1=REAL(rho1)
!c with matmulnanoc:
!c (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
!c          CALL matmulnanoc(rho,xtrans,rho1,M)
!c          rho1=basecha.GE.M,X,rho,Xtrans)
!c          rho1 = REAL(rho1)
!c The REAL part of the density matrix in the atomic.OR.ital basis is copied in RMM(1,2,3,...,MM) to compute the .OR.esponding fock matrix.
           DO j=1,M
             DO k=j,M
               IF(j.EQ.k) THEN
                 RMM(k+(M2-j)*(j-1)/2)=REAL(rho1(j,k))
               ELSE
                 RMM(k+(M2-j)*(j-1)/2)=(REAL(rho1(j,k)))*2
               ENDIF
             ENDDO
           ENDDO

! S.OR.s the density matrix each 500 steps to be able to restart the dynamics, and magnus focks
           IF (WRITEdens) CALL WRITE_RESTART(istep+ initial_step,t,rho1,f1a,f1b) !S.OR.s the density matrix each 500 steps to be able to restart the dynamics

!###################################################################!
!# DIPOLE MOMENT CALCULATION
           CALL WRITE_DIPOLE_MOMENT(istep,t, lpfrg_steps)
!c u in Debyes
!# END OF DIPOLE MOMENT CALCULATION
!c------------------------------------------------------------------------------------
!c      WRITE(*,*) 'Coulomb E',E2-Ex,Ex
           CALL g2g_timer_STOP('TD step')
           WRITE(*,*)
 999     CONTINUE
!
!##############################################################################!
! HERE FINISHES THE PROPAGATION
!##############################################################################!

 995     CONTINUE
!c
!c
         IF (memo) THEN
           DEALLOCATE (kkind,kkinds)
           DEALLOCATE(cool,cools)
         ENDIF

         IF(propaga.OR.EQ.2) DEALLOCATE (F1a,F1b)

!c calculation of energy weighted density matrix
!c
         kk=0
         DO 307 j=1,M
           DO 307 i=j,M
             kk=kk+1
             RMM(M15+kk-1)=0.D0
             IF(i.EQ.j) THEN
               ff=2.D0
             ELSE
               ff=4.D0
             ENDIF
             DO 309 k=1,NCO
               RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1) &
       *ff*X(i,M2+k)*X(j,M2+k)
 309  CONTINUE
 307   CONTINUE
!c

!despues ver este nopt, nick
!          IF (nopt.EQ.0) THEN
!c calculates Mulliken poputations


       CALL g2g_timer_sum_start('Mulliken')


! MULLIKEN POPULATION ANALYSIS (FFR - Simplified)
!--------------------------------------------------------------------!
!       IF (ipop.EQ.1) THEN
!       CALL int1(En)
!       CALL spunpack('L',M,RMM(M5),Smat)
!       CALL spunpack('L',M,RMM(M1),Rho)
!       CALL fixrho(M,Rho)
!       CALL mulliken_calc(natom,M,Rho,Smat,Nuc,Iz,q)

!       IF (ecpmode) THEN
!Modification .OR.Effective .OR. Potential, Nick
!          CALL mulliken_WRITE(85,natom,IzECP,q)
!       ELSE
!          CALL mulliken_WRITE(85,natom,Iz,q)
!       ENDIF

! NOTE: If 'mulliken_calc' is renamed as 'mulliken', the code will
! malfunction. I DON'T KNOW WHY.
!       CALL g2g_timer_sum_STOP('Mulliken')
!       ENDIF
!--------------------------------------------------------------------!
 

!         IF (ipop.EQ.1) THEN
!           CALL int1(En)
!           DO n=1,natom
!             q(n)=Iz(n)
!           ENDDO
!             DO i=1,M
!               DO j=1,i-1
!                 kk=i+(M2-j)*(j-1)/2
!                 t0=RMM(kk)*RMM(M5+kk-1)/2.D0
!                 q(Nuc(i))=q(Nuc(i))-t0
!               ENDDO
!               kk=i+(M2-i)*(i-1)/2
!               t0=RMM(kk)*RMM(M5+kk-1)
!               q(Nuc(i))=q(Nuc(i))-t0
!               DO j=i+1,M
!                 kk=j+(M2-i)*(i-1)/2
!                 t0=RMM(kk)*RMM(M5+kk-1)/2.D0
!                 q(Nuc(i))=q(Nuc(i))-t0
!              ENDDO
!             ENDDO
!             WRITE(*,*) 'MULLIKEN POPULATION ANALYSIS'
!             WRITE(*,770)
!             DO n=1,natom
!               WRITE(*,760) n,Iz(n),q(n)
!             ENDDO
!             WRITE(*,*)
!           ENDIF

!c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
!c
!c        IF (icha.GE.EQ.1) THEN
!c          Q1=-(2*NCO+Nunp)
!c         DO n=1,natom
!c          Q1=Q1+Iz(n)
!c         ENDDO
!c         CALL cha.GE.NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
!c     >            c,a,RMM,map,Q1)
!c        ENDIF
!c
!!c outputs final  MO ---------------------
!      DO l=1,M
!c      DO n=1,NCO+3
!      DO n=1,M
!        X(indexii(l),M+n)=X(l,M2+n)
!      ENDDO
!      ENDDO
!
!      DO l=1,M
!        WRITE(2,400) (X(l,M+n),n=1,NCO)
!      ENDDO
!--------------------------------------!
! Writes down MO coefficients and.OR.ital energies
!       WRITE(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
!       DO n=1,NCO
!         WRITE(29,850) n,RMM(M13+n-1)
!         WRITE(29,400) (X(l,M+n),l=1,M)
!       ENDDO
!       DO n=NCO+1,M
!         WRITE(29,851) n,RMM(M13+n-1)
!         WRITE(29,400) (X(l,M+n),l=1,M)
!       ENDDO
!       CLOSE(29)
!--------------------------------------!
!      ENDIF
#ifdef CUBLAS
           CALL CUBLAS_FREE(devPtrX)
           CALL CUBLAS_FREE(devPtrXc)
           CALL CUBLAS_FREE(devPtrY)
           CALL CUBLAS_SHUTDOWN()
#endif
!c
!c
!c---- DEBUGGINGS
!c      WRITE(*,*) 'Exc, integrated and calculated',Exc,Ex
!c      WRITE(*,*) 'Coulomb energy',E2-Ex
!c
           CALL g2g_timer_STOP('TD')
           DEALLOCATE(xnano,fock,rho)
           DEALLOCATE(factorial)
!------------------------------------------------------------------------------!
 500  FORMAT('SCF TIME ',I6,' sec')
 450  FORMAT ('FINAL ENERGY = ',F19.12)
 400  FORMAT(4(E14.7E2,2x))
 300  FORMAT(I3,E14.6,2x,F14.7)
 600  FORMAT('  ENERGY CONTRIBUTIONS IN A.U.')
 610  FORMAT(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  FORMAT(2x,'SOLVENT')
 620  FORMAT(F14.7,4x,F14.7,4x,F14.7)
 625  FORMAT(F14.7)
 760  FORMAT(I3,9x,I3,6x,F10.4)
 770  FORMAT('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
 850  FORMAT('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
 851  FORMAT('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7, &
         '(NON OCC.)')
 900  FORMAT(F15.9,2x,3(F15.9,2x),2x,F15.9)
 901  FORMAT(F15.9,2x,F15.9)
 777  FORMAT(4(F8.4,2x))
 776  FORMAT (3(F8.4,2x))

      RETURN
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

       CONTAINS

      SUBROUTINE CALC_CURRENT_TIME(propaga.OR.lpfrg_steps, initial_step &
      , istep, tdstep,t)
      IMPLICIT NONE
      INTEGER, intent(in) :: propaga.OR.lpfrg_steps, initial_step
      REAL*8, intent(inout) :: t
      INTEGER, intent(in) :: istep
      REAL*8, intent(in) :: tdstep
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
      IF (istep.GT.1) THEN
        IF ((propaga.OR.EQ.2) .AND. (initial_step+istep.LE.lpfrg_steps)) THEN
          t=t+tdstep*0.002419
        ELSE
          t=t+tdstep*0.02419
        ENDIF
        WRITE(*,*) 'evolution time (fs)  =', t
      ENDIF
!------------------------------------------------------------------------------!
       RETURN
      END SUBROUTINE



      SUBROUTINE FIELD_INT(istep,Qc2,E1)
      USE garcha_mod, ONLY : field,exter, Fx, Fy, Fz, epsilon,a0
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: istep
      REAL*8, INTENT(IN) :: Qc2
      REAL*8, INTENT(INOUT) :: E1
      DOUBLE PRECISION :: g, ux, uy, uz, factor,fxx,fyy,fzz
        IF (field) THEN
          CALL dip(ux,uy,uz)
          IF (exter) THEN
            g=1.0D0
            fac.OR.2.54D0
            fxx=fx*exp(-0.2*(REAL(istep-50))**2)
            fyy=fy*exp(-0.2*(REAL(istep-50))**2)
            fzz=fz*exp(-0.2*(REAL(istep-50))**2)
!            WRITE(*,*) fxx,fyy,fzz
          ELSE
            g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
            Fx=ux/2.54D0
            Fy=uy/2.54D0
            Fz=uz/2.54D0
            fac.OR.(2.54D0*2.00D0)
          ENDIF
          CALL intfld(g,Fxx,Fyy,Fzz) !calcula integrales del campo con mom dipolar y da fock
           E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/factor &
           0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0  !E de 1 electron .OR.el campo
        ENDIF
      RETURN
      END SUBROUTINE FIELD_INT


      SUBROUTINE WRITE_DIPOLE_MOMENT(istep,t, lpfrg_steps)
      USE garcha_mod, ONLY : propagator, tdrestart, ntdstep
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: istep, lpfrg_steps
      REAL*8, INTENT(IN) :: t
      DOUBLE PRECISION :: ux,uy,uz

        IF(istep.EQ.1) THEN
          OPEN(UNIT=134,FILE='x.dip')
          OPEN(UNIT=135,FILE='y.dip')
          OPEN(UNIT=136,FILE='z.dip')
          OPEN(UNIT=13600,FILE='abs.dip')

!aca hay q agregar q escriba ts  NCO  field en cada archivo, si o es splito propagation en NCO poner 1
        WRITE(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
        WRITE(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
        WRITE(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
        WRITE(13600,*) '#Time (fs) vs DIPOLE MOMENT (DEBYES)'
        ENDIF

              IF ((propagator.EQ.2).and.(istep.LT.lpfrg_steps) &
             .and. (.NOT.tdrestart)) THEN
                  IF(mod ((istep-1),10) == 0) THEN
                     CALL g2g_timer_start('DIPOLE')
                     CALL dip(ux,uy,uz)
                     CALL g2g_timer_STOP('DIPOLE')
                     WRITE(134,901) t,ux
                     WRITE(135,901) t,uy
                     WRITE(136,901) t,uz
                  ENDIF
              ELSE
                  CALL g2g_timer_start('DIPOLE')
                  CALL dip(ux,uy,uz)
                  CALL g2g_timer_STOP('DIPOLE')
                  WRITE(134,901) t,ux
                  WRITE(135,901) t,uy
                  WRITE(136,901) t,uz
              ENDIF

        IF(istep.EQ.ntdstep) THEN
          CLOSE (134)
          CLOSE(135)
          CLOSE(136)
          CLOSE(13600)
        ENDIF

!c u in Debyes
 901  FORMAT(F15.9,2x,F15.9)
      END SUBROUTINE WRITE_DIPOLE_MOMENT



      SUBROUTINE WRITE_RESTART(istep,t,rho1,f1a,f1b)
      USE garcha_mod, ONLY : WRITEdens, M, ntdstep
      IMPLICIT NONE
      INTEGER :: j,k
      INTEGER, INTENT(IN) :: istep
      REAL*8, INTENT(IN) :: t
      REAL*8, INTENT(IN), DIMENSION(M,M) :: f1a,f1b
#ifdef TD_SIMPLE
       COMPLEX*8,INTENT(IN), DIMENSION(M,M) :: rho1
#ELSE
       COMPLEX*16,INTENT(IN),DIMENSION(M,M) :: rho1
#ENDIF
        IF(mod (istep,500) == 0) THEN
          OPEN(UNIT=5374,FILE="rho.restart",STATUS='UNKNOWN', &
          ACCESS='STREAM')
          REWIND 5374
          WRITE (5374) istep, t
          WRITE(5374) rho1
          CLOSE(5374)
          IF(propagator.EQ.2) CALL WRITE_MAGNUS_RESTART(F1a,F1b) !WRITE_MAGNUS_RESTART(istep, fock)

        ENDIF
! In the last step density matrix is s.OR.d
        IF (istep.EQ.ntdstep) THEN
          OPEN(UNIT=235,FILE="rholast",STATUS='UNKNOWN', &
          ACCESS='STREAM')
          REWIND 235
          WRITE (235) istep, t
          WRITE(235) rho1
          CLOSE(235)
         ENDIF
       END SUBROUTINE WRITE_RESTART

      SUBROUTINE WRITE_MAGNUS_RESTART(F1a,F1b)!istep, fock)
!  s.OR.s F1a and F1b checkpoints to restart the dynamics
      USE garcha_mod, ONLY : propagator,M
      IMPLICIT NONE
      REAL*8, INTENT(IN), DIMENSION(M,M) :: F1a,F1b
      INTEGER :: kk,ii
                 OPEN(UNIT=7624,FILE="F1b.restart",STATUS='UNKNOWN', &
                 ACCESS='STREAM')
                 REWIND 7624
                 WRITE(7624) fock
                 CLOSE(7624)

                 OPEN(UNIT=7625,FILE="F1a.restart",STATUS='UNKNOWN', &
                 ACCESS='STREAM')
                 REWIND 7625
                 WRITE(7625) fock
                 CLOSE(7625)
      END SUBROUTINE WRITE_MAGNUS_RESTART


      SUBROUTINE READ_TD_RESTART(rho,RMM,F1a,F1b,M2,initial_step,t)
      USE garcha_mod, ONLY : propagator, M
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M2
      INTEGER, INTENT(INOUT) :: initial_step
      REAL*8, INTENT(INOUT), DIMENSION(M,M) :: F1a,F1b
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*) :: RMM
      INTEGER :: i,j,k
      LOGICAL :: exists
      REAL*8, intent(inout) :: t
#ifdef TD_SIMPLE
      COMPLEX*8,INTENT(INOUT), DIMENSION(M,M) :: rho
#else
      COMPLEX*16,INTENT(INOUT), DIMENSION(M,M) :: rho
#endif

         INQUIRE(FILE='rho.restart',EXIST=exists)
         IF (.NOT.exists) THEN
             WRITE(*,*) 'ERROR CANNOT FIND rho.restart'
             WRITE(*,*) "(IF you are.NOT.restarting a previous run",&
            " set tdrestart= false)"
             STOP
         ENDIF

	OPEN(UNIT=1544,FILE='rho.restart',STATUS='UNKNOWN', ACCESS='STREAM')
        READ(1544) initial_step, t
        READ(1544) rho
        CLOSE(1544)

         DO j=1,M
            DO k=j,M
               IF(j.EQ.k) THEN
                  RMM(k+(M2-j)*(j-1)/2)=REAL(rho(j,k))
               ELSE
                  RMM(k+(M2-j)*(j-1)/2)=(REAL(rho(j,k)))*2
               ENDIF
            ENDDO
         ENDDO


         IF (propaga.OR..EQ. 2) THEN
            INQUIRE(FILE='F1a.restart',EXIST=exists)
            IF (.NOT.exists) THEN
               WRITE(*,*) 'ERROR CANNOT FIND F1a.restart'
               WRITE(*,*) "(IF you are.NOT.restarting a previous run ",&
               "set tdrestart= false)"
               STOP
            ENDIF
            INQUIRE(FILE='F1b.restart',EXIST=exists)
            IF (.NOT.exists) THEN
               WRITE(*,*) 'ERROR CANNOT FIND F1b.restart'
               WRITE(*,*) "(IF you are.NOT.restarting a previous run ",&
               "set tdrestart= false)"
               STOP
            ENDIF

            OPEN(UNIT=7777,FILE="F1b.restart",STATUS='UNKNOWN', &
            ACCESS='STREAM')
            READ(7777) F1a
            CLOSE(7777)

            OPEN(UNIT=7399,FILE="F1a.restart",STATUS='UNKNOWN', &
            ACCESS='STREAM')
            READ(7399) F1b!(i,j)
         ENDIF
        RETURN
      END SUBROUTINE READ_TD_RESTART

       END SUBROUTINE TD
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
