!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine SCF(E)
!
!------------------------------------------------------------------------------!
! DIRECT VERSION
! Calls all integrals generator subroutines : 1 el integrals,
! 2 el integrals, exchange fitting , so it gets S matrix, F matrix
! and P matrix in lower storage mode (symmetric matrices)
!
! Dario Estrin, 1992
!
!------------------------------------------------------------------------------!
! Modified to f90
!
! Nick, 2017
!
!------------------------------------------------------------------------------!
! * Modified the format to comply to lio current conventions.
!   (TODO: Almost complete, details pending)
! * Re-ordered some sections to improve readability.
! * Externalization of most procedures to dedicated subroutines.
!   (TODO: This is still in progress)
!
! FFR, 03/2018
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use ehrensubs , only: ehrendyn_init
   use garcha_mod, only : M,Md, NCO,natom,Nang, number_restr, hybrid_converg,  &
                          MEMO, npas, verbose, RMM, X, SHFT, GRAD, npasw,      &
                          igrid, energy_freq, converge, noconverge, lowdin,    &
                          cubegen_only, VCINP, primera, Nunp, GOLD, igrid2,    &
                          predcoef, nsol, r, pc, DIIS, told, Etold, Enucl,     &
                          Eorbs, kkind,kkinds,cool,cools,NMAX,Dbug, idip, Iz,  &
                          nuc, doing_ehrenfest, first_step, RealRho,           &
                          total_time, MO_coef_at, MO_coef_at_b, Smat, good_cut,&
                          ndiis, ncont, nshell, rhoalpha, rhobeta, OPEN
   use ECP_mod, only : ecpmode
   use mask_ecp, only:
   use field_data, only: field, fx, fy, fz
   use field_subs, only: field_calc, field_setup_old, field_calc_all
   use td_data, only: timedep, tdrestart, tdstep
   use transport_data, only : generate_rho0
   use time_dependent, only : TD
   use faint_cpu77, only: int1, int2, intsol, int3mem, int3lu
   use dftb_data, only : dftb_calc, MDFTB, MTB, chargeA_TB, chargeB_TB,        &
                         rho_aDFTB, rho_bDFTB, TBsave, TBload
   use dftb_subs, only : dftb_init, getXY_DFTB, find_TB_neighbors,             &
                         build_chimera_DFTB, extract_rhoDFT, read_rhoDFTB,     &
                         write_rhoDFTB, construct_rhoDFTB
   use cubegen       , only: cubegen_vecin, cubegen_matin, cubegen_write
   use typedef_sop   , only: sop              ! Testing SOP
   use fockbias_subs , only: fockbias_loads, fockbias_setmat, fockbias_apply
   use tmpaux_SCF    , only: neighbor_list_2e
   use liosubs_math  , only: transform
   use liosubs_dens  , only: builds_densmat, messup_densmat, standard_coefs
   use linear_algebra, only: matrix_diagon
   use converger_subs, only: converger_init, conver
   use mask_cublas   , only: cublas_setmat, cublas_release
   use typedef_operator, only: operator !Testing operator
   use maskrmm       , only: rmmget_dens, rmmCalc_init, rmmCalc_ener &
                          &, rmmCalc_core, rmmCalc_fock_cs, rmmCalc_fock_os
#  ifdef  CUBLAS
      use cublasmath , only: cumxp_r
#  endif
   use initial_guess_subs, only: get_initial_guess

!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   logical :: do_next_iteration
   logical :: is_criteria_unmet
   logical :: is_nmax_unexceeded
   integer :: Nel
   integer :: niter
   real*8  :: sq2
   real*8  :: good
   real*8  :: del
   real*8  :: DAMP0
   real*8  :: DAMP
   integer :: igpu

!  The following two variables are in a part of the code that is never
!  used. Check if these must be taken out...
   real*8  :: factor
   integer :: IDAMP

   real*8, allocatable :: rho_test(:,:)
   real*8, allocatable :: fockat(:,:)
   real*8, allocatable :: morb_coefon(:,:)

!------------------------------------------------------------------------------!
!  carlos: variables to use as input for some subroutines instead of M and NCO
   integer :: M_in
   integer :: NCOa_in
   integer :: NCOb_in

   real*8, allocatable :: rho_a0(:,:), rho_b0(:,:), rho_t0(:,:)
   real*8, allocatable :: fock_a0(:,:), fock_b0(:,:)
   real*8, allocatable :: rho_a(:,:), rho_b(:,:)
   real*8, allocatable :: fock_a(:,:), fock_b(:,:)
   real*8, allocatable :: morb_coefat(:,:)
   real*8, allocatable :: X_min(:,:)
   real*8, allocatable :: Y_min(:,:)
   real*8, allocatable :: X_min_trans(:,:)
   real*8, allocatable :: Y_min_trans(:,:)
   real*8, allocatable :: morb_energy(:)
   integer             :: i0, ii, jj, kk, kkk

!------------------------------------------------------------------------------!
! FFR variables
   type(sop)           :: overop
   real*8, allocatable :: Xmat(:,:)
   real*8, allocatable :: Ymat(:,:)
   real*8, allocatable :: Dvec(:)
   real*8, allocatable :: sqsmat(:,:)
   real*8, allocatable :: tmpmat(:,:)

   real*8, allocatable :: Hcore(:,:), Smat_i(:,:)
   real*8, allocatable :: Hraw(:,:), Hecp(:,:), Hsol(:,:), Habp(:,:)

   real*8              :: dipxyz(3)
   real*8              :: field_vec(3)

! FIELD variables (maybe temporary)
   real*8  :: Qc, Qc2, g
   integer :: ng2

!------------------------------------------------------------------------------!
! Energy contributions and convergence
   real*8, intent(inout) :: E ! Total SCF energy

   real*8 :: E1          ! kinetic + nuclear attraction + e-/MM charge
                         !    interaction + effective core potetial
   real*8 :: E1s = 0.0D0 ! kinetic + nuclear attraction + effective core
                         !    potetial
   real*8 :: E2          ! Coulomb (e- - e-)
   real*8 :: E_ecp       ! Efective core potential
   real*8 :: E_fld       ! Energy from electric field
   real*8 :: En          ! nuclear-nuclear repulsion
   real*8 :: Ens         ! MM point charge-nuclear interaction
   real*8 :: Es          ! ???
   real*8 :: E_restrain  ! distance restrain
   real*8 :: Ex          ! exchange-correlation inside SCF loop
   real*8 :: Exc         ! exchange-correlation
   real*8 :: Etrash      ! auxiliar variable
   real*8 :: Egood       !
   real*8 :: Evieja      !


! CUBLAS
!------------------------------------------------------------------------------!
   integer*8          :: dev_Xmat, dev_Ymat


!------------------------------------------------------------------------------!
! TODO : Variables to eliminate...
   real*8, allocatable :: xnano(:,:)
   integer :: MM, MM2, MMd, Md2
   integer :: M1, M2, M3, M5, M7, M9, M11, M13, M15, M17,M18,M18b, M19, M20, M22

   real*8, allocatable :: Y(:,:)
   real*8, allocatable :: Ytrans(:,:)
   real*8, allocatable :: Xtrans(:,:)

!------------------------------------------------------------------------------!
!carlos: Operators for matrices with alpha and beta spins.
   type(operator)      :: rho_aop, fock_aop
   type(operator)      :: rho_bop, fock_bop

!------------------------------------------------------------------------------!
!carlos: Open shell, variables.
   real*8              :: ocupF
   integer             :: NCOa, NCOb


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
! 0) Initialization Procedures: This part should be reviewed and possibly
!    re-structured. Many options should probably be initialized after the
!    calculations of the overlap and core matrices.
!
!------------------------------------------------------------------------------!
   call g2g_timer_sum_start('SCF')
   call g2g_timer_sum_start('Initialize SCF')

   npas=npas+1
   E=0.0D0
   E1=0.0D0
   En=0.0D0
   E2=0.0D0
   Es=0.0D0
   E_ecp=0.d0
   Ens=0.0D0
   E_restrain=0.d0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%    Distance Restrain     %%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   IF (number_restr.GT.0) THEN
      call get_restrain_energy(E_restrain)
      WRITE(*,*) "DISTANCE RESTRAIN ADDED TO FORCES"
   END IF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!------------------------------------------------------------------------------!
! TODO: RMM should no longer exist. As a first step, maybe we should
!       put all these pointers inside the maskrmm module instead of
!       having them around all subroutines...
   sq2=sqrt(2.D0)
   MM=M*(M+1)/2
   MM2=M**2
   MMd=Md*(Md+1)/2
   Md2=2*Md
   M2=2*M

!carlos: ocupation factor
   ocupF = 2.0d0
!carlos:NCOa works in open shell and close shell
   NCOa=NCO
!Number of OM down
   if (OPEN) then
      NCOb=NCO+Nunp
      ocupF=1.0d0
      allocate(rho_b0(M,M),fock_b0(M,M))
   end if

   M1=1 ! first P
   M3=M1+MM ! now Pnew
   M5=M3+MM! now S, F also uses the same position after S was used
   M7=M5+MM! now G
   M9=M7+MMd ! now Gm
   M11=M9+MMd! now H
   M13=M11+MM! W ( eigenvalues ), also this space is used in least squares
   M15=M13+M! aux ( vector for ESSl)
   M17=M15+MM! Least squares
   M18=M17+MMd! vectors of MO
   M19=M18+M*NCO! weights (in case of using option )
   M20 = M19 + natom*50*Nang ! RAM storage of two-electron integrals (if MEMO=T)

   if (OPEN) then
      M18b=M18+M*NCOa
      M19=M18b+M*NCOb
      M20 = M19 + natom*50*Nang
      M22 = M20 +2*MM !W ( beta eigenvalues )
   end if

!------------------------------------------------------------------------------!
! TODO: I don't like ending timers inside a conditional...
   if (cubegen_only) then
      call cubegen_vecin( M, NCO, MO_coef_at )
      call g2g_timer_sum_stop('Initialize SCF')
      call g2g_timer_sum_stop('SCF')
      return
   end if
!
!
! TODO: damp and gold should no longer be here??
! TODO: Qc should probably be a separated subroutine? Apparently it is only
!       used in dipole calculation so...it only adds noise to have it here.
! TODO: convergence criteria should be set at namelist/keywords setting
!
!------------------------------------------------------------------------------!
   Nel=2*NCO+Nunp

   good=1.00D0
   Egood=1.00D0
   Evieja=0.d0

   niter=0
   DAMP0=GOLD
   DAMP=DAMP0

   Qc=0.0D0
   do ii=1,natom
      Qc=Qc+Iz(ii)
   enddo
   Qc=Qc-Nel
   Qc2=Qc**2

! switch for damping/diis convergence in hybrid
   if (hybrid_converg) DIIS=.true.
!
!
!
!  1) Initializations and calculations for the fock matrix that only depend
!     on the nuclear positions. In principle, this means the nuclei-electron
!     interactions and the overlap matrix for the basis, but extra optional
!     terms may be included. Note that rmmCalc_init will also initialize
!     many variables later used by the density-dependent calculations of
!     fock in rmmCalc_fock.
!
!     Besides the nuclei-electron interaction in Hraw, the current optional
!     terms that can be added to the core are:
!
!      - Solvent (QMMM) into Hsol.
!      - Effective Core Potentials into Hecp.
!      - Atom-centered bias potential into Habp.
!
!     Contributions to the core fock matrix are kept separated for the energy
!     contribution of each and also added together in Hcore and stored into
!     RMM in the last call to rmmCalc_core.
!
!------------------------------------------------------------------------------!
   if (allocated(Smat_i)) deallocate(Smat_i)
   if (allocated(Hcore))  deallocate(Hcore)
   if (allocated(Hraw))   deallocate(Hraw)
   if (allocated(Hecp))   deallocate(Hecp)
   if (allocated(Hsol))   deallocate(Hsol)

   allocate( Smat_i(M,M) )
   allocate( Hcore(M,M) )
   allocate( Hraw(M,M) )
   allocate( Hecp(M,M) )
   allocate( Hsol(M,M) )

   Hraw(:,:) = 0.0d0
   Hecp(:,:) = 0.0d0
   Hsol(:,:) = 0.0d0

   call rmmCalc_init()
   call rmmCalc_core( Smat_i, Hraw, En, Ens, .true., .false., .false. )
   call rmmCalc_core( Smat_i, Hecp, En, Ens, .false., .true., .false. )
   call rmmCalc_core( Smat_i, Hsol, En, Ens, .false., .false., .true. )


!  TODO: Add simplified overlap diagonalization here.
!  TODO: Add the fockbias procedures to the core here.
!  TODO: It makes sense to have Hraw and Hsol calculations in the same
!        rmmCalc_core since they are very similar and are very inter-
!        connected due to g2g/aint subs. However, ECP and fockbias
!        should probably have particular subroutines.

   Hcore(:,:) = Hraw(:,:) + Hecp(:,:) + Hsol(:,:)

   call rmmCalc_core( Smat_i, Hcore, En, Ens, .false., .false., .false. )
!
!
!
!  2) Final preparations and initializations for the self consistent cycle.
!     The basic thing to do is to generate or read the initial guess, but
!     other options for different kind of SCF calculations/configurations
!     are available:
!
!      - Convergence acceleration setup
!      - DFTB hybrid calculations
!
!     All the initialization for these options should probably be here.
!------------------------------------------------------------------------------!
!
!  TODO: Add the dftb initialization procedures here.
!  TODO: Add the starting guess generation/read here.
!

!DFTB: initialisation of DFTB variables
   allocate (fock_a0(M,M), rho_a0(M,M), rho_t0(M,M) )

   M_in = M
   if (dftb_calc) then
      M_in=MDFTB
      call dftb_init(M, OPEN)
   end if

   allocate( morb_energy(M_in), morb_coefat(M_in,M_in) )
   allocate( fock_a(M_in,M_in), rho_a(M_in,M_in) )
   if (OPEN) then
      allocate( fock_b(M_in,M_in), rho_b(M_in,M_in) )
   end if

   call converger_init( M_in, ndiis, DAMP, DIIS, hybrid_converg, OPEN )
!
!
!------------------------------------------------------------------------------!
!  Unknown stuff: should it stay or should it go?
!------------------------------------------------!

   if (predcoef.and.npas.gt.3) then
      write(*,*) 'no devería estar aca!'
   end if

!
!
!
!------------------------------------------------!
! OVERLAP DIAGONALIZATION
! TODO: Simplify and move, there is just too much
!       stuff going on here.
!      (maybe trans mats are not even necessary?)
!------------------------------------------------!
   if (allocated(X_min)) deallocate(X_min)
   if (allocated(Y_min)) deallocate(Y_min)
   if (allocated(X_min_trans)) deallocate(X_min_trans)
   if (allocated(Y_min_trans)) deallocate(Y_min_trans)

   allocate(X_min(M,M), Y_min(M,M), X_min_trans(M,M), Y_min_trans(M,M))

   if ( allocated(Xmat) ) deallocate(Xmat)
   if ( allocated(Ymat) ) deallocate(Ymat)
   allocate(Xmat(M_in,M_in), Ymat(M_in,M_in))

   call overop%Sets_smat( Smat )
   if (lowdin) then
!     TODO: inputs insuficient; there is also the symetric orthog using
!           3 instead of 2 or 1. Use integer for onbasis_id
      call overop%Gets_orthog_4m( 2, 0.0d0, X_min, Y_min, X_min_trans, Y_min_trans)
   else
      call overop%Gets_orthog_4m( 1, 0.0d0, X_min, Y_min, X_min_trans, Y_min_trans)
   end if

! TODO: replace X,Y,Xtrans,Ytrans with Xmat, Ymat, Xtrp, Ytrp

   if ( allocated(Dvec) ) deallocate(Dvec)
   allocate( Dvec(M) )
   call overop%Gets_eigens_v( 0.0d0, Dvec )
   do kk = 1, M
      RMM(M13+kk-1) = Dvec(kk)
   end do
!
!
!
!------------------------------------------------!
!  Fockbias setup
   if ( allocated(sqsmat) ) deallocate(sqsmat)
   if ( allocated(tmpmat) ) deallocate(tmpmat)
   allocate( sqsmat(M,M), tmpmat(M,M) )
   call overop%Gets_orthog_2m( 3, 0.0d0, tmpmat, sqsmat )
   call fockbias_loads( natom, nuc )
   call fockbias_setmat( tmpmat )
   deallocate( sqsmat, tmpmat )


   if (allocated(Habp))   deallocate(Habp)
   allocate( Habp(M,M) )
   Habp(:,:) = 0.0d0

   call fockbias_apply( 0.0d0, Habp )
   Hcore(:,:) = Hcore(:,:) + Habp(:,:)
!
!
!
!------------------------------------------------!
!DFTB: Dimensions of Xmat and Ymat are modified for DFTB.
!
! TODO: this is nasty, a temporary solution would be to have a Msize variable
!       be assigned M (or, even better, "basis_size") or MDFTB
!       ("basis_size_dftb") according to the case

   if (dftb_calc) then
      call getXY_DFTB(M,X_min,Y_min,Xmat,Ymat)

   else
      do jj = 1, M
      do ii = 1, M
         Xmat(ii,jj) = X_min(ii,jj)
         Ymat(ii,jj) = Y_min(ii,jj)
      enddo
      enddo

   end if
!
!
!
!------------------------------------------------!
! CUBLAS
   call cublas_setmat( M_in, Xmat, dev_Xmat)
   call cublas_setmat( M_in, Ymat, dev_Ymat)
!
!
!
!------------------------------------------------!
! Generates starting guess
!
   if ( (.not.VCINP) .and. primera ) then
      call get_initial_guess(M, MM, NCO, NCOb, Xmat(MTB+1:MTB+M,MTB+1:MTB+M),  &
                             RMM(M11:MM), RMM(M1:MM), rhoalpha, rhobeta, OPEN, &
                             natom, Iz, nshell, Nuc)
      primera = .false.
   end if
!DFTB: the density for DFTB is readed from an external file.
   if (dftb_calc.and.TBload) &
   &  call read_rhoDFTB(M, MM, RMM(M1), rhoalpha, rhobeta, OPEN)

   call rmmget_dens( rho_a0 )
   rho_t0 = rho_a0

!
!
!
!------------------------------------------------!
! TD stuff that should not be here

   if ((timedep.eq.1).and.(tdrestart)) then
      call g2g_timer_sum_start('TD')
      call TD()
      call g2g_timer_sum_stop('TD')
      return
   endif
   call g2g_timer_sum_stop('Initialize SCF')
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!  3) SELF CONSISTENT CYCLE:
!
!     This is the central part of the SCF calculation. Each of the cycles
!     consists of the following set of procedures:
!     
!      - DESCRIPTION PENDING
!
!     All the initialization for these options should probably be here.
!------------------------------------------------------------------------------!
!
! TODO: Maybe evaluate conditions for loop continuance at the end of loop
!       and condense in a single "keep_iterating" or something like that.
!       Also, the WRITE_CONV_STATUS should probably go at the end.
!
! TODO: Many sections inside of the cycle should probably be 
!
   do_next_iteration = .true.
   call g2g_timer_start('Total SCF cycles')
   do while (do_next_iteration)
      call g2g_timer_start('Total iter')
      call g2g_timer_sum_start('Iteration')
      call g2g_timer_sum_start('Fock integrals')

      niter=niter+1
!
!
!     Fock contributions that depend on the density matrix
!------------------------------------------------------------------------------!
      field_vec(:) = 0.0d0
      if (field) call field_setup_old(1.0D0, 0, fx, fy, fz)
      call field_calc_all( field_vec(1), field_vec(2), field_vec(3), 0.0d0 )

      if (OPEN) then
         call rmmCalc_fock_os( rho_a0, rho_b0, Hcore, fock_a0, fock_b0, dipxyz &
                            &, E2, Ex, field, field_vec, E_fld )
      else
         call rmmCalc_fock_cs( rho_a0, Hcore, fock_a0, dipxyz                  &
                            &, E2, Ex, field, field_vec, E_fld )
      end if

!
!
!
!------------------------------------------------------------------------------!
! REACTION FIELD CASE
!
! TODO: I'm gonna be commenting this for the following two reasons:
!
!   (1) Implementation is weird and limiting. Apparently you can only
!       SCF with a field when computing the rho0 for a transport run.
!
!   (2) field_calc does not return the contributions to the fock matrix
!       but leaves everything inside RMM. The current version of Calc_fock
!       returns this thus allowing the external function to ignore RMM.
!
!   We can later discuss if the calculation of the field contributions
!   to the fock matrix belongs in a field module or in a fock related
!   module (I would advocate for the second since it needs much more info
!   from the structure of fock whereas only needs to know the 3 values of
!   the current field), or wether it should be separated from the other
!   calculations of density-dependent fock terms, currently in rmmCalc_fock.
!
!------------------------------------------------------------------------------!
!
!        if ( generate_rho0 ) then
!           if (field) call field_setup_old(1.0D0, 0, fx, fy, fz)
!           call field_calc(E_fld, 0.0D0)
!        end if
!
!------------------------------------------------------------------------------!
!       E1 includes all contributions to core
        E1=0.0D0
        call rmmCalc_ener( M, Hcore, rho_t0, E1 )

        call g2g_timer_start('actualiza rmm')
        call g2g_timer_sum_pause('Fock integrals')
!
!
!
!------------------------------------------------------------------------------!
! DFTB: Fock and Rho for DFTB are builded.
!
! TODO: this should be wrapped inside a single dftb subroutine. Also, two
!       consecutive dftb_calc switches? really?
!
      if (dftb_calc) then

         NCOa_in = NCOa + MTB
         if (niter==1) call find_TB_neighbors(M,Nuc,natom)
         call build_chimera_DFTB (M, fock_a0, fock_a, natom, nshell, ncont)
         call construct_rhoDFTB(M, rho_a, rho_a0 ,rho_aDFTB, TBload, niter)

         if (OPEN) then
            NCOb_in = NCOb + MTB
            call build_chimera_DFTB (M, fock_b0, fock_b, natom, nshell, ncont)
            call construct_rhoDFTB(M, rho_b, rho_b0 ,rho_bDFTB, TBload, niter)
         end if

      else
         NCOa_in = NCOa
         fock_a=fock_a0
         rho_a=rho_a0

         if (OPEN) then
            NCOb_in = NCOb
            fock_b=fock_b0
            rho_b=rho_b0
         end if

      endif
!carlos: storing rho and fock in operator.

      call rho_aop%Sets_data_AO(rho_a)
      call fock_aop%Sets_data_AO(fock_a)
      if (OPEN) then
         call rho_bop%Sets_data_AO(rho_b)
         call fock_bop%Sets_data_AO(fock_b)
      end if
!------------------------------------------------------------------------------!
!  Convergence accelerator processing
      call g2g_timer_sum_start('SCF acceleration')

!carlos: this is repeted twice for open shell

!%%%%%%%%%%%%%%%%%%%%
!CLOSE SHELL OPTION |
!%%%%%%%%%%%%%%%%%%%%
#     ifdef CUBLAS
         call conver( niter, good, good_cut, M_in, rho_aop, fock_aop           &
                   &, dev_Xmat, dev_Ymat, 1)
#     else
         call conver( niter, good, good_cut, M_in, rho_aop, fock_aop           &
                   &, Xmat, Ymat, 1)
#     endif

      call g2g_timer_sum_pause('SCF acceleration')

!------------------------------------------------------------------------------!
!  Fock(ON) diagonalization
      if ( allocated(morb_coefon) ) deallocate(morb_coefon)
      allocate( morb_coefon(M_in,M_in) )

      call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
      call fock_aop%Diagon_datamat( morb_coefon, morb_energy )
      call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')
!
!
!------------------------------------------------------------------------------!
!  Base change of coeficients ( (X^-1)*C ) and construction of new density
!  matrix
      call g2g_timer_sum_start('SCF - MOC base change (sum)')
#     ifdef CUBLAS
         call cumxp_r( morb_coefon, dev_Xmat, morb_coefat, M_in)
#     else
         morb_coefat = matmul( Xmat, morb_coefon )
#     endif
      call standard_coefs( morb_coefat )
      call g2g_timer_sum_pause('SCF - MOC base change (sum)')

      if ( allocated(morb_coefon) ) deallocate(morb_coefon)
      call rho_aop%Dens_build(M_in, NCOa_in, ocupF, morb_coefat)
      call rho_aop%Gets_data_AO(rho_a)


!carlos: Alpha Energy (or Close Shell) storage in RMM

      do kk=1,M
         RMM(M13+kk-1) = morb_energy(kk)
      end do

      i0 = 0
      if (dftb_calc) i0=MTB

      kkk = 0
      do kk=1,NCOa
      do ii=1,M
         kkk = kkk+1
         MO_coef_at(kkk) = morb_coefat( i0+ii, kk )
      enddo
      enddo

!%%%%%%%%%%%%%%%%%%%%
!OPEN SHELL OPTION  |
!%%%%%%%%%%%%%%%%%%%%
      if (OPEN) then

#        ifdef CUBLAS
            call conver( niter, good, good_cut, M_in, rho_bop, fock_bop        &
                      &, dev_Xmat, dev_Ymat, 2)
#        else
            call conver( niter, good, good_cut, M_in, rho_bop, fock_bop        &
                      &, Xmat, Ymat, 2)
#        endif
         call g2g_timer_sum_pause('SCF acceleration')

!------------------------------------------------------------------------------!
!  Fock(ON) diagonalization
         if ( allocated(morb_coefon) ) deallocate(morb_coefon)
         allocate( morb_coefon(M_in,M_in) )

         call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
         call fock_bop%Diagon_datamat( morb_coefon, morb_energy )
         call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')
!
!
!------------------------------------------------------------------------------!
!  Base change of coeficients ( (X^-1)*C ) and construction of new density
!  matrix
         call g2g_timer_sum_start('SCF - MOC base change (sum)')
#        ifdef CUBLAS
            call cumxp_r( morb_coefon, dev_Xmat, morb_coefat, M_in)
#        else
            morb_coefat = matmul( Xmat, morb_coefon )
#        endif
         call standard_coefs( morb_coefat )
         call g2g_timer_sum_pause('SCF - MOC base change (sum)')

         if ( allocated(morb_coefon) ) deallocate(morb_coefon)
         call rho_bop%Dens_build(M_in, NCOb_in, ocupF, morb_coefat)
         call rho_bop%Gets_data_AO(rho_b)

!carlos: Beta Energy storage in RMM
         do kk=1,M
            RMM(M22+kk-1) = morb_energy(kk)
         end do
!carlos: Storing autovectors to create the restart
         i0 = 0
         if (dftb_calc) i0=MTB

         kkk = 0
         do kk=1,NCOb
         do ii=1,M
            kkk = kkk+1
            MO_coef_at_b(kkk) = morb_coefat( i0+ii, kk )
         enddo
         enddo

      end if
!------------------------------------------------------------------------------!
!carlos: storing matrices
!------------------------------------------------------------------------------!
!  We are not sure how to translate the sumation over molecular orbitals
!  and energies when changing from DFTB system to DFT subsystem. Forces
!  may be broken due to this.
!  This should not be affecting normal DFT calculations.

!carlos: this is not working with openshell. morb_coefat is not depending of
!        the spin factor.
      do ii=1,M
!        do jj=1,NCO
         do jj=1,M
            X( ii, M2+jj ) = morb_coefat( i0+ii, jj )
         enddo
!         do jj=NCO,M
!            X( ii, M2+jj ) = morb_coefat( i0+ii, jj+2*MTB )
!         enddo
      enddo


!------------------------------------------------------------------------------!
! carlos: added to separate from rho the DFT part
!
! TODO: again, this should be handled differently...
!       ex: extract_rhoDFT should still put rho_a0 into rho_a even if
!       no dftb calc is being executed.
!
      if (dftb_calc) then
         rho_aDFTB(:,:) = rho_a(:,:)
         call extract_rhoDFT(M, rho_a, rho_a0)
         if (OPEN) then
            rho_bDFTB(:,:) = rho_b(:,:)
            call extract_rhoDFT(M, rho_b, rho_b0)
         end if
      else
         rho_a0(:,:) = rho_a(:,:)
         if (OPEN) rho_b0(:,:) = rho_b(:,:)
      end if

      if (OPEN) then
         rho_t0(:,:) = rho_a0(:,:) + rho_b0(:,:)
      else
         rho_t0(:,:) = rho_a0(:,:)
      end if
!
!
!
!------------------------------------------------------------------------------!
! TODO: convergence criteria should be a separated subroutine.
! TODO: make xnano go away...only remains here
! TODO: finalization of the loop is still a little bit messy.
!
      allocate( xnano(M,M) )
      if (OPEN) then
         xnano(:,:) = rho_a0(:,:)
         call messup_densmat( xnano )
         call sprepack( 'L', M, rhoalpha, xnano)

         xnano(:,:) = rho_b0(:,:)
         call messup_densmat( xnano )
         call sprepack( 'L', M, rhobeta,  xnano)

         xnano(:,:) = rho_a0(:,:) + rho_b0(:,:)
         call messup_densmat( xnano )
      else
         xnano(:,:) = rho_a0(:,:)
         call messup_densmat( xnano )

      end if

      good = 0.0d0
      do jj=1,M
      do kk=jj,M
         del=xnano(jj,kk)-(RMM(kk+(M2-jj)*(jj-1)/2))
         del=del*sq2
         good=good+del**2
         RMM(kk+(M2-jj)*(jj-1)/2)=xnano(jj,kk)
      enddo
      enddo
      good=sqrt(good)/float(M)
      deallocate ( xnano )

! Damping factor update
      DAMP=DAMP0

      E=E1+E2+En
!      E=E+Es
!
!     Writes energy at every step
      if (verbose) call WRITE_E_STEP(niter, E+Ex)

      Egood=abs(E+Ex-Evieja)
      Evieja=E+Ex

!     Writes convergence criteria and the value of the dyn step
      if (verbose) call WRITE_CONV_STATUS( GOOD, TOLD, EGOOD, ETOLD )

      is_criteria_unmet = (good.ge.told).or.(Egood.ge.Etold)
      is_nmax_unexceeded = (niter.le.NMAX)
      do_next_iteration = (is_criteria_unmet).and.(is_nmax_unexceeded)

      call g2g_timer_stop('Total iter')
      call g2g_timer_sum_pause('Iteration')

   end do
   call g2g_timer_stop('Total SCF cycles')
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!
!
!  4) DESCRIPTION PENDING (finalize SCF)
!     
!------------------------------------------------------------------------------!
   call g2g_timer_sum_start('Finalize SCF')


!  Checks of convergence
   if (niter.ge.NMAX) then
      write(6,*) 'NO CONVERGENCE AT ',NMAX,' ITERATIONS'
      noconverge=noconverge + 1
      converge=0
   else
      write(6,*) 'CONVERGED AT',niter,'ITERATIONS'
      noconverge = 0
      converge=converge+1
   endif

   if (noconverge.gt.4) then
      write(6,*)  'stop for not convergion 4 times'
      stop
   endif


!DFTB: Mulliken analysis of TB part
   if (dftb_calc) then

      chargeA_TB=MTB
      chargeB_TB=MTB
      do ii=1, MTB
         chargeA_TB=chargeA_TB-rho_a(ii,ii)
         chargeB_TB=chargeB_TB-rho_a(MTB+M+ii,MTB+M+ii)
      end do

      if (OPEN) then
         do ii=1,MTB
            chargeA_TB=chargeA_TB-rho_b(ii,ii)
            chargeB_TB=chargeB_TB-rho_b(MTB+M+ii,MTB+M+ii)
         end do
      end if

      open(unit=6800,file='TB_Mulliken')
      write(6800,*) "Mulliken charge of part A", chargeA_TB
      write(6800,*) "Mulliken charge of part B", chargeB_TB

   end if


!DFTB: The last rho is stored in an output as a restart.
   if (dftb_calc.and.TBsave) call write_rhoDFTB(M_in, OPEN)
!
!
!
! TODO: Comments about a comented sections? Shouldn't it all of this go away?
!
!    CH - Why call intsol again here? with the .false. parameter,
!    E1s is not recalculated, which would be the only reason to do
!    this again; Ens isn't changed from before...
! -- SOLVENT CASE --------------------------------------
!      if (sol) then
!      call g2g_timer_sum_start('intsol 2')
!      if(nsol.gt.0) then
!        call intsol(E1s,Ens,.false.)
!        write(*,*) 'cosillas',E1s,Ens
!        call g2g_timer_sum_stop('intsol 2')
!      endif
!      call mmsol(natom,Nsol,natsol,Iz,pc,r,Em,Rm,Es)
   Es=Es+E1s+Ens
!     endif


   if (MOD(npas,energy_freq).eq.0) then
   if (GRAD) then

!     Resolve with last density to get XC energy
      call g2g_timer_sum_start('Exchange-correlation energy')
      call g2g_new_grid(igrid)
      call g2g_solve_groups(1, Exc, 0)
      call g2g_timer_sum_stop('Exchange-correlation energy')

!     COmputing the QM/MM contribution to total energy
!     Total SCF energy =
!     E1    - kinetic + nuclear attraction + QM/MM interaction + effective
!            core potential
!     E2    - Coulomb
!     En    - nuclear-nuclear repulsion
!     Ens   - MM point charge - nuclear interaction
!     Exc   - exchange-correlation
!     E_ecp - Efective core potential
!     E_restrain - distance restrain

!     NucleusQM-CHarges MM
      Es=Ens

!     One electron Kinetic (with aint >3) or Kinetic + Nuc-elec (aint >=3)
      call int1(En)

!     Computing the E1-fock without the MM atoms
      if (nsol.gt.0.and.igpu.ge.1) then
        call aint_qmmm_init(0,r,pc)
        call aint_qmmm_fock(E1s,Etrash)
        call aint_qmmm_init(nsol,r,pc)
      endif

!     E1s (here) is the 1e-energy without the MM contribution
      E1s=0.D0
      do kk=1,MM
        E1s=E1s+RMM(kk)*RMM(M11+kk-1)
      enddo

!     Es is the QM/MM energy computated as total 1e - E1s + QMnuc-MMcharges
      Es=Es+E1-E1s

!     Part of the QM/MM contrubution are in E1
      E=E1+E2+En+Ens+Exc+E_restrain

!     Write Energy Contributions
      if (npas.eq.1) npasw = 0

      if (npas.gt.npasw) then
         E_ecp = 0.0d0
         call rmmCalc_ener( M, Hecp, rho_t0, E_ecp )
         Es = Es - E_ecp
         call WriteEnergies( E1, E2, En, Ens, E_ecp, Exc, ecpmode, E_restrain )
         npasw = npas+10
      end if

   endif ! GRAD
   endif ! npas
!
!
! calculation of energy weighted density matrix
!
! TODO: this should be a separated subroutine...
!
   call g2g_timer_sum_start('energy-weighted density')
   kkk=0
   do jj=1,M
   do ii=jj,M
      kkk=kkk+1
      RMM(M15+kkk-1)=0.D0
!
      if (ii.eq.jj) then
         factor=2.D0
      else
         factor=4.D0
      endif

      do kk=1,NCO
         RMM(M15+kkk-1)= &
         RMM(M15+kkk-1)-RMM(M13+kk-1)*factor*X(ii,M2+kk)*X(jj,M2+kk)
      enddo
   enddo
   enddo

   call g2g_timer_sum_stop('energy-weighted density')

!  Variables needed for further calculations (Populations, Dip, etc).
   Enucl = En
   do kkk=1, M
       Eorbs(kkk) = RMM(M13+kkk-1)
   enddo

   call cubegen_matin( M, X )
!
!
! TODO: have ehrendyn call SCF and have SCF always save the resulting rho in
!       a module so that ehrendyn can retrieve it afterwards.
!       Remove all of this.
!
   if (doing_ehrenfest) then
      call spunpack('L',M,RMM(M1),RealRho)
      call fixrho(M,RealRho)
      call ehrendyn_init(natom, M, RealRho)
   endif
!
!
! TODO: Deallocation of variables that should be removed
! TODO: MEMO should be handled differently...
   if (MEMO) then
     deallocate(kkind,kkinds)
     deallocate(cool,cools)
   endif
!
!
! TODO: Why is TD being called again?? In any case have TD call SCF as a first
!       step to obtain the first density...
!
   if (timedep.eq.1) then
     call g2g_timer_sum_start('TD')
     call TD()
     call g2g_timer_sum_stop('TD')
   endif
!
!
! TODO: Cublas should be handled differently. Hidden behind SOP type and an
!       interface or general endstep call that takes care of cublas_shutdown
!       and all other similar stuff.
!
   call cublas_release( dev_Xmat )
   call cublas_release( dev_Ymat )
   call cublas_release( )

   call g2g_timer_sum_stop('Finalize SCF')
   call g2g_timer_sum_stop('SCF')

end subroutine SCF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
