module converger_subs

   implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine converger_init( M_in, ndiis_in, factor_in, do_diis, do_hybrid, OPshell )
   use converger_data, only: fockm, FP_PFm, conver_criter, fock_damped, &
                             hagodiis, damping_factor, bcoef, ndiis, EMAT2

   implicit none
   double precision, intent(in) :: factor_in
   integer         , intent(in) :: M_in, ndiis_in
   logical         , intent(in) :: do_diis, do_hybrid, OPshell

   hagodiis       = .false.
   damping_factor = factor_in
   ndiis          = ndiis_in

   if (do_hybrid) then
      conver_criter = 3
   else if (do_diis) then
      conver_criter = 2
   else
      conver_criter = 1
   endif

   ! Added to change from damping to DIIS. - Nick
   if (conver_criter /= 1) then
      if(OPshell) then
         if (.not. allocated(fockm)  ) allocate(fockm (M_in, M_in, ndiis, 2))
         if (.not. allocated(FP_PFm) ) allocate(FP_PFm(M_in, M_in, ndiis, 2))
         if (.not. allocated(bcoef)  ) allocate(bcoef (ndiis+1, 2) )
         if (.not.allocated(EMAT2)   ) allocate(EMAT2(ndiis+1,ndiis+1,2))
      else
         if (.not. allocated(fockm) )  allocate(fockm (M_in, M_in, ndiis, 1))
         if (.not. allocated(FP_PFm) ) allocate(FP_PFm(M_in, M_in, ndiis, 1))
         if (.not. allocated(bcoef) )  allocate(bcoef (ndiis+1, 1))
         if (.not.allocated(EMAT2) )   allocate(EMAT2(ndiis+1,ndiis+1,1))
      end if
   endif

   if(OPshell) then
      if (.not. allocated(fock_damped) ) allocate(fock_damped(M_in, M_in, 2))
   else
      if (.not. allocated(fock_damped) ) allocate(fock_damped(M_in, M_in, 1))
   end if
   fock_damped(:,:,:) = 0.0D0
end subroutine converger_init

   subroutine conver (niter, good, good_cut, M_in, rho_op, fock_op, &
#ifdef CUBLAS
                      devPtrX, devPtrY, spin)
#else
                      Xmat, Ymat, spin)
#endif
   use converger_data  , only: damping_factor, hagodiis, fockm, FP_PFm, ndiis, &
                               fock_damped, bcoef, EMAT2, conver_criter
   use typedef_operator, only: operator
   use fileio_data     , only: verbose
   use linear_algebra  , only: matmuldiag

   implicit none
   ! Spin allows to store correctly alpha or beta information. - Carlos
   integer         , intent(in)    :: niter, M_in, spin
   double precision, intent(in)    :: good, good_cut
   type(operator)  , intent(inout) :: rho_op, fock_op

#ifdef  CUBLAS
   integer*8       , intent(in) :: devPtrX, devPtrY
#else
   double precision, intent(in) :: Xmat(M_in,M_in), Ymat(M_in,M_in)
#endif

   double precision :: damp
   integer          :: ndiist, ii, jj, kk, kknew, lwork, info
   double precision, allocatable :: fock00(:,:), EMAT(:,:), diag1(:,:),      &
                                    suma(:,:), scratch1(:,:), scratch2(:,:), &
                                    fock(:,:), rho(:,:), work(:)


! INITIALIZATION
! If DIIS is turned on, update fockm with the current transformed F' (into ON
! basis) and update FP_PFm with the current transformed [F',P']
!
! (1)     Calculate F' and [F',P']
!       update fockm with F'
! now, scratch1 = A = F' * P'; scratch2 = A^T
! [F',P'] = A - A^T
! BASE CHANGE HAPPENS INSIDE OF FOCK_COMMUTS

   allocate(fock00(M_in,M_in), fock(M_in,M_in), rho(M_in,M_in), work(1000))

   ! Saving rho and the first fock AO
   call rho_op%Gets_data_AO(rho)
   call fock_op%Gets_data_AO(fock00)

   ndiist = min( niter, ndiis )
   if (conver_criter /= 1) then
      allocate( suma(M_in, M_in), diag1(M_in, M_in) )
      allocate( scratch1(M_in, M_in), scratch2(M_in, M_in) )


! If DIIS is turned on, update fockm with the current transformed F' (into ON
! basis) and update FP_PFm with the current transformed [F',P']

      do jj = ndiis-(ndiist-1), ndiis-1
         fockm(:,:,jj,spin)  = fockm(:,:,jj+1,spin)
         FP_PFm(:,:,jj,spin) = FP_PFm(:,:,jj+1,spin)
      enddo

#ifdef CUBLAS
      call rho_op%BChange_AOtoON(devPtrY, M_in, 'r')
      call fock_op%BChange_AOtoON(devPtrX, M_in, 'r')
#else
      call rho_op%BChange_AOtoON(Ymat, M_in, 'r')
      call fock_op%BChange_AOtoON(Xmat,M_in, 'r')
#endif
      call rho_op%Gets_data_ON(rho)
      call fock_op%Commut_data_r(rho, scratch1, M_in)

      FP_PFm(:,:,ndiis,spin) = scratch1(:,:)
      call fock_op%Gets_data_ON( fockm(:,:,ndiis,spin) )

   endif

   select case (conver_criter)
      ! Always do damping
      case (1)
         hagodiis = .false.

      ! Damping the first two steps, diis afterwards
      case (2)
         if (niter > 2) then
            hagodiis = .true.
         else
            hagodiis = .false.
         endif

      ! Damping until good enough, diis afterwards
      case(3)
         if (good < good_cut) then
            if ( (.not. hagodiis) .and. (verbose .gt. 3) ) then
               write(6,'(A,I4)') "  Changing to DIIS at step: ", niter
            endif
            hagodiis=.true.
         endif

      case default
         write(*,'(A,I4)') 'ERROR - Wrong conver_criter = ', conver_criter
         stop
   endselect

   ! THIS IS DAMPING
   ! If we are not doing diis this iteration, apply damping to F, save this
   ! F in fock_damped for next iteration's damping and put F' = X^T * F * X in
   ! fock the newly constructed damped matrix is stored, for next iteration in
   ! fock_damped
   if (.not.hagodiis) then
      fock = fock00

      if (niter > 1) &
         fock = (fock  + damping_factor * fock_damped(:,:,spin)) / &
                (1.0D0 + damping_factor)
      fock_damped(:,:,spin) = fock
      call fock_op%Sets_data_AO(fock)

#ifdef  CUBLAS
      call fock_op%BChange_AOtoON(devPtrX, M_in, 'r')
#else
      call fock_op%BChange_AOtoON(Xmat   , M_in, 'r')
#endif
   endif

   ! DIIS
   if (conver_criter /= 1) then
      allocate(EMAT(ndiist+1,ndiist+1))

      ! Before ndiis iterations, we just start from the old EMAT
      if ((niter .gt. 1) .and. (niter .le. ndiis)) then
         EMAT = 0.0D0
         do jj = 1, ndiist-1
         do ii = 1, ndiist-1
            EMAT(ii,jj) = EMAT2(ii,jj,spin)
         enddo
         enddo
      ! After ndiis iterations, we start shifting the oldest iteration stored
      else if (niter.gt.ndiis) then
         EMAT = 0.0D0
         do jj = 1, ndiist-1
         do ii = 1, ndiist-1
            EMAT(ii,jj) = EMAT2(ii+1,jj+1,spin)
         enddo
         enddo
      endif

      ! scratch1 and scratch2 store the commutations from different iterations.
      do kk = 1, ndiist
         kknew = kk + (ndiis - ndiist)
         scratch1(:,:) = FP_PFm(:,:,ndiis,spin)
         scratch2(:,:) = FP_PFm(:,:,kknew,spin)

         call matmuldiag( scratch1, scratch2, diag1, M_in )
         EMAT(ndiist,kk) = 0.0d0

         if (kk.ne.ndiist) EMAT(kk,ndiist) = 0.0d0
         do ii = 1, M_in
            EMAT(ndiist,kk) = EMAT(ndiist,kk) + diag1(ii,ii)
            if (kk.ne.ndiist) then
               EMAT(kk,ndiist) = EMAT(ndiist,kk)
            endif
         enddo
      enddo

      do kk = 1, ndiist
         EMAT(kk,ndiist+1) = -1.0d0
         EMAT(ndiist+1,kk) = -1.0d0
      enddo
      EMAT(ndiist+1, ndiist+1)= 0.0d0
      EMAT2(1:ndiist+1,1:ndiist+1,spin) = EMAT

      !   THE MATRIX EMAT SHOULD HAVE THE FOLLOWING SHAPE:
      !      |<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0|
      !      |<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0|
      !      |<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0|
      !      |<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0|
      !      |     .            .      ...     . |
      !      |   -1.0         -1.0     ...    0. |
      !   WHERE <E(I)*E(J)> IS THE SCALAR PRODUCT OF [F*P] FOR ITERATION I
      !   TIMES [F*P] FOR ITERATION J.

      if (hagodiis) then
         do ii = 1, ndiist
            bcoef(ii,spin) = 0.0d0
         enddo
         bcoef(ndiist+1,spin) = -1.0d0

         ! First call to DGELS sets optimal WORK size. Second call solves the
         ! A*X = B (EMAT * Ci = bCoef) problem, with bCoef also storing the
         ! result.
         LWORK = -1
         CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
                     ndiist+1, bcoef(:,spin), ndiist+1, WORK, LWORK, INFO )

         LWORK = MIN( 1000, INT( WORK( 1 ) ) )
         CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
                     ndiist+1, bcoef(:,spin), ndiist+1, WORK, LWORK, INFO )

         ! Build new Fock as a linear combination of previous steps.
         suma = 0.0D0
         do kk=1,ndiist
            kknew = kk + (ndiis - ndiist)
            do ii = 1, M_in
            do jj = 1, M_in
               suma(ii,jj) = suma(ii,jj) + bcoef(kk,spin) * &
                                           fockm(ii,jj,kknew,spin)
            enddo
            enddo
         enddo
         fock = suma
         call fock_op%Sets_data_ON(fock)

      endif
   endif
   deallocate (work)
end subroutine conver

end module converger_subs
