!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% INT3LU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! 2e integral gradients, 3 indexes: density fitting functions and wavefunction.!
! Calculates Coulomb elements for the Fock matrix, and 2e energy.              !
!                                                                              !
! EXTERNAL INPUT: system information.                                          !
!   · rho(M,M): density matrix.                                                !
!   · Fmat(M,M): Fock matrix (Fock alpha in open shell).                       !
!   · Fmat_b(M,M): Fock beta matrix (ignored in closed shell).                 !
!   · Gmat(M,M): Coulomb G matrix.                                             !
!   · Ginv(M,M): Inverted coulomb G matrix.                                    !
!   · Hmat(M,M): 1e matrix elements.                                           !
!   · open_shell: boolean indicating open-shell calculation.                   !
!                                                                              !
! INTERNAL INPUT: basis set information.                                       !
!   · M: number of basis functions (without contractions)                      !
!   · Md: number of auxiliary basis functions (without contractions)           !
!   · af(Md): variational coefficient for auxiliary function i.                !
!   · MEMO: indicates if cool/kkind/kknum are stored in memory. This is not    !
!           used when performing analytic integrals in GPU.                    !
!   · cool: precalculated 2e terms in double precision.                        !
!   · kkind: precalculated indexes for double precision Fock matrix elements.  !
!   · kknumd: number of precalculated double precision Fock matrix elements.   !
!   · cools: precalculated 2e terms in single precision.                       !
!   · kkinds: precalculated indexes for single precision Fock matrix elements. !
!   · kknums: number of precalculated single precision Fock matrix elements.   !
!                                                                              !
! EXTERNAL OUTPUTS:                                                            !
!   · E2: 2e coulomb energy.                                                   !
!                                                                              !
! Original and debugged (or supposed to): Dario Estrin Jul/1992                !
! Refactored:                             Federico Pedron Sep/2018             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module subm_int3lu
contains
subroutine int3lu(E2, rho, Fmat_b, Fmat, Gmat, Ginv, Hmat, open_shell, memo)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutines - 2e integrals, 3 index                                !
! Wavefunction and density fitting functions are calculated using the          !
! Obara-Saika recursive method.                                                !
! Inputs: G, F, standard basis and density basis.                              !
! F should already have the 1e part, and here the Coulomb part is added without!
! storing the integrals separately.                                            !
! Output: F updated with Coulomb part, also Coulomb energy.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use basis_data, only: M, Md, cool, cools, kkind, kkinds, kknumd, kknums, &
                         af, MM, MMd

   implicit none
   logical         , intent(in) :: open_shell, memo
   double precision, intent(in) :: rho(:), Gmat(:), Ginv(:), Hmat(:)
   double precision, intent(inout) :: E2, Fmat_b(:), Fmat(:)

   double precision, allocatable :: Rc(:), aux(:)
   double precision :: Ea, Eb, term
   integer          :: ll(3), iikk, k_ind, kk_ind, m_ind

   ! 16 loops for all combinations - 1-2: for wavefunction basis, 3 for the
   ! density fitting.
   ! Rc(k) is constructed adding t(i,j,k)*P(i,j), and cf(k), the variationally
   ! obtained fitting coefficient, is obtained by adding R(i)*G-1(i,k)
   ! if t(i,j,k) is not stored, they should be calculated again in order to
   ! evaluate the corresponding part of the Fock matrix.
   ! V(i,j) is obtained by adding af(k_ind) * t(i,j,k).
   allocate(Rc(Md), aux(md))
   Ea = 0.D0 ; Eb = 0.D0

   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2

   if (MEMO) then
      call g2g_timer_start('int3lu - start')

      do k_ind = 1, 3
         Ll(k_ind) = k_ind * (k_ind-1) / 2
      enddo

      do k_ind = 1, Md
         Rc(k_ind) = 0.D0
      enddo

      do kk_ind = 1, kknumd
         iikk = (kk_ind - 1) * Md
         do k_ind = 1, Md
            Rc(k_ind) = Rc(k_ind) + rho(kkind(kk_ind)) * cool(iikk + k_ind)
         enddo
      enddo

      do kk_ind = 1, kknums
         iikk = (kk_ind - 1) * Md
         do k_ind = 1, Md
            Rc(k_ind) = Rc(k_ind) + rho(kkinds(kk_ind)) * cools(iikk + k_ind)
         enddo
      enddo

      ! Calculation of variational coefficients and fitting coefficients
      do m_ind = 1, Md
         af(m_ind) = 0.0D0
         do k_ind = 1, m_ind-1
            af(m_ind) = af(m_ind) + &
                        Rc(k_ind) * Ginv(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind, Md
            af(m_ind) = af(m_ind) + &
                        Rc(k_ind) * Ginv(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo

      ! Initialization of Fock matrix elements
      do k_ind = 1, MM
         Fmat(k_ind) = Hmat(k_ind)
      enddo
      if (open_shell) then
      do k_ind = 1, MM
         Fmat_b(k_ind) = Hmat(k_ind)
      enddo
      endif

      do m_ind = 1, Md
         Ea = Ea + af(m_ind)  * Rc(m_ind)
         do k_ind = 1, m_ind
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind+1, Md
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo

      ! Calculation of all integrals again, in order to build the Fock matrix.
      aux = 0.0D0
      if (open_shell) then
         do k_ind = 1, Md
            aux(k_ind) = af(k_ind)
         enddo
      endif

      call g2g_timer_stop('int3lu - start')
      call g2g_timer_start('int3lu')
      if (open_shell) then
         do kk_ind = 1, kknumd
            iikk = (kk_ind - 1) * Md
            do k_ind = 1, Md
               Fmat(kkind(kk_ind)) = Fmat(kkind(kk_ind)) + &
                                     af(k_ind)  * cool(iikk + k_ind)
               Fmat_b(kkind(kk_ind)) = Fmat_b(kkind(kk_ind)) + &
                                       aux(k_ind) * cool(iikk + k_ind)
            enddo
         enddo
         do kk_ind = 1, kknums
            iikk = (kk_ind - 1) * Md
            do k_ind = 1, Md
               Fmat(kkinds(kk_ind)) = Fmat(kkinds(kk_ind)) + &
                                      af(k_ind)  * cools(iikk + k_ind)
               Fmat_b(kkinds(kk_ind)) = Fmat_b(kkinds(kk_ind)) + &
                                        aux(k_ind) * cools(iikk + k_ind)
            enddo
         enddo
      else
         do kk_ind = 1, kknumd
            iikk = (kk_ind - 1) * Md
            term = 0.0D0
            do k_ind = 1, Md
              term = term + af(k_ind) * cool(iikk + k_ind)
            enddo
            Fmat(kkind(kk_ind)) = Fmat(kkind(kk_ind)) + term
         enddo
         do kk_ind = 1, kknums
            iikk = (kk_ind - 1) * Md
            term = 0.0D0
            do k_ind = 1, Md
              term = term + af(k_ind) * cools(iikk + k_ind)
            enddo
            Fmat(kkinds(kk_ind)) = Fmat(kkinds(kk_ind)) + term
         enddo
      endif
      call g2g_timer_stop('int3lu')
   else
      do k_ind = 1, MM
         Fmat(k_ind) = Hmat(k_ind)
         if (open_shell) Fmat_b(k_ind) = Hmat(k_ind)
      enddo

      call aint_coulomb_fock(Ea)
      do m_ind = 1, Md
         do k_ind = 1, m_ind
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind+1, Md
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo
   endif

   E2 = Ea - Eb / 2.D0
   deallocate(Rc, aux)
   return
end subroutine int3lu
end module subm_int3lu
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
