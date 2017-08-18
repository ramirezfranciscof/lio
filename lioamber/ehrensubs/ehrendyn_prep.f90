!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn_prep
&  ( time, nbasis, natoms, nucpos, nucvel, Sinv, Uinv, Linv, Rmon,             &
&    Fmat, Dmat, force_ds, dipmom, energy )
!------------------------------------------------------------------------------!
!
!  Rmon - Density matrix enters in the ON basis
!  Fmat - Fock matrix enters in the AO basis, and exits in the ON basis
!
!  Fock and Force calculation need density and fock in AO, but the propagation
!  needs Fock in
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   real*8    , intent(in)    :: time
   integer   , intent(in)    :: nbasis
   integer   , intent(in)    :: natoms

   real*8    , intent(in)    :: nucpos(3, natoms)
   real*8    , intent(in)    :: nucvel(3, natoms)
   real*8    , intent(in)    :: Sinv(nbasis, nbasis)
   real*8    , intent(in)    :: Uinv(nbasis, nbasis)
   real*8    , intent(in)    :: Linv(nbasis, nbasis)
   complex*16, intent(in)    :: Rmon(nbasis, nbasis)

   real*8    , intent(inout) :: Fmat(nbasis, nbasis)
   real*8    , intent(inout) :: Dmat(nbasis, nbasis)
   real*8    , intent(inout) :: force_ds(3, natoms)
   real*8    , intent(inout) :: dipmom(3)
   real*8    , intent(inout) :: energy

   real*8                    :: elec_field(3)
   real*8    , allocatable   :: Bmat(:,:)
   complex*16, allocatable   :: Rmao(:,:)


   allocate( Bmat(nbasis,nbasis), Rmao(nbasis,nbasis) )

!  Fock and Force calculation (needs density and fock in AO)
!  (this should leave the right Rho in RMM for get_forces)
   Rmao = Rmon
   Rmao = matmul(Rmao, Linv)
   Rmao = matmul(Uinv, Rmao)

   call ehrenaux_setfld( time, elec_field )
   call RMMcalc3_FockMao( dens_mao, elec_field, fock_mid, dipmom, energy)
   call calc_forceDS( natoms, nbasis, nucpos, nucvel, dens_mao, fock_mid, Sinv,&
                    & Bmat, qm_forces_ds )


!  Set ups propagation cuasi-fock matrix (needs fock in ON)
   fock_mid = matmul(fock_mid, Uinv)
   fock_mid = matmul(Linv, fock_mid)
   Dmat = calc_Dmat( nbasis, Linv, Uinv, Bmat )


   deallocate( Bmat, Dmat)
end subroutine ehrenstep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
