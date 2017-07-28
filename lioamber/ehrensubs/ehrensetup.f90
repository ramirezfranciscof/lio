!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrensetup( Natom, Nbasis, real_dens, time_step )
   use liokeys,    only: edyn_steps
   use ehrendata,  only: dt_nucl, dt_elec, RhoSaveA, RhoSaveB
   implicit none
   integer, intent(in) :: Natom
   integer, intent(in) :: Nbasis
   real*8,  intent(in) :: real_dens( Nbasis, Nbasis )
   real*8,  intent(in) :: time_step

   if ( allocated(RhoSaveA) ) deallocate(RhoSaveA)
   if ( allocated(RhoSaveB) ) deallocate(RhoSaveB)
   allocate( RhoSaveA(Nbasis,Nbasis) )
   allocate( RhoSaveB(Nbasis,Nbasis) )
   RhoSaveA = DCMPLX( real_dens )
   RhoSaveB = DCMPLX( real_dens )

   dt_nucl = time_step
   dt_elec = dt_nucl / DBLE(edyn_steps)

end subroutine ehrensetup
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
