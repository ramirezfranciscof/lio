!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
!  This subroutine calculates the part of the fock matrix that depends on
!  the electronic density, but it must receive as input the core fock.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmCheckNaNs( descript )
!------------------------------------------------------------------------------!
   use garcha_mod  , only: M, RMM, rhoalpha, rhobeta, Dbug, OPEN
   implicit none
   character(len=*) :: descript
   integer          :: MM, Dens_ti, Dens_tf, Fock_bi, Fock_bf, Fock_ai, Fock_af

   if (.not.Dbug) return

   MM =  M * (M+1)/2
   Dens_ti = 1
   Dens_tf = Dens_ti - 1 + MM
   Fock_bi = Dens_tf + 1
   Fock_bf = Fock_bi - 1 + MM
   Fock_ai = Fock_bf + 1
   Fock_af = Fock_ai - 1 + MM

   if (OPEN) then
      call SEEK_NaN( rhoalpha, 1      ,  MM    , "RHO_a " //trim(descript) )
      call SEEK_NaN( rhobeta , 1      ,  MM    , "RHO_b " //trim(descript) )
      call SEEK_NaN( RMM     , Dens_ti, Dens_tf, "RHO_t " //trim(descript) )
      call SEEK_NaN( RMM     , Fock_ai, Fock_af, "FOCK_a "//trim(descript) )
      call SEEK_NaN( RMM     , Fock_bi, Fock_bf, "FOCK_b "//trim(descript) )

   else
      call SEEK_NaN( RMM     , Dens_ti, Dens_tf, "RHO "   //trim(descript) )
      call SEEK_NaN( RMM     , Fock_ai, Fock_af, "FOCK "  //trim(descript) )
      
   end if

end subroutine rmmCheckNaNs
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
