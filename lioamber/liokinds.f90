!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liokinds
!------------------------------------------------------------------------------!
!
!  This module contains all kinds declarations used to indicate precisions
!  inside the fortran subrotuines of the module.
!  
!------------------------------------------------------------------------------!
   implicit none

!  BASIC DEFINITIONS OF KINDS
   integer, parameter :: spk = selected_real_kind(  6,  37 )
   integer, parameter :: dpk = selected_real_kind( 15, 307 )

!  WORKING PRECISION FOR OTHER MODULES
   integer, parameter :: liomat_wpk = dpk

end module liokinds
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
