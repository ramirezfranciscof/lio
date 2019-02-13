!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module typedef_liomat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   use liokinds , only: wpk => liomat_wpk
   use error_log, only: check_stat

   type liomat

      private
      integer :: my_size = 0

      real(wpk), allocatable :: rvals
      real(wpk), allocatable :: ivals

      contains

!     SETUPS
      procedure, pass, public  :: unset_0
      procedure, pass, public  :: reset_0
      procedure, pass, public  :: setup_r
      procedure, pass, public  :: setup_i
      procedure, pass, public  :: setup_c

!     OPERATIONS
      procedure, pass, public  :: dot_scalar
      procedure, pass, public  :: add_direct
      procedure, pass, public  :: add_matmul
      procedure, pass, public  :: add_commut

   end type liomat

!------------------------------------------------------------------------------!
contains

#  include "resets.f90"
#  include "setups.f90"
#  include "dot_scalar.f90"
#  include "add_direct.f90"
#  include "add_matmul.f90"
#  include "add_commut.f90"

end module typedef_liomat
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
