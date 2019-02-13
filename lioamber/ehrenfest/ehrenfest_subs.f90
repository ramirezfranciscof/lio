!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrenfest_subs

   use liokinds , only: spk, dpk, wpk => liomat_wpk
   use error_log, only: check_stat

   implicit none
   contains
#  include "ehrenfest_main.f90"

end module ehrenfest_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
