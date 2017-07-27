!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrendata
!------------------------------------------------------------------------------!
   implicit none
   real*8  :: stored_energy = 0.0d0
   integer :: step_number  = 0
!   integer :: last_step    = 1+120
   integer :: rstinp_funit  = 654321
   integer :: rstout_funit  = 123456
!   logical :: restart_dyn  = .false.

   complex*16,allocatable,dimension(:,:) :: RhoSaveA
   complex*16,allocatable,dimension(:,:) :: RhoSaveB

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
