!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liounits
!------------------------------------------------------------------------------!
!
!  This module contains parameters for all changes of units the code might
!  need to do. It is highly advised that ALL variables that contains any
!  physical quantities be in atomic units: all parameters provided will be
!  to convert TO this units, or to convert from this units (with the ammount
!  of aus IN the other unit), by multiplying them. That is:
!
!  * To write some time in ps:
!
!  [TIME IN PS] = [TIME IN AU] * AU_TO_PICOSECOND
!
!  * To read some time in ps:
!
!  [TIME IN AU] = [TIME IN PS] * AU_IN_PICOSECOND
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use liokinds, only: dpk

   implicit none

!  TIME UNITS
   real(dpk), parameter :: AU_TO_PICOSECOND = 
   real(dpk), parameter :: AU_IN_PICOSECOND =

!  DISTANCE UNITS (BOHR)
   real(dpk), parameter :: AU_TO_ANGSTROM =
   real(dpk), parameter :: AU_IN_ANGSTROM =

!  ENERGY UNITS (HARTREE)
   real(dpk), parameter :: AU_TO_KCAL_MOL =
   real(dpk), parameter :: AU_TO_KCAL_MOL =

end module liounits
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
