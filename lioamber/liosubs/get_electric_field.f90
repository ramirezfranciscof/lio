!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  subroutine get_electric_field( time, field )
!
! This subroutine calculates the electric field at the time given
! according to a gaussian shape. The center, height and width of
! said gaussian are parameters that are a mixture between hardcoded
! and keyword inputs (the field max value and the time step are
! taken from the garcha_mod).
!
! Ideas for improvement:
!
!  o  Implement a general purpose (independent from any module)
!     gaussian shaper that takes as inputs the parameters of the
!     shape, and have this subroutine call the shaper with the
!     specifics of the shape (thus still hiding the parameters
!     derived from keywords from the caller of this subroutine).
!
!  o  Add the posibility of shaping the electric field in other
!     ways and not just a gaussian bell.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!
! Variables Declaration
!--------------------------------------------------------------------!
  use garcha_mod, only: Fx, Fy, Fz, tdstep
  implicit none
  real*8,intent(in)  :: time
  real*8,intent(out) :: field(3)

  real*8 :: sigma, time0, field0(3)
  real*8 :: factor, gexp

!
! Setting Parameters
!--------------------------------------------------------------------!
  time0=50*tdstep
  sigma=sqrt(5*tdstep)
  field0(1)=Fx
  field0(2)=Fy
  field0(3)=Fz

!
! Calculations
!--------------------------------------------------------------------!
  gexp=(time-time0)/sigma
  gexp=gexp*gexp
  factor=exp(-gexp)
  field=field0*factor

!
! Exit
!--------------------------------------------------------------------!
  return; end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
