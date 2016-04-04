!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  module ehrensubs
!--------------------------------------------------------------------!
!
!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
  implicit none
  include 'calc_Dmat_h.f90'
  contains
!
!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
  include 'setim.f90'
  include 'fzaDS2.f90'

  include 'calc_forceDS.f90'
  include 'calc_forceDS_dss.f90'
  include 'calc_forceDS_dds.f90'

  include 'calc_Dmat_cholesky.f90'

  include 'ehren_cholesky.f90'
  include 'ehren_magnus.f90'
  include 'ehren_verlet_e.f90'
  include 'ehren_verlet_n.f90'
  include 'ehren_masses.f90'
  include 'calc_kenergy.f90'

  end module
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
