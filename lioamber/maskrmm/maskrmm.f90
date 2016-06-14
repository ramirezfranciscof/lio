!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module maskrmm
  implicit none
  private
  public :: rmmcalc_focknucl
  public :: rmmcalc_fockdens
!
!
!--------------------------------------------------------------------!
  public :: rmmget_dens
  interface rmmget_dens
    module procedure rmmget_dens_r
    module procedure rmmget_dens_d
    module procedure rmmget_dens_c
    module procedure rmmget_dens_z
  end interface

  public :: rmmput_dens
  interface rmmput_dens
    module procedure rmmput_dens_r
    module procedure rmmput_dens_d
    module procedure rmmput_dens_c
    module procedure rmmput_dens_z
  end interface

  public :: rmmget_fock
  interface rmmget_fock
    module procedure rmmget_fock_r
    module procedure rmmget_fock_d
  end interface

  public :: rmmput_fock
  interface rmmput_fock
    module procedure rmmput_fock_r
    module procedure rmmput_fock_d
  end interface
!
!
!--------------------------------------------------------------------!
contains
# include "rmmput_dens.f90"
# include "rmmget_dens.f90"
# include "rmmput_fock.f90"
# include "rmmget_fock.f90"
# include "rmmcalc_focknucl.f90"
# include "rmmcalc_fockdens.f90"
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
