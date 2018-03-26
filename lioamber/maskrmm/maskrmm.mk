######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += rmmput_dens.f90
INCLUDES += rmmget_dens.f90
INCLUDES += rmmput_densos.f90
INCLUDES += rmmput_fock.f90
INCLUDES += rmmget_fock.f90
INCLUDES += rmmput_core.f90
INCLUDES += rmmget_core.f90

INCLUDES += rmmCheck_NaNs.f90
INCLUDES += rmmCalc_init.f90
INCLUDES += rmmCalc_core.f90
INCLUDES += rmmCalc_ener.f90
INCLUDES += rmmCalc_fock.f90

INCLUDES += rmmcalc0_init.f90
INCLUDES += rmmcalc1_overlap.f90
INCLUDES += rmmcalc2_focknuc.f90
INCLUDES += rmmcalc3_fockele.f90

$(OBJPATH)/maskrmm.o : $(INCLUDES) maskrmm.mk
######################################################################
