######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += ehrendyn.f90
INCLUDES += setim.f90
INCLUDES += calc_forceDS.f90
INCLUDES += calc_forceDS_dss.f90
INCLUDES += calc_forceDS_dds.f90

INCLUDES += ehren_magnus.f90
INCLUDES += ehren_verlet_e.f90
#INCLUDES += nuclear_verlet.f90
INCLUDES += ehren_masses.f90
INCLUDES += calc_kenergy.f90

INCLUDES += RMMcalc0_Init.f90
INCLUDES += RMMcalc1_Overlap.f90
INCLUDES += RMMcalc2_FockMao.f90


$(obj_path)/ehrenfest.o : $(INCLUDES) ehrenfest.mk
######################################################################
