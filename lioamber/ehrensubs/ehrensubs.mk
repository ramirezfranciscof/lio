######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += fzaDS2.f90
INCLUDES += setim.f90
INCLUDES += ehren_magnus.f90
INCLUDES += ehren_verlet.f90
INCLUDES += calc_forceDS.f90
INCLUDES += calc_forceDS_dss.f90
INCLUDES += calc_forceDS_dds.f90
INCLUDES += ehren_masses.f90

$(obj_path)/ehrensubs.o : $(INCLUDES) ehrensubs.mk
######################################################################
