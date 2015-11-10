######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += fzaDS2.f90
INCLUDES += setim.f90
INCLUDES += forceDS.f90

$(obj_path)/ehrendyn.o : $(INCLUDES) ehrendyn.mk
######################################################################
