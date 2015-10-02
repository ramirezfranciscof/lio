######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += fzaDS2.f90
INCLUDES += setim.f90

$(obj_path)/intmod.o : $(INCLUDES) intmod.mk
######################################################################
