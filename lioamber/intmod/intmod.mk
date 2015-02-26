######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += setim.f sumit.f
INCLUDES += fzaDS2.f

$(obj_path)/mathsubs.o : $(INCLUDES) intmod.mk
######################################################################
