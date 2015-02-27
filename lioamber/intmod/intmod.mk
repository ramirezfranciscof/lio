######################################################################
# INTERNAL DEPENDENCIES
INCLUDES :=
INCLUDES += fzaDS2.f
INCLUDES += setim.f sumit.f

$(obj_path)/intmod.o : $(INCLUDES) intmod.mk
######################################################################
