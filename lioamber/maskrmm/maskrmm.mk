######################################################################
# EXTERNAL MAKE VARIABLES
objects   += maskrmm.o
src_paths += maskrmm

# INTERNAL DEPENDENCIES
internal_files := maskrmm.mk
internal_files += maskrmm.f90
internal_files += rmmput_dens.f90  rmmget_dens.f90
internal_files += rmmput_fock.f90  rmmget_fock.f90
internal_files += rmmcalc_focknucl.f90
internal_files += rmmcalc_fockdens.f90
$(obj_path)/maskrmm.o : $(internal_files)

# EXTERNAL USERS
object_users := ehrenfest.o
object_users := SCF.o
$(object_users:%.o=$(obj_path)/%.o) : $(obj_path)/maskrmm.mod
######################################################################
