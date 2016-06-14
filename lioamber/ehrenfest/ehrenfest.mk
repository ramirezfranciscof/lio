######################################################################
# EXTERNAL MAKE VARIABLES
objects   += ehrenfest.o
src_paths += ehrenfest


# INTERNAL DEPENDENCIES
internal_files := ehrenfest.mk
internal_files += ehrenfest.f90
internal_files += ehrenstep_verlet.f90
#internal_files += ehrenstep_magnus.f90

internal_files += setim.f90
internal_files += calc_forceDS.f90
internal_files += calc_forceDS_dss.f90
internal_files += calc_forceDS_dds.f90
internal_files += calc_Dmat.f90

internal_files += ehren_cholesky.f90
internal_files += ehren_dipole.f90
internal_files += ehren_magnus.f90
internal_files += ehren_verlet_e.f90
internal_files += ehren_masses.f90
internal_files += calc_kenergy.f90
internal_files += nuclear_verlet.f90

internal_files += RMMcalc0_Init.f90
internal_files += RMMcalc1_Overlap.f90
internal_files += RMMcalc2_FockMao.f90
$(obj_path)/ehrenfest.o : $(internal_files)


# EXTERNAL USERS
object_users := SCF_in.o ehrendyn.o ehrentest.o
$(object_users:%.o=$(obj_path)/%.o) : $(obj_path)/ehrenfest.mod
######################################################################
