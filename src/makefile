##### Specify where the code will be run #####


######## Choose compilers and options ########

F90         = mpif90
CC          = mpicc
#OPTFLAGS    = -O3 
#OPTFLAGS2   = -O3
#OPTFLAGS    = -fast -g -traceback
#OPTFLAGS2   = -fast -g -traceback
   OPTFLAGS    = -g -ffree-form -ffree-line-length-none -finteger-4-integer-8
   OPTFLAGS2   = -g -ffree-form -ffree-line-length-none -finteger-4-integer-8
   LIBS        = -lmpi

############## Define dependecies ############

EFILE        = 3dh
OBJS         =  parameter_mod.of90 \
		mesh_class.of90 \
		mesh2d.of90 \
		functions.of90 \
		inject.of90 \
		dipole_field.of90 \
		dipole_field_2d.of90 \
		wrtfile.of90 \
		wrtfile_Non_MPIO.of90 \
		3dHybrid.of90 \
		get_environment_variable1.ogcc \
		get_environment_variable2.ogcc \
		integer_to_character.ogcc \
		field_2d.of90 \
		etacal.of90 \
		etacal_2d.of90
$(EFILE): $(OBJS)
	@echo "linking..."
	$(F90) $(OPTFLAGS2) -o $(EFILE) $(OBJS) $(LIBS)

############## Compile C files ###############

%.ogcc: %.c makefile 
	$(CC) -o $@ -c $<

%.c: %.c.nashi makefile 
	cp -f $< $@

############## Compile F90 files #############

%.of90: %.f90 makefile parameter_mod.f90
	$(F90)  $(OPTFLAGS) -o $@ -c $< 

clean:
	-rm  *.mod  3dh *.of90 *.ogcc
cleanall:
	-rm  restart_data.*

