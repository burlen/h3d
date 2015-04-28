# Where to find ParaView Catalyst Libraries
PARAVIEW_BUILD_DIR = /home/danlipsa/build/ParaView-ensight-crash

##### Specify where the code will be run #####


######## Choose compilers and options ########

F90         = mpif90
CC          = mpicc
OPTFLAGS    = -O3 -ffree-form -ffree-line-length-none
OPTFLAGS2   = -O3 -ffree-form -ffree-line-length-none
LIBS        = -lmpi -lUH3DAdaptor -lstdc++
VTKLIBS     = -lvtkPVCatalyst-pv4.3 -lvtkPVPythonCatalyst-pv4.3 -lvtkPVServerManagerApplication-pv4.3 -lvtkPVServerManagerCore-pv4.3 -lvtkPVServerImplementationCore-pv4.3 -lvtkPVClientServerCoreCore-pv4.3 -lvtkFiltersParallel-pv4.3 -lvtkFiltersModeling-pv4.3 -lvtkRenderingCore-pv4.3 -lvtkFiltersExtraction-pv4.3 -lvtkFiltersStatistics-pv4.3 -lvtkImagingFourier-pv4.3 -lvtkImagingCore-pv4.3 -lvtkFiltersGeometry-pv4.3 -lvtkFiltersSources-pv4.3 -lvtkFiltersGeneral-pv4.3 -lvtkCommonComputationalGeometry-pv4.3 -lvtkCommonColor-pv4.3 -lvtkFiltersProgrammable-pv4.3 -lvtkPVVTKExtensionsCore-pv4.3 -lvtkPVCommon-pv4.3 -lvtkClientServer-pv4.3 -lvtkIOXMLParser-pv4.3 -lvtkFiltersCore-pv4.3 -lvtkParallelMPI-pv4.3 -lvtkParallelCore-pv4.3 -lvtkIOLegacy-pv4.3 -lvtkIOCore-pv4.3 -lvtkCommonExecutionModel-pv4.3 -lvtkCommonDataModel-pv4.3 -lvtkCommonSystem-pv4.3 -lvtkCommonTransforms-pv4.3 -lvtkCommonMisc-pv4.3 -lvtkCommonMath-pv4.3 -lvtkPythonInterpreter-pv4.3 -lvtkWrappingPython27Core-pv4.3 -lvtkCommonCore-pv4.3 -lvtksys-pv4.3
LIBFLAGS    = -L../catalyst-build/ -L$(PARAVIEW_BUILD_DIR)/lib

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
	$(F90) $(OPTFLAGS2) -o $(EFILE) $(OBJS) $(LIBFLAGS) $(LIBS) $(VTKLIBS)

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
   
