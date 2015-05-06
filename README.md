repo for h3d globl hybrid magnetospheric code
There are 3 branches:

1. catalyst, has the catalyst bindings from our XSEDE paper
2. newcode, has bugfixes and updates post catalyst bindings
3. newcode-catalyst has the Catalyst binding for the newcode with
   ParaView 4.3
$SOURCE_DIR is the source directory for h3d.


#Build the catalyst adapter

    cd $SOURCE_DIR
    mkdir catalyst-build
    cd catalyst-build
    ccmake ../catalyst
    
Make sure you specify a valid ParaView build directory

    make
    cd ..

#Build the simulation

    cd src

Change PARAVIEW\_BUILD\_DIR variable in the makefile
to point the the same build directory for ParaView

    make
    cd ..

#Run the simulation that includes ParaView Catalyst visualization

    cd $SOURCE_DIR/catalyst-build
    ./run.sh
    
You can watch the progress of the simulation using

    tail -f 3dhout.log 

The current ParaView Catalyst script generates visualizations for:

1. a slice perpendicular on Z through the middle of the Z range colored by 'den'
2. streamlines for 'vi' colored by 'den'
3. a contour for 'den'
4. Enables connection to ParaView Live which can be used to modify the
   existing visualization pipeline. For more information about Catalyst Live
   see the ParaView Catalyst User Guide.


#Create a new visualization script for ParaView Catalyst.

To generate a new visualization script a user can edit the script or
use the Catalyst Script Generator plugin in ParaView. Execute the
following steps to use the Catalyst Script Generator plugin.

1. Start ParaView
2. Load the CatalystScriptGeneratorPlugin: select Tools / Manage Plugins...,
   select CatalystScriptGeneratorPlugin, click Load Selected, click Close.
3. Load a down-sampled version of the data produced by the simulation
4. Create desired visualization pipeline
5. Export the visualization pipeline into a Catalyst script: select
   CoProcesing / Export State, select desired source, select Live Visualization
   and desired views to be saved as images.
6. Copy the generated <new_script>.py into the catalyst folder
```cp <new_script>.py $SOURCE_DIR/catalyst```
7. Modify ```$SOURCE_DIR/catalyst/CMakeLists.txt``` such that the new script
   is used for visualization.
   ```configure_file(<new_script>.py insitu.py COPYONLY)```
   and reconfigure the run:
   ```cd $SOURCE_DIR/catalyst-build;cmake .```
   Now you can rerun the simulation and the visualization will be generated
   by the new  script.

For more information on how to create a visualization see the ParaView
User Guide.  For more information on how to generate a script for
ParaView Catalyst see the ParaView Catalyst User Guide.
