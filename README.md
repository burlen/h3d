repo for h3d globl hybrid magnetospheric code
initially there are 2 branches:

1. catalyst, has the catalyst bindings from our XSEDE paper
2. newcode, has bugfixes and updates post catalyst bindings
3. newcode-catalyst has the catalyst binding for the new code

#Build the catalyst adapter

    mkdir catalyst-build
    cd catalyst-build
    ccmake ../catalyst
    
Make sure you specify a valid ParaView build directory

    make
    cd ..

#Build the simulation

    cd src

Change PARAVIEW\_BUILD\_DIR variable in the makefile

    make
    cd ..

#Run the simulation

    cd catalyst-build
    ./run.sh
    
You can watch the progress of the simulation using

    tail -f 3dhout.log 
