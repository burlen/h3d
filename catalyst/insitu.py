
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 4.3.1-447-g2592378 64 bits


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 4.3.1-447-g2592378

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1312, 803]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.StereoType = 0
      renderView1.CameraPosition = [-2.03446004499998, 3.7509065435536666, -76.0123093466646]
      renderView1.CameraFocalPoint = [15.500000000000009, 0.4999999999999981, -3.703117506228924e-15]
      renderView1.CameraViewUp = [0.013451601339380045, 0.9436796555336437, -0.33058699634636274]
      renderView1.CameraParallelScale = 20.442252632172327
      renderView1.Background = [0.32, 0.34, 0.43]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='den_%t.png', freq=1, fittoscreen=0, magnification=1, width=1312, height=803, cinema={})

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'H3d'
      # create a producer from a simulation input
      h3dRawData = coprocessor.CreateProducer(datadescription, 'input')

      # convert b, e, vi to vectors
      # create a new 'Calculator'
      calculator1 = Calculator(Input=h3dRawData)
      calculator1.ResultArrayName = 'b'
      calculator1.Function = 'bx*iHat + by*jHat + bz*kHat'

      # create a new 'Calculator'
      calculator2 = Calculator(Input=calculator1)
      calculator2.ResultArrayName = 'e'
      calculator2.Function = 'ex*iHat + ey*jHat + ez*kHat'

      # create a new 'Calculator'
      calculator3 = Calculator(Input=calculator2)
      calculator3.ResultArrayName = 'vi'
      calculator3.Function = 'vix*iHat + viy*jHat + viz*kHat'

      h3dData = calculator3


      # create a new 'Parallel Image Data Writer'
      parallelImageDataWriter1 = servermanager.writers.XMLPImageDataWriter(Input=h3dData)

      # register the writer with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the data, etc.
      coprocessor.RegisterWriter(parallelImageDataWriter1, filename='h3d_%t.pvti', freq=1)


      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'den'
      denLUT = GetColorTransferFunction('den')
      denLUT.RGBPoints = [0.8444057760083787, 0.231373, 0.298039, 0.752941, 
                          1.3419032719439268, 0.865003, 0.865003, 0.865003, 
                          1.839400767879474, 0.705882, 0.0156863, 0.14902]
      denLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'den'
      denPWF = GetOpacityTransferFunction('den')
      denPWF.Points = [0.8444057760083787, 0.0, 0.5, 0.0, 
                       1.839400767879474, 1.0, 0.5, 0.0]
      denPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from h3d
      h3dDisplay = Show(h3dData, renderView1)
      # trace defaults for the display properties.
      h3dDisplay.ColorArrayName = ['POINTS', 'den']
      h3dDisplay.SetRepresentationType('Surface With Edges')
      h3dDisplay.LookupTable = denLUT

      # show color legend
      h3dDisplay.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for denLUT in view renderView1
      denLUTColorBar = GetScalarBar(denLUT, renderView1)
      denLUTColorBar.Title = 'den'
      denLUTColorBar.ComponentTitle = ''

      # reset view to fit data
      renderView1.ResetCamera()
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(True, 1)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
