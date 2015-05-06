from paraview.simple import *
from paraview import coprocessing
import vtk

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 4.3.1-480-g31dde1e 64 bits
# Generates a slice, volume rendering, contour of 'den' and streamlines of vi
# colored by 'den'


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 4.3.1-480-g31dde1e

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # # Create a new 'Render View'
      # renderView1 = CreateView('RenderView')
      # renderView1.ViewSize = [483, 387]
      # renderView1.AxesGrid = 'GridAxes3DActor'
      # renderView1.CenterOfRotation = [15.5, 15.5, 15.5]
      # renderView1.StereoType = 0
      # renderView1.CameraPosition = [-13.26567916460982, 59.02220746110771, 112.81887354587771]
      # renderView1.CameraFocalPoint = [15.5, 15.5, 15.5]
      # renderView1.CameraViewUp = [0.11789779112572521, 0.9190266663639157, -0.3761516946122972]
      # renderView1.CameraParallelScale = 28.578838324886476
      # renderView1.Background = [0.32, 0.34, 0.43]

      # # init the 'GridAxes3DActor' selected for 'AxesGrid'
      # renderView1.AxesGrid.Visibility = 1
      # renderView1.AxesGrid.AxesToLabel = 47

      # # register the view with coprocessor
      # # and provide it with information such as the filename to use,
      # # how frequently to write the images, etc.
      # coprocessor.RegisterView(renderView1,
      #     filename='volume_%t.png', freq=1, fittoscreen=0, magnification=1, width=483, height=387, cinema={})

      # Create a new 'Render View'
      renderView2 = CreateView('RenderView')
      renderView2.ViewSize = [483, 386]
      renderView2.AxesGrid = 'GridAxes3DActor'
      renderView2.CenterOfRotation = [15.5, 15.5, 15.5]
      renderView2.StereoType = 0
      renderView2.CameraPosition = [-13.26567916460982, 59.02220746110771, 112.81887354587771]
      renderView2.CameraFocalPoint = [15.5, 15.5, 15.5]
      renderView2.CameraViewUp = [0.11789779112572521, 0.9190266663639157, -0.3761516946122972]
      renderView2.CameraParallelScale = 28.578838324886476
      renderView2.Background = [0.32, 0.34, 0.43]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView2.AxesGrid.Visibility = 1
      renderView2.AxesGrid.AxesToLabel = 47

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView2,
          filename='contour_%t.png', freq=1, fittoscreen=0, magnification=1, width=483, height=386, cinema={})

      # Create a new 'Render View'
      renderView3 = CreateView('RenderView')
      renderView3.ViewSize = [483, 387]
      renderView3.InteractionMode = '2D'
      renderView3.AxesGrid = 'GridAxes3DActor'
      renderView3.CenterOfRotation = [15.5, 15.5, 15.5]
      renderView3.StereoType = 0
      renderView3.CameraPosition = [15.5, 14.897040729220771, 105.65767664977295]
      renderView3.CameraFocalPoint = [15.5, 14.897040729220771, 15.5]
      renderView3.CameraParallelScale = 23.33452377915607
      renderView3.Background = [0.32, 0.34, 0.43]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView3.AxesGrid.Visibility = 1
      renderView3.AxesGrid.AxesToLabel = 39

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView3,
          filename='slice_%t.png', freq=1, fittoscreen=0, magnification=1, width=483, height=387, cinema={})

      # Create a new 'Render View'
      renderView4 = CreateView('RenderView')
      renderView4.ViewSize = [483, 386]
      renderView4.AxesGrid = 'GridAxes3DActor'
      renderView4.CenterOfRotation = [15.5, 15.5, 15.5]
      renderView4.StereoType = 0
      renderView4.CameraPosition = [-13.26567916460982, 59.02220746110771, 112.81887354587771]
      renderView4.CameraFocalPoint = [15.5, 15.5, 15.5]
      renderView4.CameraViewUp = [0.11789779112572521, 0.9190266663639157, -0.3761516946122972]
      renderView4.CameraParallelScale = 28.578838324886476
      renderView4.Background = [0.32, 0.34, 0.43]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView4.AxesGrid.Visibility = 1
      renderView4.AxesGrid.AxesToLabel = 47

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView4,
          filename='streamlines_%t.png', freq=1, fittoscreen=0, magnification=1, width=483, height=386, cinema={})

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML Partitioned Image Data Reader'
      # create a producer from a simulation input
      h3d_ = coprocessor.CreateProducer(datadescription, 'input')

      # convert b, e, vi to vectors
      # create a new 'Calculator'
      calculator1 = Calculator(Input=h3d_)
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

      h3d = calculator3

      # create a new 'Text'
      text2 = Text()
      text2.Text = 'den'

      # create a new 'Annotate Time'
      annotateTime3 = AnnotateTime()
      annotateTime3.Format = 'Index: %0.0f'

      # create a new 'Text'
      text3 = Text()
      text3.Text = 'vi'

      # create a new 'Text'
      text4 = Text()
      text4.Text = 'den'

      # create a new 'Annotate Time'
      annotateTime4 = AnnotateTime()
      annotateTime4.Format = 'Index: %0.0f'

      # create a new 'Text'
      text1 = Text()
      text1.Text = 'den'

      # create a new 'Slice'
      slice3 = Slice(Input=h3d)
      slice3.SliceType = 'Plane'
      slice3.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice3.SliceType.Origin = [15.5, 15.5, 15.5]
      slice3.SliceType.Normal = [0.0, 0.0, 1.0]

      # create a new 'Contour'
      contour1 = Contour(Input=h3d)
      contour1.ContourBy = ['POINTS', 'den']
      contour1.ComputeScalars = 1
      contour1.Isosurfaces = [3.0]
      contour1.PointMergeMethod = 'Uniform Binning'

      # create a new 'Annotate Time'
      annotateTime2 = AnnotateTime()
      annotateTime2.Format = 'Index: %0.0f'

      # create a new 'Stream Tracer'
      streamTracer1 = StreamTracer(Input=h3d,
          SeedType='Point Source')
      streamTracer1.Vectors = ['POINTS', 'vi']
      streamTracer1.MaximumStreamlineLength = 33.0

      # init the 'Point Source' selected for 'SeedType'
      streamTracer1.SeedType.Center = [10.0, 15.5, 15.5]
      streamTracer1.SeedType.Radius = 4.0

      # create a new 'Tube'
      tube1 = Tube(Input=streamTracer1)
      tube1.Scalars = ['POINTS', 'AngularVelocity']
      tube1.Vectors = ['POINTS', 'Normals']
      tube1.Radius = 0.3299899423122406

      # create a new 'Annotate Time'
      annotateTime1 = AnnotateTime()
      annotateTime1.Format = 'Index: %0.0f'

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'den'
      denLUT = GetColorTransferFunction('den')
      denLUT.RGBPoints = [0.04999999999999998, 0.231373, 0.298039, 0.752941, 4.525, 0.865003, 0.865003, 0.865003, 9.0, 0.705882, 0.0156863, 0.14902]
      denLUT.LockScalarRange = 1
      denLUT.NanColor = [0.247059, 0.0, 0.0]
      denLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'den'
      denPWF = GetOpacityTransferFunction('den')
      denPWF.Points = [0.04999999999999998, 0.0, 0.5, 0.0, 4.01194270827862, 0.1184210553765297, 0.5, 0.0, 9.0, 1.0, 0.5, 0.0]
      denPWF.ScalarRangeInitialized = 1

      # # ----------------------------------------------------------------
      # # setup the visualization in view 'renderView1'
      # # ----------------------------------------------------------------

      # # show data from h3d
      # h3dDisplay = Show(h3d, renderView1)
      # # trace defaults for the display properties.
      # h3dDisplay.Representation = 'Volume'
      # h3dDisplay.ColorArrayName = ['POINTS', 'den']
      # h3dDisplay.LookupTable = denLUT
      # h3dDisplay.ScalarOpacityFunction = denPWF
      # h3dDisplay.ScalarOpacityUnitDistance = 1.7320508075688776
      # h3dDisplay.Slice = 16

      # # show color legend
      # h3dDisplay.SetScalarBarVisibility(renderView1, True)

      # # show data from text2
      # text2Display = Show(text2, renderView1)
      # # trace defaults for the display properties.
      # text2Display.Position = [0.2139, 0.9]

      # # show data from annotateTime2
      # annotateTime2Display = Show(annotateTime2, renderView1)
      # # trace defaults for the display properties.
      # annotateTime2Display.Position = [0.321784, 0.9]

      # # setup the color legend parameters for each legend in this view

      # # get color legend/bar for denLUT in view renderView1
      # denLUTColorBar = GetScalarBar(denLUT, renderView1)
      # denLUTColorBar.Position = [0.9036099585062238, 0.04145077720207259]
      # denLUTColorBar.Position2 = [0.12000000000000033, 0.4299999999999997]
      # denLUTColorBar.Title = 'den'
      # denLUTColorBar.ComponentTitle = ''
      # denLUTColorBar.TitleFontSize = 10
      # denLUTColorBar.LabelFontSize = 10
      # denLUTColorBar.AutomaticLabelFormat = 0
      # denLUTColorBar.LabelFormat = '%-6.1g'
      # denLUTColorBar.RangeLabelFormat = '%-4.3g'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView2'
      # ----------------------------------------------------------------

      # show data from contour1
      contour1Display = Show(contour1, renderView2)
      # trace defaults for the display properties.
      contour1Display.ColorArrayName = ['POINTS', 'den']
      contour1Display.DiffuseColor = [0.0, 0.65, 1.0]
      contour1Display.LookupTable = denLUT

      # show color legend
      contour1Display.SetScalarBarVisibility(renderView2, True)

      # show data from h3d
      h3dDisplay_1 = Show(h3d, renderView2)
      # trace defaults for the display properties.
      h3dDisplay_1.Representation = 'Outline'
      h3dDisplay_1.ColorArrayName = [None, '']
      h3dDisplay_1.ScalarOpacityUnitDistance = 1.7320508075688776
      h3dDisplay_1.Slice = 16

      # show data from text4
      text4Display = Show(text4, renderView2)
      # trace defaults for the display properties.
      text4Display.Position = [0.261618, 0.9]

      # show data from annotateTime4
      annotateTime4Display = Show(annotateTime4, renderView2)
      # trace defaults for the display properties.
      annotateTime4Display.Position = [0.38195, 0.9]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for denLUT in view renderView2
      denLUTColorBar_1 = GetScalarBar(denLUT, renderView2)
      denLUTColorBar_1.Position = [0.9139834024896266, 0.04675324675324674]
      denLUTColorBar_1.Position2 = [0.11999999999999955, 0.42999999999999994]
      denLUTColorBar_1.Title = 'den'
      denLUTColorBar_1.ComponentTitle = ''
      denLUTColorBar_1.TitleFontSize = 9
      denLUTColorBar_1.LabelFontSize = 9
      denLUTColorBar_1.AutomaticLabelFormat = 0
      denLUTColorBar_1.LabelFormat = '%-6.1g'
      denLUTColorBar_1.RangeLabelFormat = '%-6.1g'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView3'
      # ----------------------------------------------------------------

      # show data from slice3
      slice3Display = Show(slice3, renderView3)
      # trace defaults for the display properties.
      slice3Display.ColorArrayName = ['POINTS', 'den']
      slice3Display.DiffuseColor = [0.0, 0.65, 1.0]
      slice3Display.LookupTable = denLUT

      # show color legend
      slice3Display.SetScalarBarVisibility(renderView3, True)

      # show data from text1
      text1Display = Show(text1, renderView3)
      # trace defaults for the display properties.
      text1Display.Position = [0.269917, 0.9]

      # show data from annotateTime1
      annotateTime1Display = Show(annotateTime1, renderView3)
      # trace defaults for the display properties.
      annotateTime1Display.Position = [0.398548, 0.9]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for denLUT in view renderView3
      denLUTColorBar_2 = GetScalarBar(denLUT, renderView3)
      denLUTColorBar_2.Position = [0.8641908713692946, 0.05627123900094361]
      denLUTColorBar_2.Position2 = [0.11999999999999966, 0.42999999999999855]
      denLUTColorBar_2.Title = 'den'
      denLUTColorBar_2.ComponentTitle = ''
      denLUTColorBar_2.TitleFontSize = 10
      denLUTColorBar_2.LabelFontSize = 10
      denLUTColorBar_2.AutomaticLabelFormat = 0
      denLUTColorBar_2.LabelFormat = '%-6.1g'
      denLUTColorBar_2.RangeLabelFormat = '%-4.1g'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView4'
      # ----------------------------------------------------------------

      # show data from h3d
      h3dDisplay_2 = Show(h3d, renderView4)
      # trace defaults for the display properties.
      h3dDisplay_2.Representation = 'Outline'
      h3dDisplay_2.ColorArrayName = [None, '']
      h3dDisplay_2.ScalarOpacityUnitDistance = 1.7320508075688776
      h3dDisplay_2.Slice = 16

      # show data from tube1
      tube1Display = Show(tube1, renderView4)
      # trace defaults for the display properties.
      tube1Display.ColorArrayName = ['POINTS', 'den']
      tube1Display.DiffuseColor = [0.0, 0.65, 1.0]
      tube1Display.LookupTable = denLUT

      # show color legend
      tube1Display.SetScalarBarVisibility(renderView4, True)

      # show data from text3
      text3Display = Show(text3, renderView4)
      # trace defaults for the display properties.
      text3Display.Position = [0.311411, 0.9]

      # show data from annotateTime3
      annotateTime3Display = Show(annotateTime3, renderView4)
      # trace defaults for the display properties.
      annotateTime3Display.Position = [0.398548, 0.9]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for denLUT in view renderView4
      denLUTColorBar_3 = GetScalarBar(denLUT, renderView4)
      denLUTColorBar_3.Position = [0.8977178423236515, 0.047402597402597405]
      denLUTColorBar_3.Title = 'den'
      denLUTColorBar_3.ComponentTitle = ''
      denLUTColorBar_3.TitleFontSize = 9
      denLUTColorBar_3.LabelFontSize = 9
      denLUTColorBar_3.AutomaticLabelFormat = 0
      denLUTColorBar_3.LabelFormat = '%-6.0g'
      denLUTColorBar_3.RangeLabelFormat = '%-6.0g'
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
