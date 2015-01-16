
try: paraview.simple
except: from paraview.simple import *

# Make sure the helper methods are loaded
try:
  __cp_helper_script_loaded__
except:
  import vtkCoProcessorPython
  exec vtkCoProcessorPython.vtkCPHelperScripts.GetPythonHelperScript()

# Global variables that will hold the pipeline for each timestep
pipeline = None

# Live visualization inside ParaView
live_visu_active = False
pv_host = "localhost"
pv_port = 22222

write_frequencies    = {'input': [1]}
simulation_input_map = {'uh3d': 'input'}

# ----------------------- Pipeline definition -----------------------

def CreatePipeline(datadescription):
  class Pipeline:
    global cp_views, cp_writers
    RenderView1 = CreateCPView( CreateRenderView, "den_%t.png", 1, 0, 1, 838, 510, cp_views )
    RenderView1.LightSpecularColor = [1.0, 1.0, 1.0]
    RenderView1.UseOutlineForLODRendering = 0
    RenderView1.KeyLightAzimuth = 10.0
    RenderView1.UseTexturedBackground = 0
    RenderView1.UseLight = 1
    RenderView1.CameraPosition = [31.5, 31.5, 172.11920087683927]
    RenderView1.FillLightKFRatio = 3.0
    RenderView1.Background2 = [0.0, 0.0, 0.165]
    RenderView1.FillLightAzimuth = -10.0
    RenderView1.LODResolution = 0.5
    RenderView1.BackgroundTexture = []
    RenderView1.InteractionMode = '2D'
    RenderView1.StencilCapable = 1
    RenderView1.LightIntensity = 1.0
    RenderView1.CameraFocalPoint = [31.5, 31.5, 0.0]
    RenderView1.ImageReductionFactor = 2
    RenderView1.CameraViewAngle = 30.0
    RenderView1.CameraParallelScale = 44.54772721475249
    RenderView1.EyeAngle = 2.0
    RenderView1.HeadLightKHRatio = 3.0
    RenderView1.StereoRender = 0
    RenderView1.KeyLightIntensity = 0.75
    RenderView1.BackLightAzimuth = 110.0
    RenderView1.OrientationAxesInteractivity = 0
    RenderView1.UseInteractiveRenderingForSceenshots = 0
    RenderView1.UseOffscreenRendering = 0
    RenderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
    RenderView1.UseOffscreenRenderingForScreenshots = 0
    RenderView1.NonInteractiveRenderDelay = 0.0
    RenderView1.CenterOfRotation = [31.5, 31.5, 0.0]
    RenderView1.CameraParallelProjection = 0
    RenderView1.CompressorConfig = 'vtkSquirtCompressor 0 3'
    RenderView1.HeadLightWarmth = 0.5
    RenderView1.MaximumNumberOfPeels = 4
    RenderView1.LightDiffuseColor = [1.0, 1.0, 1.0]
    RenderView1.StereoType = 'Red-Blue'
    RenderView1.DepthPeeling = 1
    RenderView1.BackLightKBRatio = 3.5
    RenderView1.StereoCapableWindow = 1
    RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
    RenderView1.LightType = 'HeadLight'
    RenderView1.LightAmbientColor = [1.0, 1.0, 1.0]
    RenderView1.RemoteRenderThreshold = 20.0
    RenderView1.CacheKey = 0.0
    RenderView1.UseCache = 0
    RenderView1.KeyLightElevation = 50.0
    RenderView1.CenterAxesVisibility = 1
    RenderView1.MaintainLuminance = 0
    RenderView1.StillRenderImageReductionFactor = 1
    RenderView1.BackLightWarmth = 0.5
    RenderView1.FillLightElevation = -75.0
    RenderView1.MultiSamples = 0
    RenderView1.FillLightWarmth = 0.4
    RenderView1.AlphaBitPlanes = 1
    RenderView1.LightSwitch = 0
    RenderView1.OrientationAxesVisibility = 1
    RenderView1.CameraClippingRange = [170.3980088680709, 174.70098888999183]
    RenderView1.BackLightElevation = 0.0
    RenderView1.ViewTime = 0.0
    RenderView1.OrientationAxesOutlineColor = [1.0, 1.0, 1.0]
    RenderView1.LODThreshold = 5.0
    RenderView1.CollectGeometryThreshold = 100.0
    RenderView1.UseGradientBackground = 0
    RenderView1.KeyLightWarmth = 0.6
    RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]
    
    uh3d = CreateProducer( datadescription, "input" )
    
    a1_den_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
    
    a1_den_PVLookupTable = GetLookupTableForArray( "den", 1, Discretize=1, RGBPoints=[0.0, 0.23, 0.299, 0.754, 5.565536022186279, 0.706, 0.016, 0.15], UseLogScale=0, VectorComponent=0, NanColor=[0.25, 0.0, 0.0], NumberOfTableValues=256, EnableOpacityMapping=0, ColorSpace='Diverging', IndexedLookup=0, VectorMode='Magnitude', ScalarOpacityFunction=a1_den_PiecewiseFunction, HSVWrap=0, ScalarRangeInitialized=1.0, AllowDuplicateScalars=1, Annotations=[], LockScalarRange=0 )
    
    ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='den', Position2=[0.13, 0.5], TitleOpacity=1.0, TitleShadow=0, AutomaticLabelFormat=1, TitleFontSize=12, TitleColor=[1.0, 1.0, 1.0], AspectRatio=20.0, NumberOfLabels=5, ComponentTitle='', Resizable=1, TitleFontFamily='Arial', Visibility=0, LabelFontSize=12, LabelFontFamily='Arial', TitleItalic=0, Selectable=0, LabelItalic=0, Enabled=0, LabelColor=[1.0, 1.0, 1.0], Position=[0.87, 0.25], LabelBold=0, UseNonCompositedRenderer=1, LabelOpacity=1.0, TitleBold=0, LabelFormat='%-#6.3g', Orientation='Vertical', LabelShadow=0, LookupTable=a1_den_PVLookupTable, Repositionable=1 )
    GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)
    
    DataRepresentation1 = Show()
    DataRepresentation1.CubeAxesZAxisVisibility = 1
    DataRepresentation1.SelectionPointLabelColor = [0.5, 0.5, 0.5]
    DataRepresentation1.SelectionPointFieldDataArrayName = 'den'
    DataRepresentation1.SuppressLOD = 0
    DataRepresentation1.CubeAxesXGridLines = 0
    DataRepresentation1.CubeAxesYAxisTickVisibility = 1
    DataRepresentation1.Position = [0.0, 0.0, 0.0]
    DataRepresentation1.BackfaceRepresentation = 'Follow Frontface'
    DataRepresentation1.SelectionOpacity = 1.0
    DataRepresentation1.SelectionPointLabelShadow = 0
    DataRepresentation1.CubeAxesYGridLines = 0
    DataRepresentation1.CubeAxesZAxisTickVisibility = 1
    DataRepresentation1.OrientationMode = 'Direction'
    DataRepresentation1.Source.TipResolution = 6
    DataRepresentation1.ScaleMode = 'No Data Scaling Off'
    DataRepresentation1.Diffuse = 1.0
    DataRepresentation1.SelectionUseOutline = 0
    DataRepresentation1.CubeAxesZTitle = 'Z-Axis'
    DataRepresentation1.Specular = 0.1
    DataRepresentation1.SelectionVisibility = 1
    DataRepresentation1.InterpolateScalarsBeforeMapping = 1
    DataRepresentation1.CustomRangeActive = [0, 0, 0]
    DataRepresentation1.Origin = [0.0, 0.0, 0.0]
    DataRepresentation1.Source.TipLength = 0.35
    DataRepresentation1.CubeAxesVisibility = 0
    DataRepresentation1.Scale = [1.0, 1.0, 1.0]
    DataRepresentation1.SelectionCellLabelJustification = 'Left'
    DataRepresentation1.DiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation1.SelectionCellLabelOpacity = 1.0
    DataRepresentation1.CubeAxesInertia = 1
    DataRepresentation1.Source = "Arrow"
    DataRepresentation1.Source.Invert = 0
    DataRepresentation1.Masking = 0
    DataRepresentation1.Opacity = 1.0
    DataRepresentation1.LineWidth = 1.0
    DataRepresentation1.MeshVisibility = 0
    DataRepresentation1.Visibility = 1
    DataRepresentation1.SelectionCellLabelFontSize = 18
    DataRepresentation1.CubeAxesCornerOffset = 0.0
    DataRepresentation1.SelectionPointLabelJustification = 'Left'
    DataRepresentation1.OriginalBoundsRangeActive = [0, 0, 0]
    DataRepresentation1.SelectionPointLabelVisibility = 0
    DataRepresentation1.SelectOrientationVectors = ''
    DataRepresentation1.CubeAxesTickLocation = 'Inside'
    DataRepresentation1.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
    DataRepresentation1.CubeAxesYAxisVisibility = 1
    DataRepresentation1.SelectionPointLabelFontFamily = 'Arial'
    DataRepresentation1.Source.ShaftResolution = 6
    DataRepresentation1.CubeAxesUseDefaultYTitle = 1
    DataRepresentation1.SelectScaleArray = ''
    DataRepresentation1.CubeAxesYTitle = 'Y-Axis'
    DataRepresentation1.ColorAttributeType = 'POINT_DATA'
    DataRepresentation1.AxesOrigin = [0.0, 0.0, 0.0]
    DataRepresentation1.UserTransform = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    DataRepresentation1.SpecularPower = 100.0
    DataRepresentation1.Texture = []
    DataRepresentation1.SelectionCellLabelShadow = 0
    DataRepresentation1.AmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation1.MapScalars = 1
    DataRepresentation1.PointSize = 2.0
    DataRepresentation1.CubeAxesUseDefaultXTitle = 1
    DataRepresentation1.SelectionCellLabelFormat = ''
    DataRepresentation1.Scaling = 0
    DataRepresentation1.StaticMode = 0
    DataRepresentation1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
    DataRepresentation1.Source.TipRadius = 0.1
    DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
    DataRepresentation1.CubeAxesXAxisTickVisibility = 1
    DataRepresentation1.SelectionCellLabelVisibility = 0
    DataRepresentation1.NonlinearSubdivisionLevel = 1
    DataRepresentation1.CubeAxesColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Representation = 'Surface'
    DataRepresentation1.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation1.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
    DataRepresentation1.Orientation = [0.0, 0.0, 0.0]
    DataRepresentation1.CubeAxesXTitle = 'X-Axis'
    DataRepresentation1.BackfaceOpacity = 1.0
    DataRepresentation1.SelectionPointLabelFontSize = 18
    DataRepresentation1.SelectionCellFieldDataArrayName = 'vtkOriginalCellIds'
    DataRepresentation1.SelectionColor = [1.0, 0.0, 1.0]
    DataRepresentation1.Ambient = 0.0
    DataRepresentation1.CubeAxesXAxisMinorTickVisibility = 1
    DataRepresentation1.ScaleFactor = 6.300000000000001
    DataRepresentation1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
    DataRepresentation1.Source.ShaftRadius = 0.03
    DataRepresentation1.SelectMaskArray = ''
    DataRepresentation1.SelectionLineWidth = 2.0
    DataRepresentation1.CubeAxesZAxisMinorTickVisibility = 1
    DataRepresentation1.CubeAxesXAxisVisibility = 1
    DataRepresentation1.Interpolation = 'Gouraud'
    DataRepresentation1.SelectionCellLabelFontFamily = 'Arial'
    DataRepresentation1.SelectionCellLabelItalic = 0
    DataRepresentation1.CubeAxesYAxisMinorTickVisibility = 1
    DataRepresentation1.CubeAxesZGridLines = 0
    DataRepresentation1.SelectionPointLabelFormat = ''
    DataRepresentation1.SelectionPointLabelOpacity = 1.0
    DataRepresentation1.UseAxesOrigin = 0
    DataRepresentation1.CubeAxesFlyMode = 'Closest Triad'
    DataRepresentation1.Pickable = 1
    DataRepresentation1.CustomBoundsActive = [0, 0, 0]
    DataRepresentation1.CubeAxesGridLineLocation = 'All Faces'
    DataRepresentation1.SelectionRepresentation = 'Wireframe'
    DataRepresentation1.SelectionPointLabelBold = 0
    DataRepresentation1.ColorArrayName = 'den'
    DataRepresentation1.SelectionPointLabelItalic = 0
    DataRepresentation1.AllowSpecularHighlightingWithScalarColoring = 0
    DataRepresentation1.SpecularColor = [1.0, 1.0, 1.0]
    DataRepresentation1.CubeAxesUseDefaultZTitle = 1
    DataRepresentation1.LookupTable = a1_den_PVLookupTable
    DataRepresentation1.SelectionPointSize = 5.0
    DataRepresentation1.SelectionCellLabelBold = 0
    DataRepresentation1.Orient = 0
    
  return Pipeline()


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    if datadescription.GetForceOutput() == True:
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    for input_name in simulation_input_map.values():
       LoadRequestedData(datadescription, input_name)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global pipeline, cp_writers, cp_views
    timestep = datadescription.GetTimeStep()

    # Load the Pipeline if not created yet
    if not pipeline:
       pipeline = CreatePipeline(datadescription)
    else:
      # update to the new input and time
      UpdateProducers(datadescription)

    # Write output data
    WriteAllData(datadescription, cp_writers, timestep);

    # Write image capture (Last arg: rescale lookup table)
    WriteAllImages(datadescription, cp_views, timestep, True)

    # Live Visualization
    if live_visu_active:
       DoLiveInsitu(timestep, pv_host, pv_port)

