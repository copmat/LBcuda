# state file generated using paraview version 5.7.0

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1481, 808]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [33.0, 20.000000476837158, 96.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [37.25147045499608, 19.355838286686236, 142.90713232402527]
renderView1.CameraFocalPoint = [37.25147045499608, 19.355838286686236, 96.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 26.024152482097232
renderView1.Background = [0.0, 0.0, 0.0]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------
rangeINI = 0
rangeEND = 3001
rangeSTEP= 50

flvtk  = [ 'output/particle_{:06d}.vtk'.format(i) for i in range(rangeINI,rangeEND,rangeSTEP) ]
flisfl = [ 'output/isfluid_{:06d}.vti'.format(i) for i in range(rangeINI,rangeEND,rangeSTEP) ]


# create a new 'XML Image Data Reader'
isfluid_00000 = XMLImageDataReader(FileName=flisfl)
isfluid_00000.PointArrayStatus = ['isfluid']

# create a new 'Legacy VTK Reader'
particle_00000 = LegacyVTKReader(FileNames=flvtk)

# create a new 'Glyph'
glyph1 = Glyph(Input=particle_00000, GlyphType='Sphere')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 9.0
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'

# init the 'Sphere' selected for 'GlyphType'
glyph1.GlyphType.ThetaResolution = 16
glyph1.GlyphType.PhiResolution = 16

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)

# get color transfer function/color map for 'Type'
typeLUT = GetColorTransferFunction('Type')
typeLUT.RGBPoints = [1.0, 0.278431372549, 0.278431372549, 0.858823529412, 1.143, 0.0, 0.0, 0.360784313725, 1.285, 0.0, 1.0, 1.0, 1.429, 0.0, 0.501960784314, 0.0, 1.571, 1.0, 1.0, 0.0, 1.714, 1.0, 0.380392156863, 0.0, 1.857, 0.419607843137, 0.0, 0.0, 2.0, 0.878431372549, 0.301960784314, 0.301960784314]
typeLUT.ColorSpace = 'RGB'
typeLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'Type']
glyph1Display.LookupTable = typeLUT
glyph1Display.OSPRayScaleArray = 'Type'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 2.0774353027343753
glyph1Display.SelectScaleArray = 'Type'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Type'
glyph1Display.GaussianRadius = 0.10387176513671875
glyph1Display.SetScaleArray = ['POINTS', 'Type']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Type']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]

# show data from isfluid_00000
isfluid_00000Display = Show(isfluid_00000, renderView1)

# get color transfer function/color map for 'isfluid'
isfluidLUT = GetColorTransferFunction('isfluid')
isfluidLUT.RGBPoints = [1.0, 0.278431372549, 0.278431372549, 0.858823529412, 1.572, 0.0, 0.0, 0.360784313725, 2.1399999999999997, 0.0, 1.0, 1.0, 2.716, 0.0, 0.501960784314, 0.0, 3.284, 1.0, 1.0, 0.0, 3.856, 1.0, 0.380392156863, 0.0, 4.428, 0.419607843137, 0.0, 0.0, 5.0, 0.878431372549, 0.301960784314, 0.301960784314]
isfluidLUT.ColorSpace = 'RGB'
isfluidLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'isfluid'
isfluidPWF = GetOpacityTransferFunction('isfluid')
isfluidPWF.Points = [1.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]
isfluidPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
isfluid_00000Display.Representation = 'Slice'
isfluid_00000Display.ColorArrayName = ['POINTS', 'isfluid']
isfluid_00000Display.LookupTable = isfluidLUT
isfluid_00000Display.OSPRayScaleArray = 'isfluid'
isfluid_00000Display.OSPRayScaleFunction = 'PiecewiseFunction'
isfluid_00000Display.SelectOrientationVectors = 'None'
isfluid_00000Display.ScaleFactor = 12.700000000000001
isfluid_00000Display.SelectScaleArray = 'None'
isfluid_00000Display.GlyphType = 'Arrow'
isfluid_00000Display.GlyphTableIndexArray = 'None'
isfluid_00000Display.GaussianRadius = 0.635
isfluid_00000Display.SetScaleArray = ['POINTS', 'isfluid']
isfluid_00000Display.ScaleTransferFunction = 'PiecewiseFunction'
isfluid_00000Display.OpacityArray = ['POINTS', 'isfluid']
isfluid_00000Display.OpacityTransferFunction = 'PiecewiseFunction'
isfluid_00000Display.DataAxesGrid = 'GridAxesRepresentation'
isfluid_00000Display.PolarAxes = 'PolarAxesRepresentation'
isfluid_00000Display.ScalarOpacityUnitDistance = 1.732050807568877
isfluid_00000Display.ScalarOpacityFunction = isfluidPWF
isfluid_00000Display.Slice = 95

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
isfluid_00000Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isfluid_00000Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isfluid_00000Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Type'
typePWF = GetOpacityTransferFunction('Type')
typePWF.Points = [1.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
typePWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------
