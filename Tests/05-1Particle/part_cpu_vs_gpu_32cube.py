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
renderView1.ViewSize = [1474, 389]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [16.5, 16.5, 16.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [16.5, 16.5, 120.31100408081898]
renderView1.CameraFocalPoint = [16.5, 16.5, 16.5]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 26.86826494733145
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.ViewSize = [1474, 388]
renderView2.InteractionMode = '2D'
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.CenterOfRotation = [16.5, 16.5, 16.5]
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraPosition = [16.5, 16.5, 120.31100408081898]
renderView2.CameraFocalPoint = [16.5, 16.5, 16.5]
renderView2.CameraFocalDisk = 1.0
renderView2.CameraParallelScale = 26.86826494733145
renderView2.Background = [0.32, 0.34, 0.43]
renderView2.BackEnd = 'OSPRay raycaster'
renderView2.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitVertical(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------
rangeINI = 0
rangeEND = 5001
rangeSTEP= 20

flvtk = [ 'output/particle_{:06d}.vtk'.format(i) for i in range(rangeINI,rangeEND,rangeSTEP) ]
flrho = [ 'output/rho_{:06d}.vti'.format(i) for i in range(rangeINI,rangeEND,rangeSTEP) ]

flatm = [ 'cpu/outatm00{:06d}.vtk'.format(i) for i in range(rangeINI,rangeEND,rangeSTEP) ]
flrhocpu = [ 'cpu/out_rho1_00{:06d}.vti'.format(i) for i in range(rangeINI,rangeEND,rangeSTEP) ]

# create a new 'Legacy VTK Reader'
legacyVTKReader2 = LegacyVTKReader(FileNames=flatm)

# create a new 'Legacy VTK Reader'
legacyVTKReader1 = LegacyVTKReader(FileNames=flvtk)

# create a new 'XML Image Data Reader'
xMLImageDataReader2 = XMLImageDataReader(FileName=flrhocpu)
xMLImageDataReader2.PointArrayStatus = ['rho1    ']

# create a new 'XML Image Data Reader'
xMLImageDataReader1 = XMLImageDataReader(FileName=flrho)
xMLImageDataReader1.PointArrayStatus = ['rho']

# create a new 'Glyph'
glyph2 = Glyph(Input=legacyVTKReader2,
    GlyphType='Sphere')
glyph2.OrientationArray = ['POINTS', 'No orientation array']
glyph2.ScaleArray = ['POINTS', 'No scale array']
glyph2.ScaleFactor = 9.0
glyph2.GlyphTransform = 'Transform2'
glyph2.GlyphMode = 'All Points'

# create a new 'Glyph'
glyph1 = Glyph(Input=legacyVTKReader1,
    GlyphType='Sphere')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 9.0
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)

# get color transfer function/color map for 'Type'
typeLUT = GetColorTransferFunction('Type')
typeLUT.RGBPoints = [1.0, 0.278431372549, 0.278431372549, 0.858823529412, 1.000034912109375, 0.0, 0.0, 0.360784313725, 1.000069580078125, 0.0, 1.0, 1.0, 1.000104736328125, 0.0, 0.501960784314, 0.0, 1.000139404296875, 1.0, 1.0, 0.0, 1.00017431640625, 1.0, 0.380392156863, 0.0, 1.000209228515625, 0.419607843137, 0.0, 0.0, 1.000244140625, 0.878431372549, 0.301960784314, 0.301960784314]
typeLUT.ColorSpace = 'RGB'
typeLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'Type']
glyph1Display.LookupTable = typeLUT
glyph1Display.OSPRayScaleArray = 'Type'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 0.45
glyph1Display.SelectScaleArray = 'Type'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Type'
glyph1Display.GaussianRadius = 0.0225
glyph1Display.SetScaleArray = ['POINTS', 'Type']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Type']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# show data from xMLImageDataReader1
xMLImageDataReader1Display = Show(xMLImageDataReader1, renderView1)

# get color transfer function/color map for 'rho'
rhoLUT = GetColorTransferFunction('rho')
rhoLUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.1683459868431091, 0.0, 0.0, 0.360784313725, 0.3355147290229797, 0.0, 1.0, 1.0, 0.5050379605293274, 0.0, 0.501960784314, 0.0, 0.6722067027091979, 1.0, 1.0, 0.0, 0.8405526895523071, 1.0, 0.380392156863, 0.0, 1.0088986763954162, 0.419607843137, 0.0, 0.0, 1.1772446632385254, 0.878431372549, 0.301960784314, 0.301960784314]
rhoLUT.ColorSpace = 'RGB'
rhoLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'rho'
rhoPWF = GetOpacityTransferFunction('rho')
rhoPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.1772446632385254, 1.0, 0.5, 0.0]
rhoPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
xMLImageDataReader1Display.Representation = 'Volume'
xMLImageDataReader1Display.ColorArrayName = ['POINTS', 'rho']
xMLImageDataReader1Display.LookupTable = rhoLUT
xMLImageDataReader1Display.OSPRayScaleFunction = 'PiecewiseFunction'
xMLImageDataReader1Display.SelectOrientationVectors = 'None'
xMLImageDataReader1Display.ScaleFactor = 3.1
xMLImageDataReader1Display.SelectScaleArray = 'None'
xMLImageDataReader1Display.GlyphType = 'Arrow'
xMLImageDataReader1Display.GlyphTableIndexArray = 'None'
xMLImageDataReader1Display.GaussianRadius = 0.155
xMLImageDataReader1Display.SetScaleArray = [None, '']
xMLImageDataReader1Display.ScaleTransferFunction = 'PiecewiseFunction'
xMLImageDataReader1Display.OpacityArray = [None, '']
xMLImageDataReader1Display.OpacityTransferFunction = 'PiecewiseFunction'
xMLImageDataReader1Display.DataAxesGrid = 'GridAxesRepresentation'
xMLImageDataReader1Display.PolarAxes = 'PolarAxesRepresentation'
xMLImageDataReader1Display.ScalarOpacityUnitDistance = 1.7320508075688772
xMLImageDataReader1Display.ScalarOpacityFunction = rhoPWF
xMLImageDataReader1Display.Slice = 15

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
xMLImageDataReader1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for rhoLUT in view renderView1
rhoLUTColorBar = GetScalarBar(rhoLUT, renderView1)
rhoLUTColorBar.Title = 'rho'
rhoLUTColorBar.ComponentTitle = ''

# set color bar visibility
rhoLUTColorBar.Visibility = 1

# show color legend
xMLImageDataReader1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from glyph2
glyph2Display = Show(glyph2, renderView2)

# trace defaults for the display properties.
glyph2Display.Representation = 'Surface'
glyph2Display.ColorArrayName = ['POINTS', 'Type']
glyph2Display.LookupTable = typeLUT
glyph2Display.OSPRayScaleArray = 'Type'
glyph2Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph2Display.SelectOrientationVectors = 'None'
glyph2Display.ScaleFactor = 0.45
glyph2Display.SelectScaleArray = 'Type'
glyph2Display.GlyphType = 'Arrow'
glyph2Display.GlyphTableIndexArray = 'Type'
glyph2Display.GaussianRadius = 0.0225
glyph2Display.SetScaleArray = ['POINTS', 'Type']
glyph2Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph2Display.OpacityArray = ['POINTS', 'Type']
glyph2Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph2Display.DataAxesGrid = 'GridAxesRepresentation'
glyph2Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph2Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph2Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# show data from xMLImageDataReader2
xMLImageDataReader2Display = Show(xMLImageDataReader2, renderView2)

# get color transfer function/color map for 'rho1'
rho1LUT = GetColorTransferFunction('rho1')
rho1LUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.16834804952144622, 0.0, 0.0, 0.360784313725, 0.3355188399553299, 0.0, 1.0, 1.0, 0.5050441485643387, 0.0, 0.501960784314, 0.0, 0.6722149389982223, 1.0, 1.0, 0.0, 0.8405629885196686, 1.0, 0.380392156863, 0.0, 1.0089110380411148, 0.419607843137, 0.0, 0.0, 1.177259087562561, 0.878431372549, 0.301960784314, 0.301960784314]
rho1LUT.ColorSpace = 'RGB'
rho1LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'rho1'
rho1PWF = GetOpacityTransferFunction('rho1')
rho1PWF.Points = [0.0, 0.0, 0.5, 0.0, 1.177259087562561, 1.0, 0.5, 0.0]
rho1PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
xMLImageDataReader2Display.Representation = 'Volume'
xMLImageDataReader2Display.ColorArrayName = ['POINTS', 'rho1    ']
xMLImageDataReader2Display.LookupTable = rho1LUT
xMLImageDataReader2Display.OSPRayScaleArray = 'rho1    '
xMLImageDataReader2Display.OSPRayScaleFunction = 'PiecewiseFunction'
xMLImageDataReader2Display.SelectOrientationVectors = 'None'
xMLImageDataReader2Display.ScaleFactor = 3.1
xMLImageDataReader2Display.SelectScaleArray = 'None'
xMLImageDataReader2Display.GlyphType = 'Arrow'
xMLImageDataReader2Display.GlyphTableIndexArray = 'None'
xMLImageDataReader2Display.GaussianRadius = 0.155
xMLImageDataReader2Display.SetScaleArray = ['POINTS', 'rho1    ']
xMLImageDataReader2Display.ScaleTransferFunction = 'PiecewiseFunction'
xMLImageDataReader2Display.OpacityArray = ['POINTS', 'rho1    ']
xMLImageDataReader2Display.OpacityTransferFunction = 'PiecewiseFunction'
xMLImageDataReader2Display.DataAxesGrid = 'GridAxesRepresentation'
xMLImageDataReader2Display.PolarAxes = 'PolarAxesRepresentation'
xMLImageDataReader2Display.ScalarOpacityUnitDistance = 1.7320508075688772
xMLImageDataReader2Display.ScalarOpacityFunction = rho1PWF
xMLImageDataReader2Display.Slice = 15

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
xMLImageDataReader2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
xMLImageDataReader2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.177259087562561, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
xMLImageDataReader2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.177259087562561, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for rho1LUT in view renderView2
rho1LUTColorBar = GetScalarBar(rho1LUT, renderView2)
rho1LUTColorBar.WindowLocation = 'UpperRightCorner'
rho1LUTColorBar.Title = 'rho1    '
rho1LUTColorBar.ComponentTitle = ''

# set color bar visibility
rho1LUTColorBar.Visibility = 1

# show color legend
xMLImageDataReader2Display.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Type'
typePWF = GetOpacityTransferFunction('Type')
typePWF.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]
typePWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(glyph2)
# ----------------------------------------------------------------
