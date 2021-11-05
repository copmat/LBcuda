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
renderView1.ViewSize = [1521, 808]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [8.5, 8.5, 8.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [8.5, 8.5, 58.69097822426848]
renderView1.CameraFocalPoint = [8.5, 8.5, 8.5]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 12.99038105676658
renderView1.Background = [0.32, 0.34, 0.43]
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

# create a new 'XML Image Data Reader'
rho_0000000000 = XMLImageDataReader(FileName=['/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000001.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000002.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000003.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000004.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000005.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000006.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000007.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000008.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000009.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000010.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000011.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000012.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000013.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000014.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000015.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000016.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000017.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000018.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000019.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000020.vti'])
rho_0000000000.PointArrayStatus = ['rho1']

# create a new 'Legacy VTK Reader'
particle_0000 = LegacyVTKReader(FileNames=['/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000001.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000002.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000003.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000004.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000005.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000006.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000007.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000008.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000009.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000010.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000011.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000012.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000013.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000014.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000015.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000016.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000017.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000018.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000019.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000020.vtk'])

# create a new 'XML Image Data Reader'
isfluid_000000_0000 = XMLImageDataReader(FileName=['/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000001.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000002.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000003.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000004.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000005.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000006.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000007.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000008.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000009.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000010.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000011.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000012.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000013.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000014.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000015.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000016.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000017.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000018.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000019.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/isfluid_000000_000020.vti'])
isfluid_000000_0000.PointArrayStatus = ['isfluid']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from rho_0000000000
rho_0000000000Display = Show(rho_0000000000, renderView1)

# trace defaults for the display properties.
rho_0000000000Display.Representation = 'Outline'
rho_0000000000Display.ColorArrayName = [None, '']
rho_0000000000Display.OSPRayScaleArray = 'rho1'
rho_0000000000Display.OSPRayScaleFunction = 'PiecewiseFunction'
rho_0000000000Display.SelectOrientationVectors = 'None'
rho_0000000000Display.ScaleFactor = 1.5
rho_0000000000Display.SelectScaleArray = 'None'
rho_0000000000Display.GlyphType = 'Arrow'
rho_0000000000Display.GlyphTableIndexArray = 'None'
rho_0000000000Display.GaussianRadius = 0.075
rho_0000000000Display.SetScaleArray = ['POINTS', 'rho1']
rho_0000000000Display.ScaleTransferFunction = 'PiecewiseFunction'
rho_0000000000Display.OpacityArray = ['POINTS', 'rho1']
rho_0000000000Display.OpacityTransferFunction = 'PiecewiseFunction'
rho_0000000000Display.DataAxesGrid = 'GridAxesRepresentation'
rho_0000000000Display.PolarAxes = 'PolarAxesRepresentation'
rho_0000000000Display.ScalarOpacityUnitDistance = 1.7320508075688776
rho_0000000000Display.Slice = 7

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
rho_0000000000Display.OSPRayScaleFunction.Points = [1.2156107231930946e-09, 0.0, 0.5, 0.0, 0.00017424611723981798, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
rho_0000000000Display.ScaleTransferFunction.Points = [0.9997104406356812, 0.0, 0.5, 0.0, 1.000414490699768, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
rho_0000000000Display.OpacityTransferFunction.Points = [0.9997104406356812, 0.0, 0.5, 0.0, 1.000414490699768, 1.0, 0.5, 0.0]

# show data from isfluid_000000_0000
isfluid_000000_0000Display = Show(isfluid_000000_0000, renderView1)

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
isfluid_000000_0000Display.Representation = 'Volume'
isfluid_000000_0000Display.ColorArrayName = ['POINTS', 'isfluid']
isfluid_000000_0000Display.LookupTable = isfluidLUT
isfluid_000000_0000Display.OSPRayScaleArray = 'isfluid'
isfluid_000000_0000Display.OSPRayScaleFunction = 'PiecewiseFunction'
isfluid_000000_0000Display.SelectOrientationVectors = 'None'
isfluid_000000_0000Display.ScaleFactor = 1.5
isfluid_000000_0000Display.SelectScaleArray = 'None'
isfluid_000000_0000Display.GlyphType = 'Arrow'
isfluid_000000_0000Display.GlyphTableIndexArray = 'None'
isfluid_000000_0000Display.GaussianRadius = 0.075
isfluid_000000_0000Display.SetScaleArray = ['POINTS', 'isfluid']
isfluid_000000_0000Display.ScaleTransferFunction = 'PiecewiseFunction'
isfluid_000000_0000Display.OpacityArray = ['POINTS', 'isfluid']
isfluid_000000_0000Display.OpacityTransferFunction = 'PiecewiseFunction'
isfluid_000000_0000Display.DataAxesGrid = 'GridAxesRepresentation'
isfluid_000000_0000Display.PolarAxes = 'PolarAxesRepresentation'
isfluid_000000_0000Display.ScalarOpacityUnitDistance = 1.7320508075688776
isfluid_000000_0000Display.ScalarOpacityFunction = isfluidPWF
isfluid_000000_0000Display.Slice = 7

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
isfluid_000000_0000Display.OSPRayScaleFunction.Points = [1.2156107231930946e-09, 0.0, 0.5, 0.0, 0.00017424611723981798, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isfluid_000000_0000Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isfluid_000000_0000Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]

# show data from particle_0000
particle_0000Display = Show(particle_0000, renderView1)

# get color transfer function/color map for 'Type'
typeLUT = GetColorTransferFunction('Type')
typeLUT.RGBPoints = [1.0, 0.278431372549, 0.278431372549, 0.858823529412, 1.286, 0.0, 0.0, 0.360784313725, 1.5699999999999998, 0.0, 1.0, 1.0, 1.858, 0.0, 0.501960784314, 0.0, 2.142, 1.0, 1.0, 0.0, 2.428, 1.0, 0.380392156863, 0.0, 2.714, 0.419607843137, 0.0, 0.0, 3.0, 0.878431372549, 0.301960784314, 0.301960784314]
typeLUT.ColorSpace = 'RGB'
typeLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
particle_0000Display.Representation = 'Point Gaussian'
particle_0000Display.ColorArrayName = ['POINTS', 'Type']
particle_0000Display.LookupTable = typeLUT
particle_0000Display.OSPRayScaleArray = 'Type'
particle_0000Display.OSPRayScaleFunction = 'PiecewiseFunction'
particle_0000Display.SelectOrientationVectors = 'vel'
particle_0000Display.SelectScaleArray = 'Type'
particle_0000Display.GlyphType = 'Arrow'
particle_0000Display.GlyphTableIndexArray = 'Type'
particle_0000Display.GaussianRadius = 3.0
particle_0000Display.SetScaleArray = ['POINTS', 'Type']
particle_0000Display.ScaleTransferFunction = 'PiecewiseFunction'
particle_0000Display.OpacityArray = ['POINTS', 'Type']
particle_0000Display.OpacityTransferFunction = 'PiecewiseFunction'
particle_0000Display.DataAxesGrid = 'GridAxesRepresentation'
particle_0000Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
particle_0000Display.OSPRayScaleFunction.Points = [1.2156107231930946e-09, 0.0, 0.5, 0.0, 0.00017424611723981798, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
particle_0000Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
particle_0000Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for typeLUT in view renderView1
typeLUTColorBar = GetScalarBar(typeLUT, renderView1)
typeLUTColorBar.Title = 'Type'
typeLUTColorBar.ComponentTitle = ''

# set color bar visibility
typeLUTColorBar.Visibility = 1

# show color legend
particle_0000Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Type'
typePWF = GetOpacityTransferFunction('Type')
typePWF.Points = [1.0, 0.0, 0.5, 0.0, 3.0, 1.0, 0.5, 0.0]
typePWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(particle_0000)
# ----------------------------------------------------------------