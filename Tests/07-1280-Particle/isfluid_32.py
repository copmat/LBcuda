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
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [16.5, 16.5, 16.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [13.988794436664154, 70.8977033368818, 104.78412199547739]
renderView1.CameraFocalPoint = [16.5, 16.5, 16.5]
renderView1.CameraViewUp = [0.006584003535752539, 0.8514262999569452, -0.5244329381713809]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 26.846787517317598
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

# create a new 'XML Image Data Reader'
step__isfluid_00 = XMLImageDataReader(FileName=['/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000000.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000010.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000020.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000030.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000040.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000050.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000060.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000070.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000080.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000090.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000100.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000110.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000120.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000130.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000140.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000150.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000160.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000170.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000180.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000190.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000200.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000210.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000220.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000230.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000240.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000250.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000260.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000270.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000280.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000290.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000300.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000310.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000320.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000330.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000340.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000350.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000360.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000370.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000380.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000390.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000400.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000410.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000420.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000430.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000440.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000450.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000460.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000470.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000480.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000490.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000500.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000510.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000520.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000530.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000540.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000550.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000560.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000570.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000580.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000590.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000600.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000610.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000620.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000630.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000640.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000650.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000660.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000670.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000680.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000690.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000700.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000710.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000720.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000730.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000740.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000750.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000760.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000770.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000780.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000790.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000800.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000810.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000820.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000830.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000840.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000850.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000860.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000870.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000880.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000890.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000900.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000910.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000920.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000930.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000940.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000950.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000960.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000970.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000980.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_000990.vti', '/home/ziko/IIT/hackathon/lbcuda/Tests/07-1280-Particle/output/step__isfluid_001000.vti'])
step__isfluid_00.PointArrayStatus = ['step__isfluid  ']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from step__isfluid_00
step__isfluid_00Display = Show(step__isfluid_00, renderView1)

# get color transfer function/color map for 'step__isfluid'
step__isfluidLUT = GetColorTransferFunction('step__isfluid')
step__isfluidLUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.715, 0.0, 0.0, 0.360784313725, 1.4249999999999998, 0.0, 1.0, 1.0, 2.145, 0.0, 0.501960784314, 0.0, 2.8549999999999995, 1.0, 1.0, 0.0, 3.57, 1.0, 0.380392156863, 0.0, 4.285, 0.419607843137, 0.0, 0.0, 5.0, 0.878431372549, 0.301960784314, 0.301960784314]
step__isfluidLUT.ColorSpace = 'RGB'
step__isfluidLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'step__isfluid'
step__isfluidPWF = GetOpacityTransferFunction('step__isfluid')
step__isfluidPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 0.05147058889269829, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]
step__isfluidPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
step__isfluid_00Display.Representation = 'Volume'
step__isfluid_00Display.ColorArrayName = ['POINTS', 'step__isfluid  ']
step__isfluid_00Display.LookupTable = step__isfluidLUT
step__isfluid_00Display.OSPRayScaleArray = 'step__isfluid  '
step__isfluid_00Display.OSPRayScaleFunction = 'PiecewiseFunction'
step__isfluid_00Display.SelectOrientationVectors = 'None'
step__isfluid_00Display.ScaleFactor = 3.1
step__isfluid_00Display.SelectScaleArray = 'None'
step__isfluid_00Display.GlyphType = 'Arrow'
step__isfluid_00Display.GlyphTableIndexArray = 'None'
step__isfluid_00Display.GaussianRadius = 0.155
step__isfluid_00Display.SetScaleArray = ['POINTS', 'step__isfluid  ']
step__isfluid_00Display.ScaleTransferFunction = 'PiecewiseFunction'
step__isfluid_00Display.OpacityArray = ['POINTS', 'step__isfluid  ']
step__isfluid_00Display.OpacityTransferFunction = 'PiecewiseFunction'
step__isfluid_00Display.DataAxesGrid = 'GridAxesRepresentation'
step__isfluid_00Display.PolarAxes = 'PolarAxesRepresentation'
step__isfluid_00Display.ScalarOpacityUnitDistance = 1.7320508075688772
step__isfluid_00Display.ScalarOpacityFunction = step__isfluidPWF
step__isfluid_00Display.Slice = 15

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
step__isfluid_00Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.000000238418579, 1.0, 0.5, 0.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(step__isfluid_00)
# ----------------------------------------------------------------