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
renderView1.ViewSize = [1471, 808]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [64.53834116458893, 64.29839783906937, 65.02098071575165]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [239.84314544757265, -113.28438782235219, 340.24852681455013]
renderView1.CameraFocalPoint = [56.002644929073945, 53.15429428969192, 63.98716535208355]
renderView1.CameraViewUp = [0.8686204592091055, 0.2716927418037635, -0.4143447259161344]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 116.26177328818973
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

# create a new 'Legacy VTK Reader'
particle_0 = LegacyVTKReader(FileNames=['/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_000900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_001900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_002900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_003900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_004900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_005900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_006900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_007900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_008900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009000.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009100.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009200.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009300.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009400.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009500.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009600.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009700.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009800.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_009900.vtk', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/particle_010000.vtk'])

# create a new 'XML Image Data Reader'
rho_0000000 = XMLImageDataReader(FileName=['/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.000900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.001900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.002900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.003900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.004900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.005900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.006900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.007900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.008900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009000.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009100.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009200.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009300.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009400.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009500.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009600.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009700.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009800.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.009900.vti', '/home/ziko/IIT/lbsoft-cudafor/Tests/07-1280-Particle/output/rho_000000.010000.vti'])
rho_0000000.PointArrayStatus = ['rho1']

# create a new 'Contour'
contour1 = Contour(Input=rho_0000000)
contour1.ContourBy = ['POINTS', 'rho1']
contour1.ComputeGradients = 1
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Glyph'
glyph1 = Glyph(Input=particle_0,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'Vectx']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 5.0
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from particle_0
particle_0Display = Show(particle_0, renderView1)

# trace defaults for the display properties.
particle_0Display.Representation = 'Point Gaussian'
particle_0Display.AmbientColor = [0.32941176470588235, 0.32941176470588235, 0.32941176470588235]
particle_0Display.ColorArrayName = ['POINTS', '']
particle_0Display.DiffuseColor = [0.32941176470588235, 0.32941176470588235, 0.32941176470588235]
particle_0Display.OSPRayScaleArray = 'Type'
particle_0Display.OSPRayScaleFunction = 'PiecewiseFunction'
particle_0Display.SelectOrientationVectors = 'vel'
particle_0Display.ScaleFactor = 12.774504500627518
particle_0Display.SelectScaleArray = 'Type'
particle_0Display.GlyphType = 'Arrow'
particle_0Display.GlyphTableIndexArray = 'Type'
particle_0Display.GaussianRadius = 2.0
particle_0Display.SetScaleArray = ['POINTS', 'Type']
particle_0Display.ScaleTransferFunction = 'PiecewiseFunction'
particle_0Display.OpacityArray = ['POINTS', 'Type']
particle_0Display.OpacityTransferFunction = 'PiecewiseFunction'
particle_0Display.DataAxesGrid = 'GridAxesRepresentation'
particle_0Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
particle_0Display.OSPRayScaleFunction.Points = [1.2156107231930946e-09, 0.0, 0.5, 0.0, 0.00017424611723981798, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
particle_0Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1000.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
particle_0Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1000.0, 1.0, 0.5, 0.0]

# show data from glyph1
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [0.6666666666666666, 0.0, 0.0]
glyph1Display.ColorArrayName = ['POINTS', '']
glyph1Display.DiffuseColor = [0.6666666666666666, 0.0, 0.0]
glyph1Display.OSPRayScaleArray = 'Type'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 16.485483551025393
glyph1Display.SelectScaleArray = 'Type'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'Type'
glyph1Display.GaussianRadius = 0.8242741775512695
glyph1Display.SetScaleArray = ['POINTS', 'Type']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'Type']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [1.2156107231930946e-09, 0.0, 0.5, 0.0, 0.00017424611723981798, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1000.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1000.0, 1.0, 0.5, 0.0]

# show data from contour1
contour1Display = Show(contour1, renderView1)

# get color transfer function/color map for 'rho1'
rho1LUT = GetColorTransferFunction('rho1')
rho1LUT.RGBPoints = [0.0, 0.278431372549, 0.278431372549, 0.858823529412, 0.15482816469669342, 0.0, 0.0, 0.360784313725, 0.3085736149549484, 0.0, 1.0, 1.0, 0.46448449409008025, 0.0, 0.501960784314, 0.0, 0.6182299443483352, 1.0, 1.0, 0.0, 0.7730581090450287, 1.0, 0.380392156863, 0.0, 0.9278862737417221, 0.419607843137, 0.0, 0.0, 1.0827144384384155, 0.878431372549, 0.301960784314, 0.301960784314]
rho1LUT.ColorSpace = 'RGB'
rho1LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'rho1']
contour1Display.LookupTable = rho1LUT
contour1Display.Specular = 1.0
contour1Display.OSPRayScaleArray = 'rho1'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'Gradients'
contour1Display.ScaleFactor = 12.700000000000001
contour1Display.SelectScaleArray = 'rho1'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'rho1'
contour1Display.GaussianRadius = 0.635
contour1Display.SetScaleArray = ['POINTS', 'rho1']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'rho1']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [1.2156107231930946e-09, 0.0, 0.5, 0.0, 0.00017424611723981798, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for rho1LUT in view renderView1
rho1LUTColorBar = GetScalarBar(rho1LUT, renderView1)
rho1LUTColorBar.Title = 'rho1'
rho1LUTColorBar.ComponentTitle = ''

# set color bar visibility
rho1LUTColorBar.Visibility = 1

# show color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'rho1'
rho1PWF = GetOpacityTransferFunction('rho1')
rho1PWF.Points = [0.0, 0.0, 0.5, 0.0, 1.0827144384384155, 1.0, 0.5, 0.0]
rho1PWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(particle_0)
# ----------------------------------------------------------------