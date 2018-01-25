import harfang as hg
import helper_2d
from harfang_shortcut import *
import os, sys
import navier_stokes

if getattr(sys, 'frozen', False):
    # frozen
	exe_path = os.path.dirname(sys.executable)
	dir_ = sys._MEIPASS
	hg.LoadPlugins(dir_)
else:
    # unfrozen
	exe_path = dir_ = os.path.dirname(os.path.realpath(__file__))
	hg.LoadPlugins()

hg.SetLogIsDetailed(True)
hg.SetLogLevel(hg.LogAll)

plus = hg.GetPlus()
#plus.CreateWorkers()

# mount the system file driver
hg.MountFileDriver(hg.StdFileDriver())
hg.MountFileDriver(hg.StdFileDriver(dir_))
hg.MountFileDriver(hg.StdFileDriver(os.path.join(dir_, "assets")), "@assets")

plus.RenderInit(1600, 900, 4, hg.Windowed)

plus.SetBlend2D(hg.BlendAlpha)
plus.SetBlend3D(hg.BlendAlpha)

plus.SetRenderWindow()
plus.GetRendererAsync().SetVSync(False)
	

scn = plus.NewScene()
scn.AddComponent(hg.Environment())

simple_graphic_scene_overlay = hg.SimpleGraphicSceneOverlay(False)
scn.AddComponent(simple_graphic_scene_overlay)

cam = plus.AddCamera(scn, mat4.TranslationMatrix(vec3(0, 0, -5)))
cam.GetCamera().SetZNear(1)
cam.GetCamera().SetZFar(100000)  # 100km

fps = hg.FPSController(0, 0, -2)
fps.SetSmoothFactor(0.3, 0.3)

font = hg.RasterFont("@assets/fonts/Handel Gothic D Bold.ttf", 16)		

navier_stokes.setup()

while not plus.IsAppEnded():
	dt_sec = plus.UpdateClock()
	#fps.UpdateAndApplyToNode(cam, dt_sec)
	
	navier_stokes.simulation_step(simple_graphic_scene_overlay, hg.time_to_sec_f(dt_sec))

	fps.ApplyToNode(cam)	
	plus.GetRenderer().SetProjectionMatrix(cam.GetCamera().GetProjectionMatrix(plus.GetRenderer().GetAspectRatio()))
	plus.UpdateScene(scn, dt_sec)
	
	plus.Flip()
	plus.EndFrame()
	

os._exit(0)
