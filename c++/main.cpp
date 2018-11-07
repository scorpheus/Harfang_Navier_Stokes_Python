// Harfang - Copyright 2001-2017 Emmanuel Julien. All Rights Reserved.

#include "engine/camera.h"
#include "engine/engine.h"
#include "engine/environment.h"
#include "engine/init.h"
#include "engine/plugin_system.h"
#include "engine/imgui.h"
#include "engine/plus.h"
#include "engine/renderer_async.h"
#include "engine/scene.h"
#include "engine/simple_graphic_scene_overlay.h"
#include "engine/transform.h"
#include "foundation/filesystem.h"
#include "foundation/matrix3.h"
#include "foundation/std_file_driver.h"

#include "navier_stokes.h"

using namespace hg;

int main(int argc, char *argv[]) {
	Init(argv[0]);
	LoadPlugins();
	g_fs.get().Mount(std::make_shared<StdFileDriver>("assets"), "@assets/");

	SetLogIsDetailed(true);
	SetLogLevel(LogAll);

	auto &plus = g_plus.get();
	plus.RenderInit(1600, 900, 4, Window::Windowed);
	plus.GetRendererAsync()->SetVSync(false);

	plus.SetBlend2D(BlendAlpha);
	plus.SetBlend3D(BlendAlpha);

	auto scn = plus.NewScene(false, false);
	scn->AddComponent(std::make_shared<Environment>());

	auto simple_graphic_scene_overlay = std::make_shared<SimpleGraphicSceneOverlay>(false);
	scn->AddComponent(simple_graphic_scene_overlay);

	auto cam = plus.AddCamera(*scn, Matrix4::TranslationMatrix(Vector3(0, 0.0f, -2)));
	cam->GetComponent<Camera>()->SetZNear(0.1);
	cam->GetComponent<Camera>()->SetZFar(10000);

	OrbitalController controller;
	controller.Reset(Vector3(0, 0.f, 0), 4.f);

	setup_navier_stokes();
	scn->GetFrameCompleteSignal().Connect(&on_frame_complete);

	while (1) {
		plus.Clear(Color(0.1f, 0.1f, 0.1f));
		auto dt_sec = plus.UpdateClock();
		if (!plus.KeyDown(KeyLCtrl) && !ImGui::GetIO().WantCaptureMouse)
			controller.UpdateAndApplyToNode(*cam, dt_sec);

		auto mat = Matrix3::LookAt(cam->GetComponent<Transform>()->GetPosition().Normalized() * -1.f);
		simulation_step(simple_graphic_scene_overlay, mat, time_to_sec_f(dt_sec));

		plus.UpdateScene(*scn, dt_sec);
		plus.Flip();
		plus.EndFrame();
	}

	plus.RenderUninit();
}
