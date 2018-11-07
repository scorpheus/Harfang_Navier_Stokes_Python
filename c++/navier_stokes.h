#pragma once
#include <memory>

namespace hg {
	class SimpleGraphicSceneOverlay;
	class Matrix3;
	class Scene;
	class RenderSystem;
}

void simulation_step(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, hg::Matrix3 billboard_mat, const float &dt);
void setup_navier_stokes();
void on_frame_complete(hg::Scene &scn, hg::RenderSystem &render_system);
