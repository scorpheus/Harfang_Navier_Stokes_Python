#pragma once
#include <memory>

namespace hg {
	class SimpleGraphicSceneOverlay;
	class Matrix3;
}
void simulation_step(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, hg::Matrix3 billboard_mat, const float &dt);