#include "navier_stokes.h"
#include "engine/imgui.h"
#include "engine/iso_surface.h"
#include "engine/plus.h"
#include "engine/render_geometry.h"
#include "foundation/random.h"
#include "engine/simple_graphic_scene_overlay.h"
#include "foundation/binary_data.h"
#include "foundation/matrix4.h"
#include <array>

static const int navier_width = 20;
static const int navier_height = 30;
static const int navier_depth = 20;

static const int navier_width_plus2 = navier_width + 2;
static const int navier_height_plus2 = navier_height +2;
static const int navier_depth_plus2 = navier_depth + 2;
static const int size = navier_width_plus2 * navier_height_plus2 * navier_depth_plus2;

typedef std::array<float, size> farray;

int first_use = 10;
bool ui_draw_line = false;
bool ui_draw_iso = false;
bool ui_draw_billboard = false;
bool ui_draw_raymarch = true;
float ui_iso_level = 0.005f;
farray iso_array;
farray raymarch_array;

float visc = 0.0;
float diff = 0.0;
float base_force_u = 5.f;
float base_force_v = 60.f; 
float base_force_w = 60.f;
float base_force_dens = 1000.f;

farray u{};
farray u_prev{};
farray v{};
farray v_prev{};
farray w{};
farray w_prev{};
farray dens{};
farray dens_prev{};

#define IX(i, j, k) (int(i) + navier_width_plus2 * int(j) + navier_width_plus2 * navier_height_plus2 * int(k))

void add_source(farray &x, farray &s, const float &dt) {
#pragma omp parallel for
	for (int i = 0; i < size; ++i)
		x[i] += dt * s[i];
}

void set_bnd(const int &N, const int &O, const int &P, const int &b, farray &x) {
#pragma omp parallel for
	for (int i = 1; i <= O + 1; ++i)
		for (int j = 1; j <= P + 1; ++j) {
			if (b == 1) {
				x[IX(0, i, j)] = x[IX(1, i, j)] * -1.f;
				x[IX(N + 1, i, j)] = x[IX(N, i, j)] * -1.f;
			}
			else {
				x[IX(0, i, j)] = x[IX(1, i, j)];
				x[IX(N + 1, i, j)] = x[IX(N, i, j)];
			}
		}
#pragma omp parallel for
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= P + 1; ++j) {
			if (b == 2) {
				x[IX(i, 0, j)] = x[IX(i, 1, j)] * -1.f;
				x[IX(i, O + 1, j)] = x[IX(i, O, j)] * -1.f;
			}
			else {
				x[IX(i, 0, j)] = x[IX(i, 1, j)];
				x[IX(i, O + 1, j)] = x[IX(i, O, j)];
			}
		}
#pragma omp parallel for
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= O + 1; ++j) {
			if (b == 3) {
				x[IX(i, j, 0)] = x[IX(i, j, 1)] * -1.f;
				x[IX(i, j, P + 1)] = x[IX(i, j, P)] * -1.f;
			}
			else {
				x[IX(i, j, 0)] = x[IX(i, j, 1)];
				x[IX(i, j, P + 1)] = x[IX(i, j, P)];
			}
		}

	x[IX(0, 0, 0)] = (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]) / 3.f;
	x[IX(0, O + 1, 0)] = (x[IX(1, O + 1, 0)] + x[IX(0, O, 0)] + x[IX(0, O + 1, 1)]) / 3.f;
	x[IX(N + 1, 0, 0)] = (x[IX(N, 0, 0)] + x[IX(N + 1, 1, 0)] + x[IX(N + 1, 0, 1)]) / 3.f;
	x[IX(N + 1, O + 1, 0)] = (x[IX(N, O + 1, 0)] + x[IX(N + 1, O, 0)] + x[IX(N + 1, O + 1, 1)]) / 3.f;

	x[IX(0, 0, P + 1)] = (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, P)]) / 3.f;
	x[IX(0, O + 1, P + 1)] = (x[IX(1, O + 1, 0)] + x[IX(0, O, 0)] + x[IX(0, O + 1, P)]) / 3.f;
	x[IX(N + 1, 0, P + 1)] = (x[IX(N, 0, 0)] + x[IX(N + 1, 1, 0)] + x[IX(N + 1, 0, P)]) / 3.f;
	x[IX(N + 1, O + 1, P + 1)] = (x[IX(N, O + 1, 0)] + x[IX(N + 1, O, 0)] + x[IX(N + 1, O + 1, P)]) / 3.f;
}

void lin_solve(const int &N, const int &O, const int &P, const int &b, farray &x, farray &x0, const float &a, const float &c) {
	for (int k = 0; k < 20; ++k) {
#pragma omp parallel for
		for (int i = 1; i <= N + 1; ++i)
			for (int j = 1; j <= O + 1; ++j)
				for (int k = 1; k <= P + 1; ++k)
					x[IX(i, j, k)] = (x0[IX(i, j, k)] + a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] + x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] + x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / c;

		set_bnd(N, O, P, b, x);
	}
}

void diffuse(const int &N, const int &O, const int &P, const int &b, farray &x, farray &x0, const float &diff, const float &dt) {
	float a = dt * diff * N * O * P;
	lin_solve(N, O, P, b, x, x0, a, 1 + 6 * a);
}

void advect(const int &N, const int &O, const int &P, const int &b, farray &d, farray &d0, farray &u, farray &v, farray &w, const float &dt) {
	float dt0 = dt * N;
#pragma omp parallel for
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= O + 1; ++j)
			for (int k = 1; k <= P + 1; ++k) {
				float x = i - dt0 * u[IX(i, j, k)];
				float y = j - dt0 * v[IX(i, j, k)];
				float z = k - dt0 * w[IX(i, j, k)];
				if (x < 0.5f)
					x = 0.5f;
				if (x > N + 0.5f)
					x = N + 0.5f;
				float i0 = floor(x);
				float i1 = i0 + 1.f;

				if (y < 0.5f)
					y = 0.5f;
				if (y > O + 0.5f)
					y = O + 0.5f;
				float j0 = floor(y);
				float j1 = j0 + 1;

				if (z < 0.5f)
					z = 0.5f;
				if (z > P + 0.5f)
					z = P + 0.5f;
				float k0 = floor(z);
				float k1 = k0 + 1;

				float s1 = x - i0;
				float s0 = 1 - s1;

				float t1 = y - j0;
				float t0 = 1 - t1;

				float r1 = z - k0;
				float r0 = 1 - r1;
				d[IX(i, j, k)] = r0 * (s0 * (t0 * d0[IX(i0, j0, k0)] + t1 * d0[IX(i0, j1, k0)]) + s1 * (t0 * d0[IX(i1, j0, k0)] + t1 * d0[IX(i1, j1, k0)])) +
								 r1 * (s0 * (t0 * d0[IX(i0, j0, k1)] + t1 * d0[IX(i0, j1, k1)]) + s1 * (t0 * d0[IX(i1, j0, k1)] + t1 * d0[IX(i1, j1, k1)]));
			}

	set_bnd(N, O, P, b, d);
}

void dens_step(const int &N, const int &O, const int &P, farray &x, farray &x0, farray &u, farray &v, farray &w, const float &diff, const float &dt) {
	add_source(x, x0, dt);
	diffuse(N, O, P, 0, x0, x, diff, dt);
	advect(N, O, P, 0, x, x0, u, v, w, dt);
}

void project(const int &N, const int &O, const int &P, farray &u, farray &v, farray &w, farray &p, farray &div) {
#pragma omp parallel for
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= O + 1; ++j)
			for (int k = 1; k <= P + 1; ++k)
				div[IX(i, j, k)] = (u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)] + v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)] + w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]) / (N * -3.f);
	p.fill(0);

	set_bnd(N, O, P, 0, div);
	set_bnd(N, O, P, 0, p);

	lin_solve(N, O, P, 0, p, div, 1, 6);
#pragma omp parallel for
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= O + 1; ++j)
			for (int k = 1; k <= P + 1; ++k) {
				u[IX(i, j, k)] -= 0.5f * N * (p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
				v[IX(i, j, k)] -= 0.5f * O * (p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
				w[IX(i, j, k)] -= 0.5f * P * (p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
			}
	set_bnd(N, O, P, 1, u);
	set_bnd(N, O, P, 2, v);
	set_bnd(N, O, P, 3, w);
}

void vel_step(const int &N, const int &O, const int &P, farray &u, farray &v, farray &w, farray &u0, farray &v0, farray &w0, const float &visc, const float &dt) {
	add_source(u, u0, dt);
	add_source(v, v0, dt);
	add_source(w, w0, dt);

	diffuse(N, O, P, 1, u0, u, visc, dt);
	diffuse(N, O, P, 2, v0, v, visc, dt);
	diffuse(N, O, P, 3, w0, w, visc, dt);
	project(N, O, P, u0, v0, w0, u, v);

	advect(N, O, P, 1, u, u0, u0, v0, w0, dt);
	advect(N, O, P, 2, v, v0, u0, v0, w0, dt);
	advect(N, O, P, 3, w, w0, u0, v0, w0, dt);
	project(N, O, P, u, v, w, u0, v0);
}

void get_from_UI(farray &dens_prev, farray &u_prev, farray &v_prev, farray &w_prev) {
	dens_prev.fill(0.f);
	u_prev.fill(0.f);
	v_prev.fill(0.f);
	w_prev.fill(0.f);

	auto &plus = hg::g_plus.get();

	hg::Vector2 mouse_pos, mouse_dt;
	plus.GetMousePos(mouse_pos.x, mouse_pos.y);
	plus.GetMouseDt(mouse_dt.x, mouse_dt.y);

	auto size_window = plus.GetRenderWindowSize(plus.GetRenderWindow());
	if (plus.MouseButtonDown(hg::Button0) && plus.KeyDown(hg::KeyLCtrl))
		dens_prev[IX((mouse_pos.x / size_window.x) * navier_width, (mouse_pos.y / size_window.y) * navier_height, navier_depth / 2.f)] = base_force_dens;
	if (plus.MouseButtonDown(hg::Button1)) {
		u_prev[IX((mouse_pos.x / size_window.x) * navier_width, (mouse_pos.y / size_window.y) * navier_height, navier_depth / 2.f)] = base_force_u * mouse_dt.x;
		v_prev[IX((mouse_pos.x / size_window.x) * navier_width, (mouse_pos.y / size_window.y) * navier_height, navier_depth / 2.f)] = base_force_v * mouse_dt.y;
	}
	if (plus.KeyDown(hg::KeySpace) || first_use > 0) {
		first_use--;
		dens_prev[IX(navier_width / 2.f, 1, navier_depth / 2.f)] = base_force_dens;
		v_prev[IX(navier_width / 2.f, 3, navier_depth / 2.f)] = base_force_v * 10.0f;
		u_prev[IX(navier_width / 2.f, 3, navier_depth / 2.f)] = base_force_u * hg::FRand(2.0f);
		w_prev[IX(navier_width / 2.f, 3, navier_depth / 2.f)] = base_force_w * hg::FRand(2.0f);
	}
}

void draw_line(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &scene_simple_graphic, const hg::Vector3 &a, const hg::Vector3 &b, const hg::Color &color = hg::Color::White) { scene_simple_graphic->Line(a.x, a.y, a.z, b.x, b.y, b.z, color, color); }

void draw_quad(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &scene_simple_graphic, const hg::Matrix4 &mat, const float &width, const float &height, std::shared_ptr<hg::Texture> texture, const hg::Color &color = hg::Color::White) {
	auto pos = mat.GetTranslation();
	auto axis_x = mat.GetX() * width / 2;
	auto axis_y = mat.GetY() * height / 2;
	auto a = pos - axis_x - axis_y;
	auto b = pos - axis_x + axis_y;
	auto c = pos + axis_x + axis_y;
	auto d = pos + axis_x - axis_y;

	scene_simple_graphic->Quad(a.x, a.y, a.z,
		b.x, b.y, b.z,
		c.x, c.y, c.z,
		d.x, d.y, d.z,
		0, 0, 1, 1,
		texture,
		color, color, color, color);
}

void draw_iso_surface(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, farray &dens) {
	auto &plus = hg::g_plus.get();
	auto iso = std::make_shared<hg::IsoSurface>();

	#pragma omp parallel for
	for (int i = 0; i < size; ++i)
		iso_array[i] = dens[i] / 100.f;

	hg::PolygoniseIsoSurface(navier_width, navier_depth, navier_height, iso_array.data(), ui_iso_level, *iso, hg::Vector3(1.f / navier_width, 1.f / navier_width, 1.f / navier_width));

	auto mat = plus.LoadMaterial("@assets/default.mat");
	auto geo = std::make_shared<hg::RenderGeometry>();
	hg::IsoSurfaceToRenderGeometry(plus.GetRenderSystem(), iso, geo, mat);

	simple_graphic_scene_overlay->SetDepthTest(false);
	simple_graphic_scene_overlay->Geometry(0.5f, -0.5f, -0.5f, -1.57f, -1.57f, -1.57f, 1, 1, 1, geo);
}

void draw_dens(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, hg::Matrix3 billboard_mat, const int &N, const int &O, const int &P, farray &dens, farray &u, farray &v, farray &w) {
	float scale = 1.f;
	if (ui_draw_billboard) {
		simple_graphic_scene_overlay->SetDepthTest(false);
		simple_graphic_scene_overlay->SetBlendMode(hg::BlendAdditive);
		for (int i = 1; i <= N + 1; ++i)
			for (int j = 1; j <= O + 1; ++j)
				for (int k = 1; k <= P + 1; ++k) {
					hg::Vector3 p(i / float(N) - 0.5f, j / float(N) - 0.5f, k / float(N) - 0.5f);
					float density = dens[IX(i, j, k)]/100.f;
					draw_quad(simple_graphic_scene_overlay, hg::Matrix4::TransformationMatrix(p, billboard_mat), scale / N, scale / N, nullptr, hg::Color(0, density, density, 1));
				}
		simple_graphic_scene_overlay->SetBlendMode(hg::BlendOpaque);
	}

	if (ui_draw_line)
		for (int i = 1; i <= N + 1; ++i)
			for (int j = 1; j <= O + 1; ++j)
				for (int k = 1; k <= P + 1; ++k) {
					hg::Vector3 p(i / float(N) - 0.5f, j / float(N) - 0.5f, k / float(N) - 0.5f);
					draw_line(simple_graphic_scene_overlay, p, p + hg::Vector3(u[IX(i, j, k)], v[IX(i, j, k)], w[IX(i, j, k)]));
				}
}

#include <engine/renderer.h>
#include <engine/render_system.h>
#include <engine/set_shader_engine_values.h>

std::shared_ptr<hg::GpuShader> shader;
std::shared_ptr<hg::Texture> tex;
bool loaded = false;

void setup_raymarch()
{
	hg::Plus &plus = hg::g_plus.get();
	std::shared_ptr<hg::Renderer> renderer = plus.GetRenderer();

	tex = renderer->NewTexture();
	// Load shader.
	shader = renderer->LoadShader("@assets/raymarch_shader.isl");
}

void setup_navier_stokes()
{
	setup_raymarch();
}

void on_frame_complete(hg::Scene &scn, hg::RenderSystem &render_system)
{
	if (!ui_draw_raymarch)
		return;

	hg::Plus &plus = hg::g_plus.get();
	// Retrieve renderer.
	std::shared_ptr<hg::Renderer> renderer = plus.GetRenderer();

	// Render a fullscreen quad and mix frames.
	renderer->EnableDepthTest(false);
	renderer->EnableBlending(true);
	renderer->SetShader(shader);

	float width = navier_width_plus2;
	float height = navier_height_plus2;
	float depth = navier_depth_plus2;
	float size_array = size;
	
	renderer->SetShaderFloat(renderer->GetShaderVariable("width"), &width);
	renderer->SetShaderFloat(renderer->GetShaderVariable("height"), &height);
	renderer->SetShaderFloat(renderer->GetShaderVariable("depth"), &depth);
	renderer->SetShaderFloat(renderer->GetShaderVariable("size_array"), &size_array);

	// update the raymarch texture
#pragma omp parallel for
	for (int i = 0; i < size; ++i)
		raymarch_array[i] = dens[i] / 100.f;

	renderer->CreateTexture(*tex, raymarch_array.data(), raymarch_array.size() * sizeof(float), raymarch_array.size(), 1, hg::TextureDepthF, hg::TextureNoAA, hg::TextureIsShaderResource, false);
	renderer->SetShaderTexture(renderer->GetShaderVariable("array_tex"), *tex);

	render_system.DrawFullscreenQuad(hg::Vector2(0, 0));
}

void simulation_step(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, hg::Matrix3 billboard_mat, const float &dt) {
	if (ImGui::Begin("params")) {
		ImGui::Text("Press Space for small puff");
		ImGui::SliderFloat("visc", &visc, 0.0, 1.f);
		ImGui::SliderFloat("diff", &diff, 0.0, 1.f);
		ImGui::SliderFloat("base_force_u", &base_force_u, 0.0, 1000.0f);
		ImGui::SliderFloat("base_force_v", &base_force_v, 0.0, 1000.0f);
		ImGui::SliderFloat("base_force_w", &base_force_w, 0.0, 1000.0f);
		ImGui::SliderFloat("base_force_dens", &base_force_dens, 0.0, 1000.0f);
		ImGui::Checkbox("Draw Vector Field", &ui_draw_line);
		ImGui::Checkbox("Draw iso surface", &ui_draw_iso);
		if (ui_draw_iso)
			ImGui::SliderFloat("iso level", &ui_iso_level, 0, 0.1f);
		
		ImGui::Checkbox("Draw billboard", &ui_draw_billboard); 
		ImGui::Checkbox("Draw shader raymarch", &ui_draw_raymarch);

		if(ImGui::Button("Clear"))
		{
			#pragma omp parallel for
			for (int i = 0; i < size; ++i)
				dens[i] = 0;
		}
	}
	ImGui::End();

	get_from_UI(dens_prev, u_prev, v_prev, w_prev);
	vel_step(navier_width, navier_height, navier_depth, u, v, w, u_prev, v_prev, w_prev, visc/1000.f, dt);
	dens_step(navier_width, navier_height, navier_depth, dens, dens_prev, u, v, w, diff/1000.f, dt);
	draw_dens(simple_graphic_scene_overlay, billboard_mat, navier_width, navier_height, navier_depth, dens, u, v, w);
	if (ui_draw_iso)
		draw_iso_surface(simple_graphic_scene_overlay, dens);
}

