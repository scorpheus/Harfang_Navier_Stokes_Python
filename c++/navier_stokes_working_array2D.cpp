#include "navier_stokes.h"
#include "engine/imgui.h"
#include "engine/plus.h"
#include "engine/simple_graphic_scene_overlay.h"
#include "foundation/matrix4.h"
#include <array>

bool ui_draw_line = false;

float visc = 0.0;
float diff = 0.0;
float base_force_u = 5.f;
float base_force_v = 5.f;
float base_force_dens = 100.f;

static const int N = 128;
static const int NPlus2 = N+2;
static const int size = (N + 2) * (N + 2);

typedef std::array<float, size> farray;

farray u{};
farray u_prev{};
farray v{};
farray v_prev{};
farray dens{};
farray dens_prev{};

//int IX(const int &i, const int &j) { return ((i)+(N + 2) * (j)); }
#define IX(i, j) (int(i) + NPlus2 * int(j))

void add_source(const int &N, farray &x, farray &s, const float &dt) {
	for (int i = 0; i < size; ++i)
		x[i] += dt * s[i];
}

void set_bnd(const int &N, const int &b, farray &x) {
	for (int i = 1; i <= N + 1; ++i) {
		if (b == 1) {
			x[IX(0, i)] = x[IX(1, i)] * -1.f;
			x[IX(N + 1, i)] = x[IX(N, i)] * -1.f;
		} else {
			x[IX(0, i)] = x[IX(1, i)];
			x[IX(N + 1, i)] = x[IX(N, i)];
		}

		if (b == 2) {
			x[IX(i, 0)] = x[IX(i, 1)] * -1.f;
			x[IX(i, N + 1)] = x[IX(i, N)] * -1.f;
		} else {
			x[IX(i, 0)] = x[IX(i, 1)];
			x[IX(i, N + 1)] = x[IX(i, N)];
		}
	}

	x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
	x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
	x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
}

void lin_solve(const int &N, const int &b, farray &x, farray &x0, const float &a, const float &c) {
	for (int k = 0; k < 20; ++k) {
		for (int i = 1; i <= N + 1; ++i)
			for (int j = 1; j <= N + 1; ++j)
				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;

		set_bnd(N, b, x);
	}
}

void diffuse(const int &N, const int &b, farray &x, farray &x0, const float &diff, const float &dt) {
	float a = dt * diff * N * N;
	lin_solve(N, b, x, x0, a, 1 + 4 * a);
}

void advect(const int &N, const int &b, farray &d, farray &d0, farray &u, farray &v, const float &dt) {
	float dt0 = dt * N;
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= N + 1; ++j) {
			float x = i - dt0 * u[IX(i, j)];
			float y = j - dt0 * v[IX(i, j)];
			if (x < 0.5f)
				x = 0.5f;
			if (x > N + 0.5f)
				x = N + 0.5f;
			float i0 = floor(x);
			float i1 = i0 + 1.f;
			if (y < 0.5f)
				y = 0.5f;
			if (y > N + 0.5f)
				y = N + 0.5f;
			float j0 = floor(y);
			float j1 = j0 + 1;
			float s1 = x - i0;
			float s0 = 1 - s1;
			float t1 = y - j0;
			float t0 = 1 - t1;
			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}

	set_bnd(N, b, d);
}

void dens_step(const int &N, farray &x, farray &x0, farray &u, farray &v, const float &diff, const float &dt) {
	add_source(N, x, x0, dt);
	diffuse(N, 0, x0, x, diff, dt);
	advect(N, 0, x, x0, u, v, dt);
}

void project(const int &N, farray &u, farray &v, farray &p, farray &div) {
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= N + 1; ++j) 
			div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
	p.fill(0);

	set_bnd(N, 0, div);
	set_bnd(N, 0, p);

	lin_solve(N, 0, p, div, 1, 4);
	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= N + 1; ++j) {
			u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
			v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
		}
	set_bnd(N, 1, u);
	set_bnd(N, 2, v);
}

void vel_step(const int &N, farray &u, farray &v, farray &u0, farray &v0, const float &visc, const float &dt) {
	add_source(N, u, u0, dt);
	add_source(N, v, v0, dt);

	diffuse(N, 1, u0, u, visc, dt);
	diffuse(N, 2, v0, v, visc, dt);
	project(N, u0, v0, u, v);

	advect(N, 1, u, u0, u0, v0, dt);
	advect(N, 2, v, v0, u0, v0, dt);
	project(N, u, v, u0, v0);
}

void get_from_UI(farray &dens_prev, farray &u_prev, farray &v_prev) {
	dens_prev.fill(0.f);
	u_prev.fill(0.f);
	v_prev.fill(0.f);

	auto &plus = hg::g_plus.get();

	hg::Vector2 mouse_pos, mouse_dt;
	plus.GetMousePos(mouse_pos.x, mouse_pos.y);
	plus.GetMouseDt(mouse_dt.x, mouse_dt.y);

	auto size_window = plus.GetRenderWindowSize(plus.GetRenderWindow());
	if (plus.MouseButtonDown(hg::Button0))
		dens_prev[IX((mouse_pos.x / size_window.x) * N, (mouse_pos.y / size_window.y) * N)] = base_force_dens;
	if (plus.MouseButtonDown(hg::Button1))
	{
		u_prev[IX((mouse_pos.x / size_window.x) * N, (mouse_pos.y / size_window.y) * N)] = base_force_u * mouse_dt.x;
		v_prev[IX((mouse_pos.x / size_window.x) * N, (mouse_pos.y / size_window.y) * N)] = base_force_v * mouse_dt.y;
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

void draw_dens(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, const int &N, farray &dens, farray &u, farray &v) {
	float scale = 1.f;
	if (ui_draw_line)
		for (int i = 1; i <= N + 1; ++i)
			for (int j = 1; j <= N + 1; ++j) {
				hg::Vector3 p(i / float(N) - 0.5f, j / float(N) - 0.5f, 0);
				draw_line(simple_graphic_scene_overlay, p, p + hg::Vector3(u[IX(i, j)], v[IX(i, j)], 0));
			}

	for (int i = 1; i <= N + 1; ++i)
		for (int j = 1; j <= N + 1; ++j) {
			hg::Vector3 p(i / float(N) - 0.5f, j / float(N) - 0.5f, 0);
			draw_quad(simple_graphic_scene_overlay, hg::Matrix4::TransformationMatrix(p + hg::Vector3(0, 0, 0.1f), hg::Vector3(0, 0, 0)), scale / N, scale / N, nullptr, hg::Color(dens[IX(i, j)], 0, 0));
		}
}

void simulation_step(std::shared_ptr<hg::SimpleGraphicSceneOverlay> &simple_graphic_scene_overlay, const float &dt) {
	if (ImGui::Begin("params")) {
		ImGui::SliderFloat("visc", &visc, 0.0, 10.0);
		ImGui::SliderFloat("diff", &diff, 0.0, 1000.0);
		ImGui::SliderFloat("base_force_u", &base_force_u, 0.0, 1000.0);
		ImGui::SliderFloat("base_force_v", &base_force_v, 0.0, 1000.0);
		ImGui::SliderFloat("base_force_dens", &base_force_dens, 0.0, 1000.0);
		ImGui::Checkbox("Draw Vector Field", &ui_draw_line);
	}
	ImGui::End();

	get_from_UI(dens_prev, u_prev, v_prev);
	vel_step(N, u, v, u_prev, v_prev, visc, dt);
	dens_step(N, dens, dens_prev, u, v, diff, dt);
	draw_dens(simple_graphic_scene_overlay, N, dens, u, v);
}