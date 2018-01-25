import harfang as hg
import math
import bspline
from harfang_shortcut import *

font = hg.RasterFont("@assets/fonts/DS-DIGI.TTF", 16)


def get_spline_val(t, v1, v2, v3, v4):
	v = [v1, v2, v3, v4.z]

	C = bspline.C_factory(v, 1, "clamped")
	if C:
		val = C(t * C.max)
		return val
	return v2

def get_spline_val_vec(t, p1, p2, p3, p4):
	P = [(p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z), (p3.x, p3.y, p3.z), (p4.x, p4.y, p4.z)]

	C = bspline.C_factory(P, 3, "clamped")
	if C:
		val = C(t * C.max)
		return vec3(val[0], val[1], val[2])
	return p2


def draw_spline(scene_simple_graphic, p1, p2, p3, p4, color=hg.Color.White):
	P = [(p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z), (p3.x, p3.y, p3.z), (p4.x, p4.y, p4.z)]

	C = bspline.C_factory(P, 3, "clamped")
	if C:
		step = 50
		prev_value = [p1.x, p1.y, p1.z]
		for i in range(step):
			val = C(float(i)/step * C.max)
			scene_simple_graphic.Line(prev_value[0], prev_value[1], prev_value[2], val[0], val[1], val[2], color, color)
			prev_value = val


def kmh_to_mtrs(v):
	return v / 3.6


def rangeadjust(k, a, b, u, v):
	if u == v or b == a:
		return 0
	return (k - a) / (b - a) * (v - u) + u


def clamp(v, v_min, v_max):
	return max(v_min, min(v, v_max))


def rangeadjust_clamp(k, a, b, u, v):
	return rangeadjust(clamp(k, a, b), a, b, u, v)


def lerp(k, a, b):
	return a + (b - a) * k

# Range : [0.0, 1.0]
def	EaseInOutQuick(x):
	x = clamp(x, 0.0, 1.0)
	return	(x * x * (3 - 2 * x))

# Range : [0.0, 1.0]
def EaseInOutByPow(x, p = 2.0):
	x = clamp(x, 0.0, 1.0)
	y = pow(x, p) / (pow(x, p) + pow(1 - x, p))
	return y


def get_face_matrix(position, cam_pos):
		vec_look_at = (position - cam_pos)
		vec_look_at = vec_look_at.Normalized()
		return mat4.TransformationMatrix(position, hg.Matrix3.LookAt(vec_look_at))


def get_poly_from_minmax(minmax):
	return [vec3(minmax.GetMax(hg.AxisX), 0, minmax.GetMax(hg.AxisZ)),
			vec3(minmax.GetMax(hg.AxisX), 0, minmax.GetMin(hg.AxisZ)),
			vec3(minmax.GetMin(hg.AxisX), 0, minmax.GetMin(hg.AxisZ)),
			vec3(minmax.GetMin(hg.AxisX), 0, minmax.GetMax(hg.AxisZ))]


def get_poly_from_obb(obb):
	return [obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,
		 obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 - obb.rotation.GetZ() * obb.scale.z * 0.5,
		 obb.position - obb.rotation.GetX() * obb.scale.x * 0.5 - obb.rotation.GetZ() * obb.scale.z * 0.5,
		 obb.position - obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5]


def get_cube_from_obb(obb):
	return [
		obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,
		obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
		obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
		obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,

		obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,
		obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
		obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
		obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5
	 ]


def point_project_to_line(p, a, b, clamp_val=1.0):
	ap = p-a
	ab = b-a

	l = (ap.Dot(ab)/ab.Dot(ab))
	if clamp_val == -1:
		return a + ab * l, l
	else:
		return a + ab * clamp(l, 0.0, clamp_val), l


def point_in_poly_2d(point, poly):
	odd_nodes = False
	x2 = poly[3].x
	z2 = poly[3].z

	for p in poly:
		# vertex a
		x1 = p.x
		z1 = p.z
		if ((z1 < point.z) and (z2 >= point.z)) or (z1 >= point.z) and (z2 < point.z):
			if (point.z - z1) / (z2 - z1) * (x2 - x1) < (point.x - x1):
				odd_nodes = not odd_nodes

		x2 = x1
		z2 = z1

	return odd_nodes


def overlap_obb_2d(obb1, obb2):
	# test the sat for the 3 axis of the first obb
	# test for the cross product between the 3 first axis with the 3 other axis
	obb1_poly = get_cube_from_obb(obb1)
	obb2_poly = get_cube_from_obb(obb2)

	# 2D
	if hg.TestOverlap(obb1.rotation.GetX(), obb1_poly, obb2_poly) and \
		hg.TestOverlap(obb1.rotation.GetZ(), obb1_poly, obb2_poly) and \
		hg.TestOverlap(obb1.rotation.GetX().Cross(obb2.rotation.GetX()), obb1_poly, obb2_poly) and \
		hg.TestOverlap(obb1.rotation.GetX().Cross(obb2.rotation.GetZ()), obb1_poly, obb2_poly) and \
		hg.TestOverlap(obb1.rotation.GetZ().Cross(obb2.rotation.GetX()), obb1_poly, obb2_poly) and \
		hg.TestOverlap(obb1.rotation.GetZ().Cross(obb2.rotation.GetZ()), obb1_poly, obb2_poly):
		return True
	# 3D
	# if hg.TestOverlap(obb1.rotation.GetX(), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetY(), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetZ(), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetX().Cross(obb2.rotation.GetX()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetX().Cross(obb2.rotation.GetY()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetX().Cross(obb2.rotation.GetZ()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetY().Cross(obb2.rotation.GetX()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetY().Cross(obb2.rotation.GetY()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetY().Cross(obb2.rotation.GetZ()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetZ().Cross(obb2.rotation.GetX()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetZ().Cross(obb2.rotation.GetY()), obb1_poly, obb2_poly) and \
	# 	hg.TestOverlap(obb1.rotation.GetZ().Cross(obb2.rotation.GetZ()), obb1_poly, obb2_poly):
	# 	return True

	return False


def overlap_circles_2d(circles1, circles2):
	for id1, c1 in enumerate(circles1):
		for id2, c2 in enumerate(circles2):
			if hg.Dist2(c1["p"], c2["p"]) < (c1["r"] + c2["r"])**2:
				return True, id1, id2
	return False, None, None


def overlap_min_max_2d(minmax1, minmax2):
	return point_in_poly_2d(minmax1.mn, get_poly_from_minmax(minmax2)) or point_in_poly_2d(minmax1.mx, get_poly_from_minmax(minmax2)) or \
		   point_in_poly_2d(vec3(minmax1.mn.x, 0, minmax1.mx.z), get_poly_from_minmax(minmax2)) or point_in_poly_2d(vec3(minmax1.mx.x, 0, minmax1.mn.z), get_poly_from_minmax(minmax2))


def draw_minmax(scene_simple_graphic, minmax, color=hg.Color.White):
	scene_simple_graphic.Line(minmax.mn.x, minmax.mn.y, minmax.mn.z, minmax.mx.x, minmax.mn.y, minmax.mn.z, color, color)
	scene_simple_graphic.Line(minmax.mn.x, minmax.mn.y, minmax.mn.z, minmax.mn.x, minmax.mx.y, minmax.mn.z, color, color)
	scene_simple_graphic.Line(minmax.mn.x, minmax.mn.y, minmax.mn.z, minmax.mn.x, minmax.mn.y, minmax.mx.z, color, color)

	scene_simple_graphic.Line(minmax.mx.x, minmax.mx.y, minmax.mx.z, minmax.mn.x, minmax.mx.y, minmax.mx.z, color, color)
	scene_simple_graphic.Line(minmax.mx.x, minmax.mx.y, minmax.mx.z, minmax.mx.x, minmax.mn.y, minmax.mx.z, color, color)
	scene_simple_graphic.Line(minmax.mx.x, minmax.mx.y, minmax.mx.z, minmax.mx.x, minmax.mx.y, minmax.mn.z, color, color)

	scene_simple_graphic.Line(minmax.mn.x, minmax.mx.y, minmax.mn.z, minmax.mx.x, minmax.mx.y, minmax.mn.z, color, color)
	scene_simple_graphic.Line(minmax.mn.x, minmax.mx.y, minmax.mn.z, minmax.mn.x, minmax.mx.y, minmax.mx.z, color, color)

	scene_simple_graphic.Line(minmax.mx.x, minmax.mn.y, minmax.mx.z, minmax.mx.x, minmax.mn.y, minmax.mn.z, color, color)
	scene_simple_graphic.Line(minmax.mx.x, minmax.mn.y, minmax.mx.z, minmax.mn.x, minmax.mn.y, minmax.mx.z, color, color)


def draw_obb(scene_simple_graphic, obb, color=hg.Color.White):
	cube = [
	obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,
	 obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
	 obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
	 obb.position + obb.rotation.GetX() * obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,

	 obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5,
	 obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
	 obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * -obb.scale.z * 0.5,
	 obb.position + obb.rotation.GetX() * -obb.scale.x * 0.5 + obb.rotation.GetY() * -obb.scale.y * 0.5 + obb.rotation.GetZ() * obb.scale.z * 0.5
	 ]
	draw_cube(scene_simple_graphic, cube, color)


def draw_cube(scene_simple_graphic, cube, color=hg.Color.White):
	draw_line(scene_simple_graphic, cube[0], cube[1], color)
	draw_line(scene_simple_graphic, cube[1], cube[2], color)
	draw_line(scene_simple_graphic, cube[2], cube[3], color)
	draw_line(scene_simple_graphic, cube[0], cube[3], color)

	draw_line(scene_simple_graphic, cube[4], cube[5], color)
	draw_line(scene_simple_graphic, cube[5], cube[6], color)
	draw_line(scene_simple_graphic, cube[6], cube[7], color)
	draw_line(scene_simple_graphic, cube[4], cube[7], color)

	draw_line(scene_simple_graphic, cube[0], cube[4], color)
	draw_line(scene_simple_graphic, cube[1], cube[5], color)
	draw_line(scene_simple_graphic, cube[2], cube[6], color)
	draw_line(scene_simple_graphic, cube[3], cube[7], color)


def draw_line(scene_simple_graphic, a, b, color=hg.Color.White):
	scene_simple_graphic.Line(a.x, a.y, a.z, b.x, b.y, b.z, color, color)


def draw_cross(scene_simple_graphic, pos, color=hg.Color.White, size=0.5):
	scene_simple_graphic.Line(pos.x-size, pos.y, pos.z, pos.x+size, pos.y, pos.z, color, color)
	scene_simple_graphic.Line(pos.x, pos.y-size, pos.z, pos.x, pos.y+size, pos.z, color, color)
	scene_simple_graphic.Line(pos.x, pos.y, pos.z-size, pos.x, pos.y, pos.z+size, color, color)


def draw_circle(scene_simple_graphic, world, r, color=hg.Color.White):
	step = 50
	prev = vec3(math.cos(0) * r, 0, math.sin(0) * r) * world
	for i in range(step+1):
		val = vec3(math.cos(math.pi*2*float(i)/step) * r, 0, math.sin(math.pi*2*float(i)/step) * r) * world

		scene_simple_graphic.Line(prev.x, prev.y, prev.z, val.x, val.y, val.z, color, color)
		prev = val


def draw_text(scene_simple_graphic, text, mat, size=0.01, color=hg.Color.White, text_centered=False, use_font=font):
	scene_simple_graphic.SetDepthWrite(False)
	scene_simple_graphic.SetBlendMode(hg.BlendAlpha)
	if text_centered:
		text_rect = use_font.GetTextRect(hg.GetPlus().GetRenderSystem(), text)
		mat = mat * mat4.TranslationMatrix(vec3(-(text_rect.ex - text_rect.sx) * 0.5 * size, -(text_rect.ey - text_rect.sy) * 0.5 * size, 0))
	scene_simple_graphic.Text(mat, text, color, use_font, size)
	scene_simple_graphic.SetDepthWrite(True)


def draw_quad(scene_simple_graphic, mat, width, height, texture, color = hg.Color.White):
	pos = mat.GetTranslation()
	axis = [mat.GetX() * width / 2, mat.GetY() * height / 2, mat.GetZ()]
	a = pos - axis[0] - axis[1]
	b = pos - axis[0] + axis[1]
	c = pos + axis[0] + axis[1]
	d = pos + axis[0] - axis[1]

	scene_simple_graphic.Quad(a.x, a.y, a.z,
				b.x, b.y, b.z,
				c.x, c.y, c.z,
				d.x, d.y, d.z,
				0, 0, 1, 1,
				texture,
				color, color, color, color
				)


def draw_geometry(scene_simple_graphic, mat, geo):
	p, s, r = mat.Decompose(hg.RotationOrder_Default)
	scene_simple_graphic.Geometry(p.x, p.y, p.z, r.x, r.y, r.z, s.x, s.y, s.z, geo)