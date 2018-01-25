import numpy as np
import harfang as hg
import helper_2d
from harfang_shortcut import *

visc = 0.0
diff = 0.0
base_force_u, base_force_v, base_force_dens = 5, 5, 100

N = 64
size = (N + 2) * (N + 2)

uvdens_tex, uvdens_render_target, uvdens_tex_prev, uvdens_render_target_prev = None, None, None, None
idx, vtx, vtx_layout = None, None, None
shader, shader_set_inputs, shader_add_sources, shader_draw_uvdens, shader_lin_solve, shader_advect, shader_projectA, shader_projectB = None, None, None, None, None, None, None, None

plus = hg.GetPlus()


u = np.zeros((N + 2, N + 2))
u_prev = np.zeros((N + 2, N + 2))
v = np.zeros((N + 2, N + 2))
v_prev = np.zeros((N + 2, N + 2))
dens = np.zeros((N + 2, N + 2))
dens_prev = np.zeros((N + 2, N + 2))


def add_source(N, dt):	
	renderer = plus.GetRenderer()

	renderer.SetRenderTarget(uvdens_render_target)
	renderer.SetShader(shader_add_sources)
	renderer.SetShaderFloat(renderer.GetShaderVariable("dt"), dt)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_in"), uvdens_tex_prev)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_out"), uvdens_tex)
	hg.SetShaderEngineValues(plus.GetRenderSystem())
	hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)
	

def set_bnd(N, b, x):
	for i in range(1, N + 1):
		x[0 ,i] = x[1,i] * -1 if(b == 1) else x[1,i]
		x[N + 1,i] = x[N,i] * -1 if(b == 1) else x[N,i]
		x[i,0] = x[i,1] * -1 if(b == 2) else x[i,1]
		x[i,N + 1] = x[i,N] * -1 if(b == 2) else x[i,N]
	
	x[0 ,0] = 0.5 * (x[1,0] + x[0 ,1])
	x[0 ,N + 1] = 0.5 * (x[1,N + 1] + x[0 ,N])
	x[N + 1,0] = 0.5 * (x[N,0] + x[N + 1,1])
	x[N + 1,N + 1] = 0.5 * (x[N,N + 1] + x[N + 1,N])


def conv2d(a, f):
	s = f.shape + tuple(np.subtract(a.shape, f.shape) + 1)
	strd = np.lib.stride_tricks.as_strided
	subM = strd(a, shape = s, strides = a.strides * 2)
	return np.einsum('ij,ijkl->kl', f, subM)


def lin_solve(N, a, c):
	for k in range(20):
		#for i in range(1, N+1):
		#	for j in range(1, N+1):
		#		x[i, j] = (x0[i, j] + a * (x[i-1, j] + x[i+1, j] + x[i, j-1] + x[i, j+1]))/ c
		
		# double the frame rate, but it's not exactly the same, never find, close
		# enough
		#x[1:-1, 1:-1] = (x0[1:-1, 1:-1] + a * (conv2d(x, np.array([[0,1,0],[1,0,1],[0,1,0]])))) / c

		#set_bnd(N, b, x)

		renderer = plus.GetRenderer()

		renderer.SetRenderTarget(uvdens_render_target)
		renderer.SetShader(shader_lin_solve)
		renderer.SetShaderFloat2(renderer.GetShaderVariable("size_pixel"), 1.0/size, 1.0/size)
		renderer.SetShaderFloat(renderer.GetShaderVariable("a"), a)
		renderer.SetShaderFloat(renderer.GetShaderVariable("c"), c)
		renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_in"), uvdens_tex_prev)
		renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_out"), uvdens_tex)
		hg.SetShaderEngineValues(plus.GetRenderSystem())
		hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)


def diffuse(N, diff, dt):	
	a = dt * diff * N * N
	lin_solve(N, a, 1 + 4 * a)
	

def advect(N, b, d, d0, u, v, dt):
	dt0 = dt * N
	#for i in range(1, N + 1):
	#	for j in range(1, N + 1):
	#		x = i - dt0 * u[i, j]
	#		y = j - dt0 * v[i, j]
	#		if x < 0.5:
	#			x = 0.5 
	#		if x > N + 0.5:
	#		   x = N + 0.5
	#		i0 = int(x)
	#		i1 = i0 + 1
	#		if y < 0.5: 
	#			y = 0.5
	#		if y > N + 0.5:
	#		   y = N + 0.5 
	#		j0 = int(y)
	#		j1 = j0 + 1
	#		s1 = x - i0 
	#		s0 = 1 - s1 
	#		t1 = y - j0 
	#		t0 = 1 - t1
	#		d[i, j] = s0 * (t0 * d0[i0,j0] + t1 * d0[i0,j1]) + s1 * (t0 * d0[i1,j0] + t1 * d0[i1,j1])
			

	#set_bnd(N, b, d)
	

	renderer = plus.GetRenderer()

	renderer.SetRenderTarget(uvdens_render_target)
	renderer.SetShader(shader_advect)
	renderer.SetShaderFloat2(renderer.GetShaderVariable("size_pixel"), 1.0/size, 1.0/size)
	renderer.SetShaderFloat(renderer.GetShaderVariable("dt0"), dt0)
	renderer.SetShaderFloat(renderer.GetShaderVariable("N"), N)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_in"), uvdens_tex_prev)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_out"), uvdens_tex)
	hg.SetShaderEngineValues(plus.GetRenderSystem())
	hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)
	

def dens_step(N, x, x0, u, v, diff, dt):
	#add_source(N, x, x0, dt)
	x0, x = x, x0
	diffuse(N, 0, x, x0, diff, dt)
	x0, x = x, x0
	advect(N, 0, x, x0, u, v, dt)


def project(N, u, v, p, div):
	renderer = plus.GetRenderer()

	#for i in range(1, N + 1):
	#	for j in range(1, N + 1):
	#		div[i, j] = -0.5 * (u[i+1, j] - u[i-1, j] + v[i, j+1] - v[i, j-1]) / N
	#		p[i, j] = 0
	#div[1:-1, 1:-1] = -0.5 * (conv2d(u, np.array([[0,-1,0],[0,0,0],[0,1,0]])) + conv2d(v, np.array([[0,0,0],[-1,0,1],[0,0,0]]))) / N
	#p[:] = 0
		
	renderer.SetRenderTarget(uvdens_render_target_prev)
	renderer.SetShader(shader_projectA)
	renderer.SetShaderFloat2(renderer.GetShaderVariable("size_pixel"), 1.0/size, 1.0/size)
	renderer.SetShaderFloat(renderer.GetShaderVariable("N"), N)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_in"), uvdens_tex)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_out"), uvdens_tex_prev)
	hg.SetShaderEngineValues(plus.GetRenderSystem())
	hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)
	

	#set_bnd(N, 0, div)
	#set_bnd(N, 0, p)	
	
	lin_solve(N, 1, 4)

	#for i in range(1, N + 1):
	#	for j in range(1, N + 1):
	#		u[i, j] -= 0.5 * N * (p[i+1, j] - p[i-1, j])
	#		v[i, j] -= 0.5 * N * (p[i, j+1] - p[i, j-1])
	#v[1:-1, 1:-1] -= 0.5 * N * (conv2d(p, np.array([[0,0,0],[-1,0,1],[0,0,0]])))
	#u[1:-1, 1:-1] -= 0.5 * N * (conv2d(p, np.array([[0,-1,0],[0,0,0],[0,1,0]])))

	
	renderer.SetRenderTarget(uvdens_render_target)
	renderer.SetShader(shader_projectB)
	renderer.SetShaderFloat2(renderer.GetShaderVariable("size_pixel"), 1.0/size, 1.0/size)
	renderer.SetShaderFloat(renderer.GetShaderVariable("N"), N)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_in"), uvdens_tex_prev)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_out"), uvdens_tex)
	hg.SetShaderEngineValues(plus.GetRenderSystem())
	hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)

	#set_bnd(N, 1, u)
	#set_bnd(N, 2, v)

	
def vel_step(N, u, v, u0, v0, visc, dt):
	global uvdens_render_target_prev, uvdens_render_target, uvdens_tex_prev, uvdens_tex
	add_source(N, dt)
	
	uvdens_render_target_prev, uvdens_render_target = uvdens_render_target, uvdens_render_target_prev
	uvdens_tex_prev, uvdens_tex = uvdens_tex, uvdens_tex_prev

	diffuse(N, visc, dt)

	
	project(N, u, v, u0, v0)
	
	uvdens_render_target_prev, uvdens_render_target = uvdens_render_target, uvdens_render_target_prev
	uvdens_tex_prev, uvdens_tex = uvdens_tex, uvdens_tex_prev

	advect(N, 1, u, u0, u0, v0, dt) 	
	project(N, u, v, u0, v0)

	#add_source(N, u, u0, dt)
	#add_source(N, v, v0, dt)
	#u0, u = u, u0

	#diffuse(N, 1, u, u0, visc, dt)
	#v0, v = v, v0
	#diffuse(N, 2, v, v0, visc, dt)
	#project(N, u, v, u0, v0)
	#u0, u = u, u0
	#v0, v = v, v0
	#advect(N, 1, u, u0, u0, v0, dt) 
	#advect(N, 2, v, v0, u0, v0, dt)
	#project(N, u, v, u0, v0)

	


def get_from_UI(dens_prev, u_prev, v_prev):
	#dens_prev[:] = 0.0
	#u_prev[:] = 0.0
	#v_prev[:] = 0.0

	#if  hg.GetPlus().KeyDown(hg.KeyF1):
	#	dens_prev[1, 1] = base_force_dens
	#if  hg.GetPlus().KeyDown(hg.KeyF2):
	#	for i in range(2, 5):
	#		v_prev[i, 5] = base_force_v
	#if  hg.GetPlus().KeyDown(hg.KeyF3):
	#	u_prev[3, 3] = base_force_u
	#if  hg.GetPlus().KeyDown(hg.KeyF4):
	#	u_prev[int(N / 2.0), int(N / 2.0)] = base_force_u
	#	v_prev[int(N / 2.0), int(N / 2.0)] = base_force_v
	
	#for i in range(1, N+1):
	#	v_prev[IX(1, i)] = base_force_v if int(N/3) < i < int(N/3*2) else 0
	#	u_prev[IX(1, i)] = 0

	
	renderer = plus.GetRenderer()
	prev_viewport = renderer.GetViewport()

	renderer.SetRenderTarget(uvdens_render_target_prev)
	renderer.SetViewport(hg.Rect(0, 0, size, size))
	renderer.EnableDepthTest(False)
	renderer.Set2DMatrices()
	renderer.SetShader(shader_set_inputs)
	renderer.SetShaderFloat2(renderer.GetShaderVariable("size_pixel"), 1.0/size, 1.0/size)
	
	mouse_pos = plus.GetMousePos()
	renderer.SetShaderFloat3(renderer.GetShaderVariable("left_button_pressed"), mouse_pos[0]/prev_viewport.ex, 1.0-mouse_pos[1]/prev_viewport.ey, 1 if plus.MouseButtonDown(hg.Button0) else 0)
	hg.SetShaderEngineValues(plus.GetRenderSystem())
	hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)


false_texture = hg.GetPlus().LoadTexture("")

def draw_dens(simple_graphic_scene_overlay, N, dens, u, v):	
	#scale = 2
	#for i in range(1, N + 1):
	#	for j in range(1, N + 1):
	#		p = vec3(i / (N) - 0.5, j / (N) - 0.5, 0)
	#		#helper_2d.draw_line(simple_graphic_scene_overlay, p, p + vec3(u[i, j], v[i,j], 0))
	#		helper_2d.draw_quad(simple_graphic_scene_overlay, mat4.TransformationMatrix(p + vec3(0, 0, 0.1), vec3(0, 0, 0)), 1 / N, 1 / N, false_texture, col(dens[i, j], 0, 0))
		
	renderer = plus.GetRenderer()

	renderer.SetRenderTarget(uvdens_render_target)
	renderer.SetShader(shader_add_sources)
	renderer.SetShaderTexture(renderer.GetShaderVariable("tex_xydens_in"), uvdens_tex)
	hg.SetShaderEngineValues(plus.GetRenderSystem())
	hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)

	
def simulation_step(simple_graphic_scene_overlay, dt):
	global visc, diff, base_force_u, base_force_v, base_force_dens
	hg.ImGuiLock()
	if hg.ImGuiBegin("params"):
		visc = hg.ImGuiSliderFloat("visc", visc, 0.0, 10.0)[1]
		diff = hg.ImGuiSliderFloat("diff", diff, 0.0, 1000.0)[1]
		base_force_u = hg.ImGuiSliderFloat("base_force_u", base_force_u, 0.0, 1000.0)[1]		
		base_force_v = hg.ImGuiSliderFloat("base_force_v", base_force_v, 0.0, 1000.0)[1]		
		base_force_dens = hg.ImGuiSliderFloat("base_force_dens", base_force_dens, 0.0, 1000.0)[1]		
	hg.ImGuiEnd()
	hg.ImGuiUnlock()

	renderer = plus.GetRenderer()
	prev_viewport = renderer.GetViewport()

	get_from_UI(dens_prev, u_prev, v_prev)
	vel_step(N, u, v, u_prev, v_prev, visc, dt)
	#dens_step(N, dens, dens_prev, u, v, diff, dt)
	draw_dens(simple_graphic_scene_overlay, N, dens, u, v)

	# compute the shader

	#renderer.SetRenderTarget(uvdens_render_target_prev)
	#renderer.SetViewport(hg.Rect(0, 0, size, size))
	#renderer.Clear(col.Blue)
	#renderer.EnableDepthTest(False)
	#renderer.Set2DMatrices()
	#renderer.SetShader(shader)
	#renderer.SetShaderFloat2(renderer.GetShaderVariable("size_pixel"), 1.0/size, 1.0/size)
	#renderer.SetShaderTexture(renderer.GetShaderVariable("tex_in"), uvdens_tex)
	#hg.SetShaderEngineValues(plus.GetRenderSystem())
	#hg.DrawBuffers(renderer, 6, idx, vtx, vtx_layout)
	
	
	# reset matrices and draw the texture
	renderer.SetIdentityMatrices()
	renderer.SetViewport(prev_viewport)
	renderer.ClearRenderTarget()	
	helper_2d.draw_quad(simple_graphic_scene_overlay, mat4.Identity, 1, 1, uvdens_tex, col.White)



def setup():
	global uvdens_tex, uvdens_render_target, uvdens_tex_prev, uvdens_render_target_prev
	global idx, vtx, vtx_layout
	global shader, shader_set_inputs, shader_add_sources, shader_draw_uvdens, shader_lin_solve, shader_advect, shader_projectA, shader_projectB

	renderer = plus.GetRenderer()

	# create the 2 render target
	def create_tex_and_render_target():
		tex = renderer.NewTexture()
		pic = hg.Picture(size, size, hg.PictureRGBAF)
		renderer.CreateTexture(tex, pic, hg.TextureIsRenderTarget | hg.TextureIsShaderResource)

		render_target = renderer.NewRenderTarget()
		renderer.CreateRenderTarget(render_target)
		renderer.SetRenderTargetColorTexture(render_target, tex)

		return tex, render_target

	uvdens_tex, uvdens_render_target = create_tex_and_render_target()
	uvdens_tex_prev, uvdens_render_target_prev = create_tex_and_render_target()

	# Create index buffer
	data = hg.BinaryData()
	data.WriteInt16s([0, 1, 2, 0, 2, 3])
	idx = renderer.NewBuffer()
	renderer.CreateBuffer(idx, data, hg.GpuBufferIndex)

	# Create vertex buffer
	vtx_layout = hg.VertexLayout()
	vtx_layout.AddAttribute(hg.VertexPosition, 3, hg.VertexFloat)
	vtx_layout.AddAttribute(hg.VertexUV0, 2, hg.VertexUByte)

	data = hg.BinaryData()
	x, y = 1, 1
	data.WriteFloats([-x, -y, 0.5])
	data.WriteUInt8s([0, 0])
	data.WriteFloats([-x, y, 0.5])
	data.WriteUInt8s([0, 1])
	data.WriteFloats([x, y, 0.5])
	data.WriteUInt8s([1, 1])
	data.WriteFloats([x, -y, 0.5])
	data.WriteUInt8s([1, 0])

	vtx = renderer.NewBuffer()
	renderer.CreateBuffer(vtx, data, hg.GpuBufferVertex)

	# Load shader.
	shader = renderer.LoadShader("@assets/vel.isl")
	shader_set_inputs = renderer.LoadShader("@assets/set_inputs.isl")
	shader_add_sources = renderer.LoadShader("@assets/add_source.isl")
	shader_draw_uvdens = renderer.LoadShader("@assets/draw_uvdens.isl")
	shader_lin_solve = renderer.LoadShader("@assets/lin_solve.isl")
	shader_advect = renderer.LoadShader("@assets/advect.isl")
	shader_projectA = renderer.LoadShader("@assets/projectA.isl")
	shader_projectB = renderer.LoadShader("@assets/projectB.isl")
	
	