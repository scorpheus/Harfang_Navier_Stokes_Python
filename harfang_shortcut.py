import harfang as hg
from math import radians, pow

ivec2 = hg.IntVector2
vec2 = hg.Vector2
vec3 = hg.Vector3
mat3 = hg.Matrix3
mat4 = hg.Matrix4
col = hg.Color

def vec3rad(x,y,z):
	return vec3(radians(x), radians(y), radians(z))

# The infamous clamp function 
def Clamp(k, a, b):
	k = min(k, b)
	k = max(k, a)
	return k

# Range : [0.0, 1.0]
def	EaseInOutQuick(x):
	x = Clamp(x, 0.0, 1.0)
	return	(x * x * (3 - 2 * x))

# Range : [0.0, 1.0]
def EaseInOutByPow(x, p = 2.0):
	x = Clamp(x, 0.0, 1.0)
	y = pow(x, p) / (pow(x, p) + pow(1 - x, p))
	return y


def create_instance(scn, path):
	node = hg.Node()
	node.AddComponent(hg.Transform())
	instance = hg.Instance()
	instance.SetPath(path)
	node.AddComponent(instance)
	node.GetTransform().SetPosition(vec3(0, 0, 0))
	scn.AddNode(node)
	return node