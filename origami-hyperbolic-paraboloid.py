#!/usr/bin/python
from math import radians, cos, sin, sqrt, asinh
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-n", "--sides", type=int, default=4, help="|")
parser.add_argument("-m", "--pleats", type=int, default=16, help="|")
parser.add_argument("-a", "--pleat_angle", type=float, default=45.0, help="|")
parser.add_argument("-l", "--length", type=float, default=0.5, help="|")
parser.add_argument("-i", "--integration_steps", type=int, default=20000, help="|")
parser.add_argument("-s", "--scale", type=float, default=10.0, help="|")
parser.add_argument("-o", "--output", type=str, help="|")
arguments = parser.parse_args()

def vector3_multiply_by_scalar(a, s):
	return [a[0] * s, a[1] * s, a[2] * s]

def vector3_cross(a, b):
	return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

def vector3_normalize(v):
	length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
	if length:
		return [(v[0] / length), (v[1] / length), (v[2] / length)]
	return [0.0, 0.0, 0.0]

if arguments.output == None:
	arguments.output = "hypar.obj"
elif not arguments.output.endswith(".obj"):
	arguments.output += ".obj"

# trying to open output file for writing
output_file = None
try:
	output_file = open(arguments.output, "w")
except Exception:
	print(f"Failed opening output file: {arguments.output}")
	quit()

n = arguments.sides
pcurves = []
for index in range(n):
	angle = radians(360.0 * (index / n))
	x, y = cos(angle), sin(angle)
	z = x * x - y * y
	pcurves.append([x, y, z])

m = arguments.pleats
a = arguments.pleat_angle
delta = 1.0 / arguments.integration_steps

x = arguments.length
arc_length = sqrt(4.0*x*x+1.0)*(x/2.0) + asinh(2.0*x)*(1.0/4.0)

w = arc_length / m
l = w / cos(radians(a))
h = sqrt(l * l - w * w)

pcurve_vertices = []
vertices = [[0.0, 0.0, 0.0]]
for pcurve in pcurves:
	pcurve_vertices.append([])
	current_pcurve_vertices = pcurve_vertices[-1]

	index, length, pleat_sign = 0, 0.0, 1
	while len(current_pcurve_vertices) < m:
		t = index * delta
		dx = pcurve[0] * delta
		dy = pcurve[1] * delta
		dz = pcurve[2] * delta * (delta + 2 * t)
		length += sqrt(dx * dx + dy * dy + dz * dz)
		x = [pcurve[0] * t, pcurve[1] * t, pcurve[2] * t * t]
		normal = vector3_normalize(vector3_cross([dx, dy, dz], [-pcurve[1], pcurve[0], 0.0]))
		normal = vector3_multiply_by_scalar(normal, h * (0 if pleat_sign < 0 else 1))

		if length > w:
			x = pcurve[0] * t + normal[0]
			y = pcurve[1] * t + normal[1]
			z = pcurve[2] * t * t + normal[2]
			vertices.append([x, y, z])
			current_pcurve_vertices.append(len(vertices))
			pleat_sign *= -1
			length -= w
		index += 1

def write_obj():
	for v in vertices:
		ov = vector3_multiply_by_scalar(v, arguments.scale)
		output_file.write(f"v {ov[0]:g} {ov[2]:g} {ov[1]:g}\n")
	output_file.write("o 1\ns 1\n")
	for i2 in range(n):
		p1 = 1
		p2 = pcurve_vertices[i2][0]
		p3 = pcurve_vertices[(i2 + 1) % n][0]
		output_file.write(f"f {p1} {p2} {p3}\n")
	for i1 in range(m - 1):
		output_file.write(f"o {i1 + 2}\ns 1\n")
		for i2 in range(n):
			p1 = pcurve_vertices[i2][i1]
			p2 = pcurve_vertices[i2][i1 + 1]
			p3 = pcurve_vertices[(i2 + 1) % n][i1 + 1]
			p4 = pcurve_vertices[(i2 + 1) % n][i1]
			output_file.write(f"f {p1} {p2} {p3} {p4}\n")

def write_split_obj():
	for v in vertices:
		ov = vector3_multiply_by_scalar(v, arguments.scale)
		output_file.write(f"v {ov[0]:g} {ov[2]:g} {ov[1]:g}\n")

	objects = 0
	for i2 in range(n):
		objects += 1
		output_file.write(f"o {objects}\n")

		p1 = 1
		p2 = pcurve_vertices[i2][0]
		p3 = pcurve_vertices[(i2 + 1) % n][0]
		output_file.write(f"f {p1} {p2} {p3}\n")

		p1 = pcurve_vertices[i2][0]
		p2 = pcurve_vertices[i2][1]
		p3 = pcurve_vertices[(i2 + 1) % n][1]
		p4 = pcurve_vertices[(i2 + 1) % n][0]
		output_file.write(f"f {p1} {p2} {p3} {p4}\n")

	for i1 in range(1, m - 1, 2):
		for i2 in range(n):
			objects += 1
			output_file.write(f"o {objects}\n")

			p1 = pcurve_vertices[i2][i1]
			p2 = pcurve_vertices[i2][i1 + 1]
			p3 = pcurve_vertices[(i2 + 1) % n][i1 + 1]
			p4 = pcurve_vertices[(i2 + 1) % n][i1]
			output_file.write(f"f {p1} {p2} {p3} {p4}\n")

			p1 = pcurve_vertices[i2][i1 + 1]
			p2 = pcurve_vertices[i2][i1 + 2]
			p3 = pcurve_vertices[(i2 + 1) % n][i1 + 2]
			p4 = pcurve_vertices[(i2 + 1) % n][i1 + 1]
			output_file.write(f"f {p1} {p2} {p3} {p4}\n")

write_obj()

# closing output file
output_file.close()
