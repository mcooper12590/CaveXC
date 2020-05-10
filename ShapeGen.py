from numpy import sin, cos, pi, sqrt, fabs, arctan2, linspace
import numpy as np

# Number of points that define a cross-section
d = 1000

# Generate x, y points for an ellipse. Can set rotation angle with theta
def genEll(r1,r2,theta=0):
	"""
	Generates an ellipise consisting of x, y points.

	Parameters
	----------
	r1 : x-axis radius
	r2 : y-axis radius
	theta : rotation angle

	Returns
	-------
	x : x-coordinates of ellipse
	y : y-coordinates of ellipse
	"""

	t=linspace(0, 2*pi-2*pi/d, d-1)
	x=r1*cos(t)
	y=r2*sin(t)


	if(theta!=0):

		tx=x
		ty=y

		x=tx*cos(theta)-ty*sin(theta)
		y=tx*sin(theta)+ty*cos(theta)

	return x, y

# Generate x, y points for a circle
def genCirc(r):
	"""
	Generates a circle.

	Parameters
	----------
	r : radius of circle

	Returns
	-------
	x : x-coordinates of circle
	y : y-coordinates of circle
	"""

	return genEll(r, r)

# Generate x, y points for the lower half of a circle
def genSemiCirc(r):
	"""
	Generates an open lower semi-circle

	Parameters
	----------
	r : radius of semi-circle

	Returns
	-------
	x : x-coordinates of semi-circle
	y : y-coordinates of semi-circle
	"""

	return genSemiEll(r, r)

# Generate x, y points for the lower half of an ellipse
def genSemiEll(r1,r2):
	"""
	Generates an open lower semi-ellipse

	Parameters
	----------
	r1 : half-width of semi-ellipse
	r2 : height of semi-ellipise

	Returns
	-------
	x : x-coordinates of semi-ellipse
	y : y-coordinates of semi-ellipse
	"""

	t=linspace(-1*pi, 0, d)
	x=r1*cos(t)
	y=r2*sin(t)

	return x, y
