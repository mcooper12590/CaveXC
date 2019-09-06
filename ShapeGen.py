from numpy import sin, cos, pi, sqrt, fabs, arctan2, linspace
import numpy as np

# Number of points that define a cross-section
d = 1000

# Generate x, y points for an ellipse. Can set rotation angle with theta
def genEll(r1,r2,theta=0):

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

	return genEll(r, r)

# Generate x, y points for the lower half of a circle
def genSemiCirc(r):

	return genSemiEll(r, r)

# Generate x, y points for the lower half of an ellipse
def genSemiEll(r1,r2):

	t=linspace(-1*pi, 0, d)
	x=r1*cos(t)
	y=r2*sin(t)

	return x, y

