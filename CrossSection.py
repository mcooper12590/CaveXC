from numpy import sin, cos, pi, fabs, sign, roll, arctan2, diff, cumsum, hypot, logical_and, where, linspace
from scipy import interpolate
import matplotlib.pyplot as plt

# Cross-section class that stores x, y points for the cross-section
# and calculates various geometry data

d = 1000

class CrossSection:

	# Number of points that define the cross-section

	def __init__(self, x, y):

		self.x = x
		self.y = y
		self.roll()

	# Sets the point of maximum velocity
	def setUMPoint(self, umx, umy):

		self.umx = umx
		self.umy = umy

	# Create arrays of x+1, y+1, x-1, x+1
	def roll(self):

		self.xm = roll(self.x, 1)
		self.ym = roll(self.y, 1)
		self.xp = roll(self.x, self.x.size-1)
		self.yp = roll(self.y, self.y.size-1)

	# Calculate perimeter and area
	def calcShapeParams(self):

		self.genL()
		self.calcA()

	# l stores difference between each perimeter point
	# pp is the length along perimeter to a point
	# pp[-2] is the channel perimeter
	def genL(self):

		self.l = hypot(self.x - self.xp, self.y - self.yp)
		self.pp = cumsum(self.l)
		self.P = self.pp[-2]

	# Calculates area of the cross-section
	def calcA(self):

		self.sA = (self.xm*self.y - self.x*self.ym).sum() * 0.5
		self.A = fabs(self.sA)

	# Generate lengths from maximum velocity point to perimeter points
	def genRL(self):

		self.r_l = hypot(self.x-self.umx, self.y-self.umy)

	# Find left and right points defining a height above the cross-section
	# bottom
	def findLR(self, h):

		ymin = self.y.min()
		a_h = ymin + h

		condL = logical_and(self.y > a_h, a_h > self.yp)
		condR = logical_and(self.y < a_h, a_h < self.yp)

		L = where(condL)[0][0] + 1
		R = where(condR)[0][0]
		return L,R

	# Find centroid, maximum velocity position in phreatic cases
	def findCentroid(self):

		m = self.xm*self.y-self.x*self.ym
		cx = (1/(6*self.sA))*((self.x + self.xm)*m).sum()
		cy = (1/(6*self.sA))*((self.y + self.ym)*m).sum()

		return cx, cy

	# Redraw some length rl away normal to the perimeter
	# It may be advantageous for stability to resample using a spline fit
	# Setting dl sets the number of points defining the cross-section
	## after resampling.
	def redraw(self, rl, resample=False, dl=d):

		alpha = arctan2(self.xp-self.xm, self.yp-self.ym)

		nx = self.x + sign(self.x)*rl*cos(alpha)
		ny = self.y - sign(self.x)*rl*sin(alpha)

		# Check if we drew inside or outside..
		c = ccw(self.x, self.y, self.xm, self.ym, nx, ny)

		nx[c] = (self.x - sign(self.x)*rl*cos(alpha))[c]
		ny[c] = (self.y + sign(self.x)*rl*sin(alpha))[c]

		#Resample points by fitting spline
		if resample:
			tck, u = interpolate.splprep([nx, ny], u=None, k=1, s=0.0)
			un = linspace(u.min(), u.max(), dl if dl!=nx.size else nx.size)
			nx, ny = interpolate.splev(un, tck, der=0)


		# New coordinates
		y_roll = ny.size - ny.argmax()
		nx = roll(nx, y_roll)
		ny = roll(ny, y_roll)
		self.x = nx
		self.y = ny
		self.roll()


# Counter clockwise function to determine if we drew points in the correct
# direction
def ccw(x, y, xm, ym, nx, ny):

	        return (x - xm) * (ny - ym) > (y - ym) * (nx - xm)

# Calculate length of curve defined by points
def calcL(x,y):
	sub_lengths = hypot(x[1:] - x[:-1], y[1:] - y[:-1])
	sub_sums = cumsum(sub_lengths)
	length = sub_sums[-1]
	return length

def calcArea(x, y, l=0, r=0):
	if l and r:
		sA = (roll(x[l:r],1)*y[l:r] - x[l:r]*roll(y[l:r],1)).sum()
	else:
		sA = (roll(x,1)*y - x*roll(y,1)).sum()
	return fabs(0.5 * sA)
