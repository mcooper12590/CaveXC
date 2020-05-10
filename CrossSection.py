from numpy import sin, cos, pi, fabs, sign, roll, arctan2, diff, cumsum, hypot, logical_and, where, linspace
from scipy import interpolate
import matplotlib.pyplot as plt

# Cross-section class that stores x, y points for the cross-section
# and calculates various geometry data

d = 1000

class CrossSection:

	"""
	Class that contains geometry of and measurements of a cross-section.

	Parameters
	----------
	x : x-coordinates of cross-section
	y : y-coordinates of cross-section

	Attributes
	----------
	x : x-coordinates of cross-section
	y : y-coordinates of cross-section
	xm : coordinates x[i-1]
	xp : coordinates x[i+1]
	ym : coordinates y[i-1]
	yp : coordinates y[i+1]
	A : cross-sectional area
	P : cross-section perimeter
	l : distance between x[i-1], y[i-1] and x, y
	r_l : distance between x, y and reference point
	umx : x-coordinate of reference point
	umy : y-coordinate of reference point

	"""

	# Number of points that define the cross-section

	def __init__(self, x, y):

		self.x = x
		self.y = y
		self.roll()

	# Sets the point of maximum velocity
	def setUMPoint(self, umx, umy):
		"""
		Sets umx, umy.

		Parameters
		----------
		umx : x-coordinate of reference point
		umy : y-coordinate of reference point
		"""

		self.umx = umx
		self.umy = umy

	# Create arrays of x+1, y+1, x-1, x+1
	def roll(self):
		"""Creates xm, xp, ym, yp.
		"""

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
		"""
		Creates l and P.
		"""

		self.l = hypot(self.x - self.xp, self.y - self.yp)
		self.pp = cumsum(self.l)
		self.P = self.pp[-2]

	# Calculates area of the cross-section
	def calcA(self):
		"""
		Creates A.
		"""

		self.sA = (self.xm*self.y - self.x*self.ym).sum() * 0.5
		self.A = fabs(self.sA)

	# Generate lengths from maximum velocity point to perimeter points
	def genRL(self):
		"""
		Creates r_l.
		"""

		self.r_l = hypot(self.x-self.umx, self.y-self.umy)

	# Find left and right points defining a height above the cross-section
	# bottom
	def findLR(self, h):
		"""
		Finds left and right index given a height above the
		lowest point in the cross-section.

		Parameters
		----------
		h : height above the floor

		Returns
		-------
		L : left index of x, y coordinate h above the floor
		R : right index of x, y coordinate h above the floor

		"""
		ymin = self.y.min()
		a_h = ymin + h

		condL = logical_and(self.y > a_h, a_h > self.yp)
		condR = logical_and(self.y < a_h, a_h < self.yp)

		L = where(condL)[0][0] + 1
		R = where(condR)[0][0]
		return L,R

	# Find centroid, maximum velocity position in phreatic cases
	def findCentroid(self):
		"""
		Calculates centroid of the cross-section.

		Returns
		-------
		cx : x-coordinate of centroid
		cy : y-coordinate of centroid
		"""

		m = self.xm*self.y-self.x*self.ym
		cx = (1/(6*self.sA))*((self.x + self.xm)*m).sum()
		cy = (1/(6*self.sA))*((self.y + self.ym)*m).sum()

		return cx, cy

	# Redraw some length rl away normal to the perimeter
	# It may be advantageous for stability to resample using a spline fit
	# Setting dl sets the number of points defining the cross-section
	## after resampling.
	def redraw(self, rl, resample=False, dl=d):
		"""
		Regenerate cross-section perpendicular to current
		given a distance for each x,y point.

		Parameters
		----------
		rl : array of distances to move x, y points
		resample : [bool] option to resample points equidistantly along
		perimeter (optional)
		dl : number of points in resampled cross-section (optional)
		"""

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
	"""
	Determines if redrawn points are counter clockwise in cross-section

	Parameters
	----------
	x : x-coordinates of cross-section
	y : y-coordinates of cross-section
	xm : x[i-1]
	ym : y[i-1]
	nx : new x-coordinate
	ny : new y-coordinate

	Returns
	-------
	ccw : Array of bools indicating which new points are counter clockwise
	"""

	        return (x - xm) * (ny - ym) > (y - ym) * (nx - xm)

# Calculate length of curve defined by points
def calcL(x,y):
	"""
	Calculates length of a curve given x,y points

	Parameters
	----------
	x : x-coordinates of points
	y : y-coordinates of points

	Returns
	-------
	length : length of curve
	"""
	sub_lengths = hypot(x[1:] - x[:-1], y[1:] - y[:-1])
	sub_sums = cumsum(sub_lengths)
	length = sub_sums[-1]
	return length

def calcArea(x, y, l=0, r=0):
	"""
	Calculates area of a polygon given x,y points

	Parameters
	----------
	x : x-coordinates of points defining polygon
	y : y-coordinates of points defining polygon
	l : left index of subset of points (optional)
	r : right index of subset of points (optional)

	Returns
	-------
	A - area of polygon
	"""
	if l and r:
		sA = (roll(x[l:r],1)*y[l:r] - x[l:r]*roll(y[l:r],1)).sum()
	else:
		sA = (roll(x,1)*y - x*roll(y,1)).sum()
	return fabs(0.5 * sA)
