from numpy import log, sin, cos, pi, sqrt, roll, fabs, hypot
from CrossSection import *
import matplotlib.pyplot as plt


# Calculate the difference between prescribed Q_w and that determined by
# the height of the free-surface. This function is minimized to determine
# the free-surface location in such cases
def findH(g_h, cs, Q, S, rh):
	"""
	Calculates discharge difference for a given height of water.
	The minimum of this function is the height of water with a prescribed
	discharge.

	Parameters
	----------
	g_h : guess height [m]
	cs : cross-section
	Q : discharge [m^3/s]
	S : slope
	rh : roughness height (z_0) [m]

	Returns
	-------
	epsilon : difference between prescribed discharge and discharge for a
	given height.
	"""

	L, R = cs.findLR(g_h)
	wx = cs.x[L:R]
	wy = cs.y[L:R]
	wcs = CrossSection(wx, wy)
	wcs.calcShapeParams()
	Pw = wcs.pp[-1]
	A = wcs.A
	H = A/Pw

	C = 2.5*sqrt(9.81)*log(0.37*H/rh)

	u_bar = C*sqrt(H*S)
	Qcalc = u_bar*A

	return fabs(Q - Qcalc)

# Flow class for hydraulic info for wetted cross-section
class flow:
	"""
	Flow object for calculating boundary shear stress and other
	flow parameters.

	Parameters
	----------
	cs : cross-section object containing wetted cross-section
	Q : prescribed discharge [m^3/s]
	rh : roughness height (z_0) [m]

	Attributes
	----------
	cs : wetted cross-section object
	Q : prescribed discharge [m^3/s]
	rh : roughness height (z_0) [m]
	u_bar : average flow velocity [m/s]
	R : hydraulic radius of wetted cross-section [m]
	S : hydraulic gradient/energy slope
	umax : maximum flow velocity [m/s]
	T_b : boundary shear stress [Pa]

	Notes
	-----
	Implements method of calculating boundary shear stress for a
	cross-section from Wobus et al. [2006; 2008]
	"""

	def __init__(self, cs, Q, rh):

		self.cs = cs

		self.Q = Q
		self.rh = rh

		self.u_bar = Q/cs.A
		self.R = cs.A/cs.P


	# Calculate everything needed
	def calcFlowParams(self):

		self.calcUmax()
		self.genVGrad()
		self.calcphi()

		self.genT_b()

	# Calculate hydraulic gradient for phreatic case
	def calcChezy(self, c_rh=False):
		"""
		Calculate hydraulic gradient.

		Parameters
		----------
		c_rh : composite roughness (optional, given when per-point roughness
		is used)
		"""

		C = 2.5*sqrt(9.81)*log(0.37*self.R/(c_rh if c_rh else self.rh))
		self.S = self.u_bar**2/(C**2*self.R)

	# For free surface hydraulic gradient must be set
	def setS(self, S):
		"""
		Set hydraulic gradient to a prescribed value.

		Parameters
		----------
		S : hydraulic gradient/energy slope
		"""
		self.S = S

	# Find the value of U_max by weighted law of the wall
	def calcUmax(self, method="area"):
		"""
		Calculate maximum velocity at a given position.

		Parameters
		----------
		method : which integral method to use when calculating maximum
		velocity

		Notes
		-----
		Full details of algorithm in Supplemental 1 of
		Cooper&Covington, "Modeling cave cross-section evolution
		including sediment transport and paragenesis".
		"""

		R = self.cs.r_l
		z0 = rh = self.rh
		Ax = self.cs.umx
		Ay = self.cs.umy
		Bx = roll(self.cs.x,1)
		By = roll(self.cs.y,1)
		Cx = self.cs.x
		Cy = self.cs.y

		a_top = Ax*(By-Cy) + Bx*(Cy-Ay) + Cx*(Ay-By)
		a_i = fabs(a_top/2)

		if method=="line":
			u_i = R/(R-z0) - 1/log(R/z0)

		if method=="area":
			u_top = R**2*log(R/z0) - 3./2*R**2 + 2.*R*z0 - z0**2/2.
			u_bottom = z0**2 + R**2 - 2*z0*R
			u_i = 1./log(R/z0)*u_top/u_bottom

		self.umax = self.Q/(a_i*u_i).sum()

	# Generate velocity gradients from Wobus et al [2006], sans sine term
	def genVGrad(self):
		"""
		Calculates velocity gradient for each point along perimeter.

		Notes
		-----
		Equation 3 in Wobus et al [2006], "Self-formed bedrock channels",
		sans sine term.
		"""

		self.vgrad = (self.umax/self.rh)*(1/log(self.cs.r_l/self.rh))


	# Calculate phi from Wobus et al [2006]
	def calcphi(self):
		"""
		Calculates force balance correction.

		Notes
		-----
		Equation 5 in Wobus et al [2006], "Self-formed bedrock channels".
		"""

		sum = ((self.vgrad**2)*self.cs.l).sum()
		self.phi = (9.81*self.S)/sum

	# Generate boundary shear stress around the perimeter
	def genT_b(self):
		"""
		Calculates boundary shear stress.

		Notes:
		Equation 4 in Wobus et al [2006], "Self-formed bedrock channels".
		"""

		self.T_b = self.phi*998.2*self.cs.A*self.vgrad**2

	# Calculate gradients and T_b with sine term
	def regenT_b(self, phi, alpha):
		"""
		Calculates T_b on bed-normal vector.

		Parameters
		----------
		phi : angle of line connecting perimeter points to reference point
		alpha : angle of line tangent to perimeter points
		"""

		nvgrad = self.vgrad*fabs(sin(phi-alpha))
		self.vgrad = nvgrad
		self.calcphi()
		self.genT_b()
