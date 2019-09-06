from numpy import log, sin, cos, pi, sqrt, roll, fabs, hypot
from CrossSection import *
import matplotlib.pyplot as plt


# Calculate the difference between prescribed Q_w and that determined by
# the height of the free-surface. This function is minimized to determine
# the free-surface location in such cases
def findH(g_h, cs, Q, S, rh):

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

		C = 2.5*sqrt(9.81)*log(0.37*self.R/(c_rh if c_rh else self.rh))
		self.S = self.u_bar**2/(C**2*self.R)

	# For free surface hydraulic gradient must be set
	def setS(self, S):
		self.S = S

	# Find the value of U_max by weighted law of the wall
	def calcUmax(self, method="area"):

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

		self.vgrad = (self.umax/self.rh)*(1/log(self.cs.r_l/self.rh))


	# Calculate phi from Wobus et al [2006]
	def calcphi(self):

		sum = ((self.vgrad**2)*self.cs.l).sum()
		self.phi = (9.81*self.S)/sum

	# Generate boundary shear stress around the perimeter
	def genT_b(self):

		self.T_b = self.phi*998.2*self.cs.A*self.vgrad**2

	# Calculate gradients and T_b with sine term
	def regenT_b(self, phi, alpha):

		nvgrad = self.vgrad*fabs(sin(phi-alpha))
		self.vgrad = nvgrad
		self.calcphi()
		self.genT_b()
