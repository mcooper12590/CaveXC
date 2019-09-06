from simps2 import simps2
from numpy import log, tan, radians, any

# Functions relating to erosion by saltating bed load, Nelson&Seminara, 2011

# Calculate sediment transport capacity for each point around perimeter
def calcLSTC(Tb, rho_s, D_s, T_sc):

	R_b = rho_s/999.97 - 1
	T_s = Tb/((rho_s - 999.97)*9.81*D_s)

	# If shear stress is under critical shields stress force alluviation
	if any((T_s-T_sc)<0):
		return -1
	else:
		q_t = 5.7*((T_s - T_sc)**(3./2))*2650*((R_b*9.81*D_s**3)**(1./2))
		return q_t


# Find the left and right coords. for bed-load layer
def findBedLoad(wcs, y_min, q_t, Q_s):

	# Alluviation forcing if critical shields stress is not reached
	if (type(q_t) is int): return 0, -1
        # Force alluviation if we reach the ceiling
    	elif simps2(q_t, wcs.x) < Q_s: return 0, -1

	tQ_s = 0.
	low = y_min
	i=0

	# Integrate sediment flux until it equals sed. supply
	while True:

		xRs = low + i + 1
		xLs = low - i
		xRsn = xRs + 1
		xLsn = xLs - 1

		# Integrate using simpson's method
		tQ_s = simps2(q_t[xLs:xRs], wcs.x[xLs:xRs])
		tQ_sn = simps2(q_t[xLsn:xRsn], wcs.x[xLsn:xRsn])

		if (tQ_s < Q_s and Q_s < tQ_sn):
			return xLs, xRs - 1

		i+=1

# Calculate critical shields stress from sed. diameter
def calcTscr(rho_s, D_s):

	if D_s < 0.0005: phi = 30.
	elif D_s >= 0.0005 and D_s < 0.064:
		phi = 2.16404*log(2.64225e9*D_s) #Fit from data
	else: phi = 42.

	if D_s < 0.002:
		ssg = (rho_s-999.97)/999.97
		d_star = D_s*(((ssg-1)*9.81)/((1e-6)**2))**(1./3)
		T_cr = 0.25*d_star**(-0.6)*9.81*(rho_s-999.97)*D_s * tan(radians(phi))

	else: T_cr = 0.06*9.81*(rho_s - 999.97)*D_s * tan(radians(phi))

	T_sc = T_cr/((rho_s - 999.97)*9.81*D_s)

	return T_sc
