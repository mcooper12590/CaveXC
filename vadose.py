import sys
from numpy import zeros, linspace, flip, ones, concatenate, arange
from matplotlib.pyplot import figure, plot, axis, show
from ShapeGen import *
from CrossSection import *
from flow import *
from scipy.optimize import fmin

# Flow parameters

Q =60. #discharge, m^3/s
D_s = 0.001 #sediment diameter, m
S = 10**-3 #slope, unitless

rh = 6.8*D_s/30. #roughness height, m

g_h = 0.25 # initial free surface position guess, m

# Erosion parameters
a = 0.25 # exponent in erosion law
ff = 10 # Factor to adjust k in shear stress law

t = 10000+1 # total number of time steps

# Output file for saving data
sFile = "WTAWidths.csv"

# Create initial cross-section geometry
## This geometry is an elliptical channel bottom with vertical sides.
w = 8. # half initial width, m
d = 2 # initiail depth, m
x, y = genSemiEll(w,d)

y_s2 = linspace(0.01,5,50)
y_s1 = flip(y_s2)
x_s1 = ones(50)*-8.
x_s2 = ones(50)*8.

nx = concatenate((x_s1, x, x_s2))
ny = concatenate((y_s1, y, y_s2))

cs = CrossSection(nx, ny)

# Per timestep data
widths = zeros(t)

for ts in arange(t):
    FH = fmin(findH, g_h, (cs, Q, S, rh), disp=False) # find position of free-surface by minimizing difference between calcQ
    h = FH[0] # height of free-surface above bottom of cross-section
    g_h = h # use for next guess
    L, R = cs.findLR(h) # get indices for left and right of free-surface

    # New wetted cross-section from portion below free-surface
    wx = cs.x[L:R+1]
    wy = cs.y[L:R+1]
    wcs = CrossSection(wx, wy)
    widths[ts] = wx[-1]-wx[0]
    uMx = 0. # U_{max} ref point in center of f-s
    #uMx = (wcs.x[-1]+wcs.x[0])/2
    uMy = wcs.y[0]
    wcs.setUMPoint(uMx, uMy)
    wcs.calcShapeParams()
    wcs.genRL()

    F = flow(wcs, Q, rh) # flow object to find shear stress
    F.setS(S) # we define slope in f-s case
    F.calcFlowParams()

    # Adjust onto bed-normal vector
    phi = fabs( arctan2( (wcs.umy - wcs.y), (wcs.umx - wcs.x) ) )
    alpha = arctan2(wcs.yp-wcs.ym, wcs.xp-wcs.xm)
    alpha[0] = alpha[1]
    alpha[-1] = alpha[-2]
    F.regenT_b(phi, alpha) # adjust \tau_b with angles


    # Redraw cross-section with shear stress erosion law
    rl = zeros(cs.x.size)
    rl[L:R+1] = ff*( ( 1*10**(-1*(2*a+2)) )/50. )*F.T_b**(a) # power law erosion
    cs.redraw(rl, resample=True) # redraw cross-section w/ spline resample

    # Print width every 1k timesteps
    if ts%1000 == 0:
        print(widths[ts])

# Plot channel geometry
fi = figure(figsize=(10,6))
plot(cs.x,cs.y,'.')
axis("equal")
show()

# Ask if the channel is equilibrated
if sys.version_info[0]<3:
    yn = raw_input("Has channel reached a constant width [y/n]? ")
else:
    yn = input("Has channel reached a constant width [y/n]? ")


L, R = cs.findLR(h)
x_l = cs.x[L-15:L-5].mean()
x_r = cs.x[R+5:R+15].mean()

w_mean = x_r - x_l


if yn=="y":
    f = open(sFile, "a")
    f.write("%f,%f,%f,%f,%f\n" % (Q,a,w_mean, D_s, S))
    f.close()
