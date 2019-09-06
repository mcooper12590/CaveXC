from numpy import zeros, trim_zeros, ones, linspace, concatenate, std,\
	roll, diff, any, digitize, savetxt, append, array, savez_compressed, mean,\
	sqrt
from pickle import dump
import matplotlib.pyplot as plt
from ShapeGen import *
from CrossSection import *
from flow import *
from sediment import *
import sys
from os.path import exists
from os import makedirs, devnull

# Default parameters
params = {
	'rh':3.333e-5,
	'Q':1.,
	'Qs':30.0,
	'rho_s':2650,
	'D_s':0.001,
	'n':0.5,
	'sFile':'runParams.csv',
	'suffix':'',
	'suppress_save':0,
	'save_data':0,
	'rh_type':1,
	'suppress_print':0,
	't':10000
}

# Get parameters from the command line
for i in np.arange(1, len(sys.argv)):
	spl = sys.argv[i].split("=")
	var = spl[0]

	val=0
	if var!="sFile" and var!="suffix": val = float(spl[1])
	else: val = spl[1]
	params[var] = val

# Set parameters
suffix = params['suffix']

brh = params['rh']
print "Bedrock roughness height:", brh

Q = params['Q']
print "Discharge:", Q

Q_s = params['Qs']
print "Sediment discharge:", Q_s

n = params['n']
print "Erosional power:", n
k = 10

rho_s = params['rho_s']
print "Sediment density:", rho_s

D_s = params['D_s']
print "Sediment diameter:", D_s

S_rh = 6.8*D_s/30
print "Sediment Roughness height:", S_rh

rh_type = int(params['rh_type'])
if rh_type == 1: rh = brh; print "Roughness:", rh
elif rh_type == 2: rh = S_rh; print "Roughness:", rh
elif rh_type == 3: print "Using composite roughness"
elif rh_type == 4: print "Using per point roughness"
else: print "Invalid roughness type, defaulting to bedrock."

T_sc = calcTscr(rho_s, D_s)
print "Critical shields stress:", T_sc

# Done printing run info, do we want to print per ts info?
if params['suppress_print']: sys.stdout = open(devnull, "w")

sFile = params['sFile']
if not exists(sFile):
	f = open(sFile, "w")
	f.write("Q_w,Q_s,final_ts,W,area,D_s,n,brh,AvTau\n")
	#f.write("Q_w,Q_s,final_ts,W,area,D_s,n,brh,aspect\n")
	f.close()

# Settings
# total number of timesteps possible
t = int(params['t'])

# History arrays
###
###  Arrays to hold info from every time step.
###  Create with array = zeros(t) and add in
###  data with array[ts-1] = data.
###

widths = zeros(t+1)
Aw = zeros(t+1)
AvgTau_b = zeros(t+1)
aspect = zeros(t+1)
#xsecs = array([], dtype=object)
#alluvhs = array([])

# Generate initial shape and cross-section object

r = 0.5
x, y = genCirc(r)

## Rotate geometry so index 0 is the topmost point
y_roll = y.size - y.argmax()
x = roll(x, y_roll)
y = roll(y, y_roll)

cs = CrossSection(x, y)

# Loop variables
ts = 1
ah = 0 # alluviation height
j = 0

#Data save directory
#savedir = str(Q) + "-" + str(Q_s) + "-" + str(n) + "/"
#if not exists(savedir):
#	makedirs(savedir)
#Save x, y, ah for all time steps dict
#xyah = { 'x':[], 'y':[], 'ah':ah }

while(True):


	## Save cross-sections for every time step for movie making.
	## Comment out for mass runs.

#	xyah['x'] = cs.x
#	xyah['y'] = cs.y
#	xyah['ah'] = ah

#	f = open(savedir + str(j).zfill(6) + ".pkl", 'w')
#	dump(xyah, f)
#	f.close()
#	j+=1

	# Determine wetted cross-section
	## No alluviation layer..
	if ah == 0:
		wcs = CrossSection(cs.x, cs.y)
		# Index for bottom of cross-section
		y_min = cs.y.argmin()
		if rh_type==3 or rh_type==4: c_rh = rh = brh

	## Alluviated -- create wcs with an alluviation layer
	else:
		al_L, al_R = cs.findLR(ah)
		al_ypos = cs.y.min() + ah

		# Alluviation layer points
		al_y = al_ypos*ones(100)
		alx_diff = cs.x[al_R]/40
		al_x = linspace(cs.x[al_L]+alx_diff, \
			cs.x[al_R]-alx_diff, 100)

		# Wetted cross-section with discritized alluviaiton layer
		wx = concatenate( (cs.x[:al_L], al_x, cs.x[al_R+1:]) )
		wy = concatenate( (cs.y[:al_L], al_y, cs.y[al_R+1:]) )
		wcs = CrossSection(wx, wy)

		# Index for center of alluviation layer
		y_min = cs.y[:al_L].size + 49

		# Weighted roughness if type 3
		## In this case, we calculate roughness height by weighting
		## the bedrock roughness by the length of bedrock, and the sediment
		## roughness by the length of the alluviation layer. These weighted
		## roughnesses are divided by the total length of the wetted perimeter.
		## This is "Composite Roughness"
		if rh_type == 3 or rh_type == 4:
			left_length = calcL(cs.x[:al_L],cs.y[:al_L])
			sed_length = calcL(al_x, al_y)
			right_length = calcL(cs.x[al_R+1:], cs.y[al_R+1:])
			tot_length = calcL(wx, wy)
			c_rh = rh = ( brh*(left_length + right_length) + S_rh*sed_length ) \
				/tot_length

		if rh_type == 4:
			rh = ones(wcs.x.size)
			rh[:al_L] = brh
			rh[al_L:al_L+101] = S_rh
			rh[al_L+101:] = brh

	# Set maximum velocity location -- centroid
	wcs.calcShapeParams()
	uMx, uMy = wcs.findCentroid()
	wcs.setUMPoint(uMx, uMy)

	# Calculate flow parameters
	wcs.genRL()
	F = flow(wcs, Q, rh)
	F.calcChezy(c_rh if rh_type == 4 else False)
	F.calcFlowParams()

	## Get adjusted tau_b
	phi = arctan2( (wcs.umy - wcs.y), (wcs.umx - wcs.x) )
	alpha = arctan2(wcs.yp-wcs.ym, wcs.xp-wcs.xm)

	F.regenT_b(phi, alpha)

	# Locate bed-load layer, if not enough capacity set to 0, -1 to garuntee
	# alluviation
	q_t = calcLSTC(F.T_b, rho_s, D_s, T_sc)
	xLs, xRs = findBedLoad(wcs, y_min, q_t, Q_s)

	# Alluviate by D_s if BLL thickness > 5*D_s
	if (wcs.y[xLs] - wcs.y[y_min]) > 5*D_s:
		ah = ah+1*D_s
		continue

	# Redraw the cross-section as a function of shear stress
	rl = zeros(cs.x.size)

	## If alluviation hasn't occured yet redraw entire cross-section
	if ah==0: rl = 1e-4*F.T_b**(n)

	## If alluviation has occured extract redraw only for bedrock
	else:
		rl[:al_L+1] = k*1e-6*(F.T_b[:al_L+1])**(n)
		rl[al_R+1:] = k*1e-6*(F.T_b[al_L+100:])**(n)

	widths[ts-1] = wcs.x.max() - wcs.x.min()
	Aw[ts-1] = wcs.A
	AvgTau_b[ts-1] = F.T_b.mean()
	aspect[ts-1] = widths[ts-1] / ( wcs.y.max() - wcs.y.min() )

	if ((ts-1)%1000 == 0):
		print "Width:", wcs.x.max() - wcs.x.min()
		#print "RReynolds:", sqrt(9.81*wcs.A/wcs.P*F.S)*D_s/1e-6
		#if ts!=1: print std(aspect[ts-102:ts-2])
		#plt.plot(cs.x, cs.y, 'bo')
		#plt.plot(wcs.x, wcs.y, 'r.')
		#plt.plot(cs.x[0], cs.y[0], 'go')
		#plt.axis("equal")
		#plt.show()
		#plt.plot(phi, ".")
		#plt.plot(alpha,".")
		#plt.plot(phi-alpha)
		#plt.plot(sin(phi-alpha))
		#plt.plot(F.T_b)
		#plt.show()

		# Break if equilibrium width is reached
		#if ts != 1 and (widths[ts-1] - widths[ts-1001])<0:
		#	break

	cs.redraw(rl, resample=True)

	# Time step increment and break if time limit is reached
	ts+=1
	if (ts==t+1):
		print F.T_b.min(); break

plt.plot(cs.x, cs.y, 'bo')
plt.plot(wcs.x, wcs.y, 'r.')
plt.plot(cs.x[0], cs.y[0], 'go')
plt.axis("equal")
plt.show()

ws = widths[ts-102:ts-2]

eqWidth = round(mean(ws), 3)
print "Equilibrium width:", eqWidth

if params['save_data']==1:
	sd = zeros((ts, 4))
	sd[:,0] = widths[:ts]
	sd[:,1] = Aw[:ts]
	sd[:,2] = AvgTau_b[:ts]
	sd[:,3] = aspect[:ts]
	savetxt(str(Q)+"-"+str(Q_s)+"-"+str(n)+"-"+suffix+".csv", sd, delimiter=",", header="W,A,Tau_b,aspect")

	plt.plot(widths[:ts])
	plt.show()

	plt.plot(Aw[:ts])
	plt.show()

	plt.plot(AvgTau_b[:ts])
	plt.show()

	plt.plot(aspect[:ts])
	plt.show()


#last = {'x':cs.x, 'y':cs.y, 'ah':ah}
#dump(last, open( str(Q)+"-"+str(Q_s)+"-"+str(n)+"-last.pkl", 'w') )
#savez_compressed(str(Q)+"-"+str(Q_s)+"-"+str(n)+".npz", xsecs, alluvhs)

AreaF = Aw[ts-2]
AspectM = mean(aspect[ts-102:ts-2])
print AspectM

# Write out results to database unless supress save is enabled
if params['suppress_save']==1: sys.exit()
f = open(sFile, "a")

## Write with last Avg Tau b
f.write(str(Q) + ","+str(Q_s) + ","+str(ts-1)+","+str(eqWidth) + \
	","+str(AreaF) + ","+str(D_s) + ","+str(n) + "," + \
	str(brh) + "," + str(AvgTau_b[ts-2]) + "\n")

## With aspect ratio
#f.write("%f,%f,%d,%f,%f,%f,%f,%f,%f\n" %\
#	(Q, Q_s, ts-1, eqWidth, AreaF, D_s, n, brh, AspectM))

f.close()
