from pylab import *
from subprocess import call, Popen
from time import gmtime, asctime, time
import threading
import sys

# Default save/file options
opts = {
	'sFile':'runOutputs.csv',
	'trf':'torun.csv',
	'sedpart':0
}

# Get options from cmd
for i in np.arange(1, len(sys.argv)):
	spl = sys.argv[i].split("=")
	var = spl[0]
	val = spl[1]
	opts[var] = val

# Set save file and runs file
sFile = opts['sFile']
print sFile
trf = opts['trf']
spart = opts['sedpart']

# Run with timeout
class CallPara(threading.Thread):

	def __init__(self, cmd, to):
		threading.Thread.__init__(self)
		self.cmd = cmd
		self.timeout = to

	def run(self):
		self.p = Popen(self.cmd)
		self.p.wait()

	def Run(self):
		self.start()
		self.join(self.timeout)

		if self.is_alive():
			self.p.terminate()
			self.join()

# Save array to file
def lwrite():
	savetxt(trf, data, delimiter=',', \
		header='Qw,Qs,rh,rho_s,D_s,n,t,r_type,Ran')

# Create command to be run
def callcmd(Q, Qs, rh, rho_s, D_s, n, t, rh_type):
	return ["python", "paragenesis.py", \
		Q, Qs, rh, rho_s, D_s, n, t, rh_type,\
		"sFile=" + sFile, "sc=" + spart, "suppress_print=0"]

# Pull run parameters out of csv
data = genfromtxt(trf, delimiter=",", skip_header=1)

Qw = data[:,0]
Qs = data[:,1]
rh = data[:,2]
rho_s = data[:,3]
D_s = data[:,4]
n = data[:,5]
t = data[:,6]
r_type = data[:,7]
Ran = data[:,8]

# Loop through CSV. If ran column is not TRUE run paragenesis with those params
for i in arange(Qw.size):
	if Ran[i]:
		continue

	startt = time()
	print "Run:", Qw[i], Qs[i]
	Q_s = "Q=" + str(Qw[i])
	Qs_s = "Qs=" + str(Qs[i])
	rh_s = "rh=" + str(rh[i])
	rho_s_s = "rho_s=" + str(rho_s[i])
	D_s_s = "D_s=" + str(D_s[i])
	n_s = "n=" + str(n[i])
	t_s = "t=" + str(int(t[i]/4.))
	r_s = "rh_type=" + str(r_type[i])

	print "Starting time:", asctime(gmtime()) + "..."

	# Time out if things ran for more than 10 minutes.. possible hangup
	CallPara(callcmd(Q_s, Qs_s, rh_s, rho_s_s, D_s_s, n_s, t_s, r_s), \
		120*60).Run()
	endt = time()
	Ran[i] = True
	lwrite()
	print "Done at:", asctime(gmtime()) + \
		". Run took:", endt-startt, "seconds.."
