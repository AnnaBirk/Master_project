import matplotlib.pylab as plt

plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 10})


class SimpleNamespace(object):
	"""Contains parameters of the problem"""
	def __init__(self,**kwargs):
		self.__dict__.update(kwargs)

def make_system(a=1,W=1,L=10,Deltapos=4,Delta=0.1,t=1.,barrierpos=(3,4),barrier=1.5,mu=0.4):
	"""
	N-S 1D junction with barrier at junction.

	Input
	-----
		norbs:	number of orbitals per site. is 2, e and h orbital, for SC.

		W : 	width of NS device in the y-direction.		
		
		sites 0 to deltapos : 	Normal metal region.

		barrierpos : 			position of potential barrier. Coincides with deltapos.

		sites Deltapos to L: 	superconducting region (sites) for which the order parameter couples sites.
		
	Constructing Kwant Bulders
	--------------------------
		Scattering region : 	N - barrier - S

		Left lead : 			Normal metal

		Right lead : 			Superconducting

	Observation
		Left lead : 			Normal metal

		Right lead : 			Superconducting

	Observation
	-----------
		Varying W:
		 - W = 1: reproduces BTK.
		 - W > 1: only resembles BTK for barrier >~ 1.

	"""

	lat = kwant.lattice.square(norbs=2) #define lattice type. 'norbs' is necessary to define here if we want to look at the scattering matrix elements for electrons and holes separately.
	syst = kwant.Builder() # Initiate building system. After, add pairings.

	"""Defining scattering region:"""
	
	syst[(lat(x,y) for x in range(Deltapos) for y in range(W))] = (2*t-mu)*tau_z 	# N region
	syst[(lat(x,y) for x in range(Deltapos,L) for y in range(W))] = (2*t-mu)*tau_z + Delta*tau_x 	# SC region

	"""Making tunnel barrier:
	This overrides the previously set points at the positions given in barrierpos"""

	syst[(lat(x,y) for x in range(barrierpos[0], barrierpos[1]) for y in range(W))] = (2*t - mu + barrier)*tau_z

	"""Hoppings between neighbors:"""
	
	syst[lat.neighbors()] = -t*tau_z

	"""Defining leads:
		LEFT = NORMAL
		RIGHT = SUPERCONDUCTING
	"""
	
	sym_left = kwant.TranslationalSymmetry((-a,0))

	lead0 = kwant.Builder(sym_left, conservation_law=-tau_z, particle_hole=tau_y)
	lead0[(lat(0, j) for j in range(W))] = (2 * t - mu) * tau_z
	lead0[lat.neighbors()] = -t * tau_z

	# Right lead - superconducting, so Delta is nonzero
	sym_right = kwant.TranslationalSymmetry((a, 0))
	lead1 = kwant.Builder(sym_right)
	lead1[(lat(0, j) for j in range(W))] = (2 * t - mu) * tau_z + Delta * tau_x
	lead1[lat.neighbors()] = -t * tau_z
	
	""" Attach the leads and return the system:"""
	syst.attach_lead(lead0)
	syst.attach_lead(lead1)

	return syst


def calc_conductance(syst,energies):
	G = []		# list for appending conductance G = N - R_ee + R_he at the left lead
	for energy in energies:
		smatrix = kwant.smatrix(syst,energy)
		N = smatrix.submatrix((0, 0), (0, 0)).shape[0]
		G.append(N -\
						smatrix.transmission((0, 0), (0, 0)) +\
						smatrix.transmission((0, 1), (0, 0)))

	np.save("Conductance_B%.1f.npy"%barrieri,G)
	return G

def plot_conductance(energies,G,figfilename):
	fig = plt.figure(barrieri)
	plt.plot(energies,G)

	plt.xlabel("E")
	plt.ylabel("dI/dV")
	plt.title("barrier = %.1f, $\Delta = %.1f$" %(barrieri,0.1))
	plt.savefig(figfilename,bbox_inches='tight')
	plt.hold("on")

def plot_conductances(energies,conductances_and_barriers,transmission_N,figfilename,sizexy=(1.8,1.8),nolegend=True):
	"""
	Input
	-----
	conductances_and_barriers : 	Array. Shape: (4, 2)
									first column contains conductances, second column contains the barrier values used for the particular conductance.

	"""
	fig, ax  = plt.subplots(1,1,figsize=sizexy)
	legend = []
	i = 0
	color = 'C0' 
	alpha = 1
	
	for conductance, barrier, transmission_N_ in zip(conductances_and_barriers[:,0],conductances_and_barriers[:,1],transmission_N):
		plt.plot(energies,conductance,color=color,alpha=alpha)
		legend.append("$B=%.1f\,meV$"%barrier)
		plt.axhline(y=transmission_N_, linestyle=':',alpha=alpha)
		alpha -= 0.2
		
	plt.xlabel("$eV$")
	plt.ylabel("$dI/dV\ [e^2/h]$")
	if nolegend==False:
		plt.legend(legend)
	plt.xlim(0,0.2)
	plt.ylim(0,2)
	plt.subplots_adjust(left=0.18,bottom=0.24,right=0.94,top=0.96)
	plt.savefig(figfilename,bbox_inches="tight")
	plt.show()

def check_PHS(syst):
	# Scattering matrix
	s = kwant.smatrix(syst, energy=0)
	# Electron to electron block
	s_ee = s.submatrix((0,0), (0,0))
	# Hole to hole block
	s_hh = s.submatrix((0,1), (0,1))
	print('s_ee: \n', np.round(s_ee, 3))
	print('s_hh: \n', np.round(s_hh[::-1, ::-1], 3))
	print('s_ee - s_hh^*: \n',
	np.round(s_ee - s_hh[::-1, ::-1].conj(), 3), '\n')
	# Electron to hole block
	s_he = s.submatrix((0,1), (0,0))
	# Hole to electron block
	s_eh = s.submatrix((0,0), (0,1))
	print('s_he: \n', np.round(s_he, 3))
	print('s_eh: \n', np.round(s_eh[::-1, ::-1], 3))
	print('s_he + s_eh^*: \n',
	np.round(s_he + s_eh[::-1, ::-1].conj(), 3))


import numpy as np
import kwant
import matplotlib.pyplot as plt
tau_z = np.array([[1,0], [0,-1]])
tau_y = np.array([[0,-1j], [1j,0]])
tau_x = np.array([[0,1],[1,0]])
hbar = 4.135667e-12		# units: meV*s

barriers = np.array([0,0.5,1.5,5.0])

par = SimpleNamespace(t=0.5,L=10,Deltapos=4,barrierpos=(3,4))

par.transmission_N = 1./(1+barriers**2) 	# standard transmission probability as a function of the barrier strength, for scattering off a delta-potential

for barrieri in barriers:
	syst = make_system(barrier = barrieri,t=par.t,L=par.L,Deltapos=par.Deltapos,barrierpos=par.barrierpos)

	kwant.plot(syst)

	syst = syst.finalized()

	# check_PHS(syst)

	Energies = np.linspace(0,0.2,201)
	conductance = calc_conductance(syst,Energies)
	plot_conductance(Energies,conductance,figfilename="Conductance_Kwant_barrier%.1f_Delta%.1f_W1.pdf" %(barrieri, 0.1))

Energies = np.linspace(0,0.2,201)

G_B0 = np.load("Conductance_B0.0.npy")
G_B05 = np.load("Conductance_B0.5.npy")
G_B1 = np.load("Conductance_B1.5.npy")
G_B5 = np.load("Conductance_B5.0.npy")

conductances_and_barriers = np.array([[G_B0,0.0],[G_B05,0.5], [G_B1,1.0], [G_B5,5.0]])

plot_conductances(Energies,conductances_and_barriers,par.transmission_N,figfilename="Conductances_BAll_t%.1f_L%.0f_barrierpos%.0f_%.0f.pdf"%(par.t,par.L,par.barrierpos[0],par.barrierpos[1]))