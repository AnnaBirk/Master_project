"""
Program diagonalizing a tight-binding version of an 
Oreg-Lutchyn Hamiltonian, for a one-dimensional wire.
"""

import matplotlib.pyplot as plt
import kwant
import numpy as np
import tinyarray
import scipy.sparse.linalg

class SimpleNamespace(object):
	"""A simple container for parameters."""
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)

""" 2-by-2 Pauli matrices : """
s_0 = np.identity(2)
s_z = np.array([[1, 0], [0, -1]])
s_x = np.array([[0, 1], [1, 0]])
s_y = np.array([[0, -1j], [1j, 0]])

""" 4-by-4 Pauli matrices in electron-hole times spin space : """
tau_z = tinyarray.array(np.kron(s_z, s_0))
tau_x = tinyarray.array(np.kron(s_x, s_0))
sigma_z = tinyarray.array(np.kron(s_0, s_z))
tau_zsigma_x = tinyarray.array(np.kron(s_z, s_x))

def onsite(site, p):
	""" Onsite term in Hamiltonian : """
	return tau_z * (p.mu - 2 * p.t) + \
						sigma_z * p.B + tau_x * p.Delta

def hopping(site0, site1, p):
	""" Nearest-neighbor hopping term in Hamiltonian : """
	return tau_z * p.t + 1j * tau_zsigma_x * p.alpha

def make_system(l=70):
	""" 
		Making and finalizing a kwant.builder object describing the system 
		graph of a closed, one-dimensional wire with l number of sites. 
	"""
	sys = kwant.Builder()
	lat = kwant.lattice.chain()
	sys[(lat(x) for x in range(l))] = onsite
	sys[lat.neighbors()] = hopping
	return sys.finalized()

# Make and finalize system graph :
sys = make_system()

# Calculate and plot lowest eigenenergies in B-field :
B_values = np.linspace(0, 0.6, 251)
energies = []

# SimpleNamespace object holding physical parameters :
params = SimpleNamespace(	t=1,		# tight-binding parameter 
							mu=-0.1, 	# chemical potential (a.u.)
							alpha=0.05, # spin-orbit coupling strength
							Delta=0.2)	# superconducting gap (setting energy scale)

# Looping over different values of B-field :
for params.B in B_values:
	H = sys.hamiltonian_submatrix(args=[params], sparse=True)
	H = H.tocsc()
	# Diagonalizing the Hamiltonian :
	eigs = scipy.sparse.linalg.eigsh(H, k=20, sigma=0)
	energies.append(np.sort(eigs[0]))

def cm2inch(num_cm):
	"""
		Conversion function from centimeters to inches, 
		used in plotting function.
	"""
    inch = 2.54
    return num_cm/inch

def plot_energies(B_values, energies, Delta=params.Delta):
	""" 
		Plot energies in units of Delta.
	"""
	plt.rc('text', usetex=True)
	font_size = 10
	plt.rcParams.update({'font.size': font_size})
	fig, ax = plt.subplots(nrows=1, ncols=1)
	fig.set_size_inches(cm2inch(4),cm2inch(4))
	fig.subplots_adjust(left=0.17,bottom=0.28,top=1,right=1)

	# Plotting energies :
	color="gray"

	plt.plot(B_values/Delta, [energies[i]/Delta for i in range(len(energies))],color=color)

	# Energy axis settings :
	legend_00 = ax.legend(	[r"$E_n\ [\Delta]$"],
							frameon=True,
							loc="upper center",
							borderpad=0.1,
							handletextpad=0.1,
							edgecolor="white")
	ax.set_xlabel(r"$\rm{B-field}\ [\Delta]$")

	plt.show()

plot_energies(B_values, energies)