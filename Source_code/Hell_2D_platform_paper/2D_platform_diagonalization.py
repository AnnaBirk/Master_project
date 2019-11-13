from __future__ import division
import matplotlib.pyplot as plt
import kwant
import numpy as np
import tinyarray 					# array type optimized for small array sizes.
import scipy.sparse.linalg
import datetime

plt.rc('text', usetex=True)

print(datetime.date.today())

"""Pauli matrices:"""

s_0 = np.identity(2)
s_z = np.array([[1, 0], [0, -1]])
s_x = np.array([[0, 1], [1, 0]])
s_y = np.array([[0, -1j], [1j, 0]])

"""Paui matrices and tensor products in spin, e-h space:"""

# spin space first, then e-h space
sigma_y = tinyarray.array(np.kron(s_y, s_0))
sigma_z = tinyarray.array(np.kron(s_z, s_0))
# np.kron product is ~ a tensor product.
tau_z = tinyarray.array(np.kron(s_0, s_z))
tau_x = tinyarray.array(np.kron(s_0, s_x))
tau_y = tinyarray.array(np.kron(s_0, s_y))
sigma_xtau_z = tinyarray.array(np.kron(s_x, s_z))
sigma_ytau_z = tinyarray.array(np.kron(s_y, s_z))

"""Lattice type:"""

lat = kwant.lattice.square()

###

"""Classes and functions:"""

class SimpleNamespace(object):
	"""Contains parameters of the problem"""

	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)


def write_file_and_plot(f, par, Ez_values, energies, parameters, parameters_values):
	""" Writing log to file outfile, containing the system parameters, 
	the EZ_values and the eigenenergies, energies. Plots energies vs. Ez.

		Notes
		-----
		 - naming convention : 		
		 		EvsEz_ 	+ "Asym_" if asymmetryc/"Sym_" if symmetric device
						+ str(datetime.date.today())
						+ the following parameters, in order : 
								par.Nx,par.Ny,par.LM,par.LL,par.LR,par.LY,
								par.ax,par.ay,par.mu,par.Gamma,par.phiL,par.phiR
	"""
	# If asymmetric device :
	if par.asymmetric == True:
		# Physical parameters in title : 
		plt.title("$N_x$:%g,$N_x^A$:%g,$N_y$:%g,$\mu$:%g,$\Gamma$:%g,$\phi_L$:%g,$\phi_R$:%g,$a_x$:%.0f,$a_y$:%.0f" % (
			par.Nx, par.Nxasymmetric, par.Ny, par.mu, par.Gamma, par.phiL, par.phiR, par.ax, par.ay), fontsize=16)
		plt.xlabel("$E_z\, [meV]$", fontsize=16)
		plt.ylabel("$E_n\, [meV]$", fontsize=16)
		# Saving figure : 
		f.savefig("EvsEz_Asym_%s_%g_%g_%g_%g_%g_%g_%g_%g_%g_'alphabar'_%g_%g_%g.pdf" % (str(datetime.date.today(
		)), par.Nx, par.Ny, par.LM, par.LL, par.LR, par.LY, par.ax, par.ay, par.mu, par.Gamma, par.phiL, par.phiR))
		outfile = open(str("%s_datalog_Asym_EvsEz_%g_%g_%g_%g_%g_%g_%g_'alphabar'_%g_%g_%g.txt" % (str(
			datetime.date.today()), par.Nx, par.Ny, par.LM, par.LL, par.LR, par.LY, par.mu, par.Gamma, par.phiL, par.phiR)), "w")
	# If symmetric device :
	else:
		# Physical parameters in title : 
		plt.title("$N_x$:%g,$N_y$:%g,$\mu$:%g,$\Gamma$:%g,$\phi_L$:%g,$\phi_R$:%g,$a_x$:%.0f,$a_y$:%.0f" % (
			par.Nx, par.Ny, par.mu, par.Gamma, par.phiL, par.phiR, par.ax, par.ay), fontsize=16)
		plt.xlabel("$E_z\, [meV]$", fontsize=16)
		plt.ylabel("$E_n\, [meV]$", fontsize=16)
		# Saving figure :
		f.savefig("EvsEz_Sym_%s_%g_%g_%g_%g_%g_%g_%g_%g_%g_'alphabar'_%g_%g_%g.pdf" % (str(datetime.date.today(
		)), par.Nx, par.Ny, par.LM, par.LL, par.LR, par.LY, par.ax, par.ay, par.mu, par.Gamma, par.phiL, par.phiR))
		outfile = open(str("%s_datalog_Sym_EvsEz_%g_%g_%g_%g_%g_%g_%g_'alphabar'_%g_%g_%g.txt" % (str(
			datetime.date.today()), par.Nx, par.Ny, par.LM, par.LL, par.LR, par.LY, par.mu, par.Gamma, par.phiL, par.phiR)), "w")

	Ez_values_str = "np.linspace(%g,%g,%g)" % (
		Ez_values[0], Ez_values[-1], len(Ez_values))
	# Writing parameters, Zeeman field values and eigenenergies to file :
	outfile.write("Parameters=%s;\nEz_values=...;\nenergies=...: \n%s\n%s\n%s.\n\n" % (
		parameters, str(parameters_values), str(Ez_values_str), str(energies)))
	outfile.close()


def hopping_lw(site0, site1):
	"""
	Hopping linewidth between site0 and site1.
	Implements smaller linewidth in the y-direction.

	Uses the fact that to site naming convention, 
	sites above/below each other and next to each 
	other must have the y coupling. Else, they have the x coupling."""
	return 0.04 if np.abs(site0-site1) == 1 else 0.1


def H_nxy_as000110_print(sys, par, n_xystart, printtype):
	""" Function printing particular blocks of the Hamiltonian depending in the input "printtype".

		Parameters
		----------
		 - n_xystart = [to_sites,from_sites] for H00==H_n_xystart.
		 - printtype = 	off-diagonal (print term at n_xystart, and its neighbors) 
		 				LMboundary (print term of last site of Left region, and first site of Middle region) 
						MRboundary (print term of last site of Middle region, and first site of Right region)
		"""
	print("in Hnxy...")
	if printtype == 'off-diagonal':
		"""From this, also the neighboring off-diagonal blocks are printed."""
		# Diagonal block :
		H00 = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                                n_xystart[0]], from_sites=[n_xystart[1]])
		H00 = H00.tocsc()  
		H00 = H00.todense()
		print("H%.0f,%.0f for Ez=%g:" % (n_xystart[0], n_xystart[1], par.Ez))
		print(H00)

		# Off-diagonal block :
		H01 = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                                n_xystart[0]], from_sites=[n_xystart[1]+1])
		H01 = H01.tocsc()  
		H01 = H01.todense()
		print("H%.0f,%.0f for Ez=%g:" % (n_xystart[0], n_xystart[1]+1, par.Ez))
		print(H01)

		# Other off-diagonal block :
		H10 = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                                n_xystart[0]+1], from_sites=[n_xystart[1]])
		H10 = H10.tocsc()  
		H10 = H10.todense()
		print("H%.0f,%.0f Ez=%g:" % (n_xystart[0]+1, n_xystart[1], par.Ez))
		print(H10)

	elif printtype == "LMboundary":
		"""	Accesses hamiltonian for last site of L and first site of M 
			in order to check that Gamma(x) and phi(x) change as they 
			should on the boundary :"""
		LMboundary_left_Ny = (par.left[-1] + 1)*(par.Ny) - 1
		LMboundary_leftp1_0 = (par.left[-1] + 1)*(par.Ny)
		HL = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                               LMboundary_left_Ny], from_sites=[LMboundary_left_Ny])
		HM = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                               LMboundary_leftp1_0], from_sites=[LMboundary_leftp1_0])
		HLyym1 = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                                   LMboundary_left_Ny], from_sites=[LMboundary_left_Ny-1])
		HL = HL.todense()
		HM = HM.todense()
		HLyym1 = HLyym1.todense()
		print("H%.0f,%.0f Ez=%g:" % (LMboundary_left_Ny, LMboundary_left_Ny, par.Ez))
		print(HL)
		print("H%.0f,%.0f Ez=%g:" %
		      (LMboundary_leftp1_0, LMboundary_leftp1_0, par.Ez))
		print(HM)
		print("H%.0f,%.0f Ez=%g:" %
		      (LMboundary_left_Ny, LMboundary_left_Ny-1, par.Ez))
		print(HLyym1)

	elif printtype == "MRboundary":
		"""Same as for LM, but now last site of M and first site of R."""
		MRboundary_right_Ny = (par.middle[-1] + 1)*(par.Ny) - 1
		MRboundary_rightp1_0 = (par.middle[-1] + 1)*(par.Ny)
		HM = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                               MRboundary_right_Ny], from_sites=[MRboundary_right_Ny])
		HR = sys.hamiltonian_submatrix(args=[par], sparse=True, to_sites=[
		                               MRboundary_rightp1_0], from_sites=[MRboundary_rightp1_0])
		HM = HM.todense()
		HR = HR.todense()
		print("H%.0f,%.0f Ez=%g:" %
		      (MRboundary_right_Ny, MRboundary_right_Ny, par.Ez))
		print(HM)
		print("H%.0f,%.0f Ez=%g:" %
		      (MRboundary_rightp1_0, MRboundary_rightp1_0, par.Ez))
		print(HR)

def Gamma_phi_fn(site, p):
	"""
	Returns Gamma and phi for a given site, for the parameters specified by p.

	Parameters
	----------
	 - site : 		kwant.builder.site object, specifying current site 
	 				for which computations are made.
	
	 - p : 			SimpleNamespace object, containing system parameters.
	"""
	# If in the Left region, Gamma is GammaL, and phi is phiL :
	if site.pos[0] <= par.left[-1]:  # site.pos[0] is the x value of the position of the site.
		Gamma = p.GammaL
		phi = p.phiL
	# If in the Middle region, there is no superconductivity : 
	elif par.middle[0] <= site.pos[0] <= par.middle[-1]:
		Gamma = 0
		phi = 0
	# If in the Right region, Gamma is GammaR, and phi is phiR :
	elif par.right[0] <= site.pos[0] <= par.right[-1]:
		Gamma = p.GammaR
		phi = p.phiR
	else:
		raise ValueError(
				"In Gamma_phi_fn: site.pos[0] was in neither part of the system. \
				Cannot assign Gamma- and phi-values.")
	return [Gamma, phi]


def onsite(site, p):
	""" Onsite term in Hamiltonian : """
	[Gamma, phi] = Gamma_phi_fn(site, p)
	return (2*(p.tx+p.ty) - const.hbar**2/(2*p.m_star)*p.mu)*tau_z + p.Ez*sigma_y/2 + Gamma*(np.cos(phi)*tau_x - np.sin(phi)*tau_y)


def hoppingx(site0, site1, p):
	""" Nearest-neighbor hopping term in Hamiltonian, x-direction : """
	return -p.tx*tau_z + 1j*p.alphahbar/(2*p.ax)*sigma_ytau_z


def hoppingy(site0, site1, p):
	""" Nearest-neighbor hopping term in Hamiltonian, y-direction : """
	return -p.ty*tau_z - 1j*p.alphahbar/(2*p.ay)*sigma_xtau_z


def make_system(p):
	"""Function building system containing all relevant regions, 
	for parameters from SimpleNamespace object p.
	Returns the kwant.builder object "sys. """
	sys = kwant.Builder()

	# if no right superconducting 2D layer, the device is asymmetric.
	if len(par.right) == 0: 	# if asymmetric
		par.asymmetric = True  	# to be used when printing last blocks of the Hamiltonian
		par.Nxasymmetric = len(par.left) + len(par.middle)
	else:
		par.asymmetric = False

	"""Building the different parts by adding onsite/hoppings to sys:"""
	# Looping over the three regions : 
	for i, part_ in zip(('L', 'M', 'R'), (par.left, par.middle, par.right)):

		"""Printing which part is being built at the current loop iteration:"""
		if len(part_) == 0:  # if asymmetric device, the length of the <part_> is zero 
							 # for either the Left or the Right region.
			if part_ != par.right:
				print("%s: %s (site %s to %s)" %
				      (str(datetime.datetime.now()), i, part_[0], part_[-1]))
			else:
				print("%s: R abscent" % str(datetime.datetime.now()))

		else:
			print("%s: %s (site %s to %s)" %
			      (str(datetime.datetime.now()), i, part_[0], part_[-1]))

		""" Tight-binding hamiltonian for all sites inside <part_>:"""
		sys[(lat(x, y) for x in part_ for y in range(p.Ny))] = onsite
		sys[kwant.builder.HoppingKind((1, 0), lat)] = hoppingx
		sys[kwant.builder.HoppingKind((0, 1), lat)] = hoppingy

	return sys

###


def main(**kwargs):
	"""Making system:"""
	sys = make_system(par)

	"""Finalizing system:"""
	sys = sys.finalized()

	"""Plotting system graph:"""
	plt.rcParams.update({'font.size': 14})
	fig = kwant.plot(sys, site_symbol='o', hop_color='pink',
	                 site_edgecolor='red', site_lw=0.01, hop_lw=hopping_lw)

	"""Calculating energies for a range of applied fields, Ez:"""
	Ez_values = np.linspace(0, 0.8, 80)  # [Ez] = meV
	energies = []

	"""Printing parts of the Hamiltonian the first time main() is called:"""
	if counter_Hprint == 0:  # if test here won't slow down program significantly because loop we're inside is only over different Nx, which there aren't usually that many of.
		par.Ez = Ez_values[0]

		"""	printing different parts of H, specified by printtype.
			printtype=off-diagonal: Printing H00, H01 and H10 for given par.Ez."""
		H_nxy_as000110_print(sys, par, n_xystart=[0, 0], printtype="off-diagonal")

	# Looping over Zeeman field values :
	for par.Ez in Ez_values:

		H = sys.hamiltonian_submatrix(args=[par], sparse=True)
		H = H.tocsc()  

		# Find k eigenvalues and eigenvectors of the real symmetric square matrix 
		# or complex hermitian matrix A.[...] 
		# k: number of eigenvalues/vecs desired.[...]. 
		# sigma: Find eigenvalues near sigma using shift-invert mode :
		eigs = scipy.sparse.linalg.eigsh(H, k=50, sigma=0) 	# these are the energies of the k first eigenstates. 
															# Element eigs[1] are the eigenvectors. After the loop, 
															# energies are collected in a nested array, which can be 
															# accessed as follows: energies[0:len(Ez_values)][0:k]
		energies.append(np.sort(eigs[0]))

	"""Plotting energies vs Ez:"""
	f = plt.figure()
	plt.plot(Ez_values, energies)
	plt.show()

	"""Writing energies to log file and save plot:"""
	write_file_and_plot(f, par, Ez_values, energies,
                     parameters="Nx,LM,LL,LR,LY,mu,'alphahbar',Gamma,phiL,phiR",
                     parameters_values=(par.Nx, par.LM, par.LL, par.LR, par.LY, par.mu, 'alphahbar', par.Gamma, par.phiL, par.phiR))

###

""" Parameters """
const = SimpleNamespace(c=2.99792458e17, m_e=0.5109989461e9,
                        hbar=6.582119514e-13) 	# [c]=nm/s, [m_e]=meV/c^2, [hbar]=meV*s
par = SimpleNamespace(LM=250., LL=1.e3, LR=1.e3, LY=4.e3, mu=0,
                      	alphahbar=1.42e-4*const.c*const.hbar, Gamma=180e-3, phiL=0, phiR=0)
						# [L]=nm, [mu]=meV, [alpha]=nm*meV, [Gamma]=meV, [phi]=1
# GammaL/R, the effective superconducting gap, used in gamma_phi_fn(site,p) :
par.GammaL = par.Gamma
par.GammaR = par.Gamma

counter_Hprint = 0 		# If set to zero, parts of the Hamiltonian is printed
						# the first time main() is run. Can be used as a check.

# Looping over asymmetric (Nx=160) ands symmetric (Nx=260) devices :
for par.Nx in [120, 260]:
	# rounding off to nearest int. 1.625 defines the aspect ratio Nx/Ny.
	par.Ny = int(round(float(par.Nx/1.625)))
	par.ax = float(par.LM+par.LL+par.LR)/par.Nx
	par.ay = float(par.LY)/par.Ny
	par.m_star = 0.023*const.m_e/(const.c)**2
	par.tx = const.hbar**2/(2*0.023*const.m_e/(const.c)**2*par.ax**2)
					# [t]=meV(=(meV*s)^2/(meV/(nm/s)^2*nm^2))
	par.ty = const.hbar**2/(2*par.m_star*par.ay**2)
	# Left and Middle regions (defined by the following domains along the x-axis) :
	par.left = np.arange(0, round(float(par.LL)/par.ax))
	par.middle = np.arange(round(float(par.LL)/par.ax),
	                       round(float(par.LL+par.LM)/par.ax))
	# Right region :
	for par.right in (np.arange(	round(float(par.LL+par.LM)/par.ax), 				# asymmetric device
									round(float(par.LL+par.LM+par.LR)/par.ax)), []): 	# symmetric device
		print("Nx=%s, Ny=%.0f; ax=%g, ay=%g" % (par.Nx, par.Ny, par.ax, par.ay))

		main()

	counter_Hprint += 1

print("Program finished %s" % str(datetime.datetime.now()))