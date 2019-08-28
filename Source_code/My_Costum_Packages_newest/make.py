



import kwant
from My_Costum_Packages_newest import misc
import datetime
import numpy as np
import tinyarray		#Tinyarrays are similar to NumPy arrays, but optimized for small sizes. Common operations on very small arrays are to 3-7 times faster than with NumPy (with NumPy 1.6 it used to be up to 35 times), and 3 times less memory is used to store them. Tinyarrays are useful if you need many small arrays of numbers, and cannot combine them into a few large ones.



"""Pauli matrices:"""

s_0 = np.identity(2)
s_z = np.array([[1,0],[0,-1]])
s_x = np.array([[0,1],[1,0]])
s_y = np.array([[0,-1j],[1j,0]])

"""Paui matrices and tensor products in spin, p-h space:"""

sigma_y = tinyarray.array(np.kron(s_y,s_0)) # spin space first, then e-h space
sigma_z = tinyarray.array(np.kron(s_z,s_0))
tau_z = tinyarray.array(np.kron(s_0,s_z)) # kron product is ~ a tensor product that works on matrices. So for my use, they're the same.
tau_x = tinyarray.array(np.kron(s_0,s_x))
tau_y = tinyarray.array(np.kron(s_0,s_y))
sigma_xtau_z = tinyarray.array(np.kron(s_x,s_z))
sigma_ytau_z = tinyarray.array(np.kron(s_y,s_z))

sigma_ytau_y = tinyarray.array(np.kron(s_y,s_y))

I_x = np.array([[1,0],[0,-1]])	# rotation around x-axis


"""Lattice type:"""

lat = kwant.lattice.square(norbs=4)




def Gamma_phi_fn(site,p):
	"""
	Calculates the value of Gamma and phi given the input site and the input parameters.
	
	Parameters
	----------
	p : 	Namespace class containing the parameters of the system

	Notes
	-----
	site.pos[0] is the x value of the position. !!! Could also implement as heavyside step function, but this works as it is.
	"""

	if site.pos[0] <= p.left[-1]:	# 
		Gamma = p.GammaL; phi = p.phiL
	elif p.middle[0] <= site.pos[0] <= p.middle[-1]:
		Gamma = 0; phi = 0
	elif p.right[0] <= site.pos[0] <= p.right[-1]:
		Gamma = p.GammaR; phi = p.phiR
	else:
		raise ValueError("In Gamma_phi_fn: site.pos[0] was in neither parts of the system. Cannot assign Gamma- and phi-values.")
	return [Gamma,phi]


def onsite_1D_semiconductor(site,p):
	"""
	Onsite energy of one-dimensional semiconductor, whose Hamiltonian is given as
			H_k = (hbar^2*k^2/2m - mu)c_k,sigma^dagger c_k,sigma. 

	Parameters
	----------
	p.t_N : 	tunneling amplitude inside the semiconductor.

	p.mu : 		chemical potential of semiconductor.

	Notes
	-----
	 - have not defined own parameter mu_N for the semimetal, in order for it to work better with different parts of the code that already makes use of p.mu, e.g. when running calculations over a range of different p.mu. So bear in mind that now that onsite_1D_semiconductor is used, p.mu means the chemical potential of the semiconductor, and not that of e.g. a system that is a proximitized semiconductor. In other words, p.mu is chosen to be the general chemical potential of what is considered the system, independently of what type of system is considered.

	??? : + p.Ez*sigma_y/2 ok?
	"""
	return (2*p.t_N - p.mu)*tau_z + p.Ez*sigma_y/2.


def hoppingy_1D_semiconductor(site,p):
	"""
	Hopping energy of one-dimensional semiconductor, whose Hamiltonian is given as in onsite_1D_semiconductor.

	Parameters
	----------
	p.t_N : 	tunneling amplitude inside the semiconductor.

	"""
	return -p.t_N*tau_z + 1j*p.alphahbar/(2*p.ay)*sigma_ytau_z 


def hoppingy_2D_superconductor_pincher(site0,site1,p):
	"""
	Hopping energy in y-direction of two-dimensional superconductor that is e.g. attached to a SC lead. In Kwant, it becomes part of the system Builder object, not that of the SC lead, while in the real model, it is thought of as part of the lead.

	Parameters
	----------
	p.ty_SC_pincher : 	hopping in y-direction of pincher sites.
	ay_eff : 			effective lattice constant corresponding to chosen bandwidth p.ty_SC_pincher

	Notes
	-----
	 - Pincher can be effected by altering p.tx_SC_pincher and p.ty_SC_pincher such that they are not equal to p.tx and p.ty.

	"""

	ay_eff = p.ay# !!!???[1] np.abs(np.sqrt(p.hbar**2/(2*p.m_star*p.ty_SC_pincher)))
	return - p.ty_SC_pincher*tau_z - 1j*p.alphahbar/(2.*ay_eff)*sigma_xtau_z


def hoppingx_2D_superconductor_pincher(site0,site1,p):
	"""
	Hopping energy in y-direction of two-dimensional superconductor that is e.g. attached to a SC lead. In Kwant, it becomes part of the system Builder object, not that of the SC lead, while in the real model, it is thought of as part of the lead.

	Parameters
	----------
	p.tx_SC_pincher : 	hopping in x-direction of pincher sites.
	ax_eff : 			effective lattice constant corresponding to chosen bandwidth p.tx_SC_pincher
	
	Notes
	-----
	 - Pincher can be effected by altering p.tx_SC_pincher and p.ty_SC_pincher such that they are not equal to p.tx and p.ty.

	"""

	ax_eff = p.ax # !!!???[1] np.abs(np.sqrt(p.hbar**2/(2*p.m_star*p.tx_SC_pincher)))
	return - p.tx_SC_pincher*tau_z + 1j*p.alphahbar/(2.*ax_eff)*sigma_ytau_z


def onsite_2D_superconductor_pincher(site,p):
	"""
	Onsite energy of 2D superconductor, e.g. used for same purpose as hoppingy_2D_superconductor_pincher.

	Parameters
	----------
	p.Delta_pincher : 		Superconducting gap at pincher site.
							If e.g. whole pincher is effected by the ts, then one may set p.Delta_pincher equal to Gamma.

	p.mu_SC_pincher : 		Chemical potential of pincher site.

	p.Delta_pincher : 		Superconducting gap of pincher sites. May be made the same as the superconducting lead to be attached, while altering the ts to encompas the uneven interface.

	??? : including p.Ez*sigma_y/2 ok? Assumes also that p.Ez uniform/same as inside the semiconductor.

	p.Ez_2D_S : 			Effective Ez in sc

	"""

	return (2*(p.tx_SC_pincher+p.ty_SC_pincher) - p.mu_SC_pincher)*tau_z + p.Ez_2D_S*sigma_y/2. + p.Delta_pincher*tau_x

def hoppingy_1D_superconductor_pincher(site0,site1,p):
	"""
	Hopping y dir sc 1S
	"""

	ay_eff = p.ay# !!!???[1] np.abs(np.sqrt(p.hbar**2/(2*p.m_star*p.ty_SC_pincher)))
	return - p.ty_SC_pincher*tau_z - 1j*p.alphahbar/(2.*ay_eff)*sigma_xtau_z


def onsite_1D_superconductor_pincher(site,p):
	"""
	Onsite energy of 1D superconductor, e.g. used for same purpose as hoppingy_2D_superconductor_pincher.
	"""
	Ez_eff = p.Ez
	return (2*p.ty_SC_pincher - p.mu_SC_pincher)*tau_z + Ez_eff*sigma_y/2. + p.Delta_pincher*tau_x




def hoppingy_2D_superconductor_lead(site0,site1,p):
	"""
	Hopping energy in y-direction of two-dimensional SC lead.

	Parameters
	----------
	p.ty_SC_lead: 	hopping in y-direction of pincher sites.
	ay_eff : 		the lattice constant in the y-direction that corresponds to the bandwidth p.ty_SC_lead that is chosen in set_params.py. Taking ay_eff to be a positive length.

	Notes
	-----
	 - Pincher can be effected by altering p.tx_SC_pincher and p.ty_SC_pincher such that they are not equal to p.tx and p.ty.

	"""

	ay_eff = np.abs(np.sqrt(p.hbar**2/(2*p.m_star*p.ty_SC_lead)))
	return - p.ty_SC_lead*tau_z + 1j*p.alphahbar/(2.*ay_eff)*sigma_ytau_z 


def hoppingx_2D_superconductor_lead(site0,site1,p):
	"""
	Hopping energy in y-direction of two-dimensional SC lead.

	Parameters
	----------
	p.tx_SC_lead : 	hopping in x-direction of pincher sites.
	ax_eff : 		the lattice constant in the x-direction that corresponds to the bandwidth p.tx_SC_lead that is chosen in set_params.py. Taking ax_eff to be a positive length.
	
	Notes
	-----
	 - Pincher can be effected by altering p.tx_SC_pincher and p.ty_SC_pincher such that they are not equal to p.tx and p.ty.

	"""

	ax_eff = np.abs(np.sqrt(p.hbar**2/(2*p.m_star*p.tx_SC_lead)))
	return - p.tx_SC_lead*tau_z + 1j*p.alphahbar/(2.*ax_eff)*sigma_ytau_z 


def onsite_2D_superconductor_lead(site,p):
	"""
	Onsite energy of 2D superconducting lead, e.g. used for same purpose as hoppingy_2D_superconductor_lead.

	Parameters
	----------
	p.Delta_lead : 		Superconducting gap at pincher site.
							May be made the same as the superconducting lead to be attached, while altering the ts to encompas the uneven interface. If e.g. whole pincher is effected by the ts, then one may set p.Delta_pincher equal to Gamma.

	p.mu_SC_lead : 		Chemical potential of pincher site.

	p.Ez_2D_S_lead : 	Ez in SC lead, determined from its g-factor

	??? : including p.Ez*sigma_y/2 ok? Assumes also that p.Ez uniform/same as inside the semiconductor.

	"""


	return (2*(p.tx_SC_lead+p.ty_SC_lead) - p.mu_SC_lead)*tau_z + p.Ez_2D_S_lead*sigma_y/2. + p.Delta_SC_lead*tau_x


def onsite(site,p):
	"""
	Proximitized superconductor onsite energy for a 2D rectangular lattice.

	Parameters
	----------
	p.tx : 	tx = hbar^2/((m*/c^2)*ax^2); [tx]=meV=(meV*s)^2/(meV/(nm/s)^2*nm^2).
	p.ty : 	(x<->y)
	p.mu :	chemical potential of the system whose onsite Hamiltonian is being calculated; [mu] = meV.
	p.Ez :	Energy corresponding to magnetic field along the z-direction; [Ez] = meV.
	Gamma :	proximitized superconducting tunneling. Functional dependence on spacial coordinate(s) is given by the function Gamma_phi_fn(site,p); [Gamma] = meV
	phi : 	Superconducting phase.
	"""

	[Gamma,phi] = Gamma_phi_fn(site,p)
	return (2*(p.tx+p.ty) - p.mu)*tau_z + p.Ez*sigma_y/2. + Gamma*(np.cos(phi)*tau_x - np.sin(phi)*tau_y)


def onsite_1D(site,p):
	"""
	Onsite energy for a 1D lattice.

	Parameters - see onsite function
	----------
	"""

	[Gamma,phi] = Gamma_phi_fn(site,p)
	return (2*p.ty - p.mu)*tau_z + p.Ez*sigma_y/2. + Gamma*(np.cos(phi)*tau_x - np.sin(phi)*tau_y)


def onsite_1D_pincher(site,p):
	"""
	Parameters - see onsite function
	----------
	p.pincher : 	True if pincher is attached.
	p.t_pincher : 	Value of ty for pincher at current site. Note that if you want a different value of t_pincher, you need to call onsite_1D() separately after re-defining p.t_pincher = (e.g.) p.t_pincherL or p.t_pincherR.
	"""

	[Gamma,phi] = Gamma_phi_fn(site,p)
	if p.pincher == True:
		return (2*p.t_pincher - p.mu)*tau_z + p.Ez*sigma_y/2. + Gamma*(np.cos(phi)*tau_x - np.sin(phi)*tau_y)
	else:
		ValueError("p.pincher != True in onsite_1D_pincher. Need it to be True in order to compute onsite_1D Hamiltonian for the t-value p.t_pincher.")


def hoppingx(site0,site1,p):
	"""
	Hopping energy in x-direction of 1D system with proximitized superconductivity.

	Parameters - see onsite function
	----------
	p.alphabar : 	velocity parameter of Rashba spin-orbit interaction.
	
	Notes
	-----
	This is the same functional form as for a 2D/1D system.
	"""

	return -p.tx*tau_z + 1j*p.alphahbar/(2.*p.ax)*sigma_ytau_z 


def hoppingy(site0,site1,p):
	"""
	Hopping energy i y-direction of 1D system with proximitized superconductivity.

	Parameters and Notes - see onsite and hoppingx functions
	--------------------
	"""

	return -p.ty*tau_z - 1j*p.alphahbar/(2.*p.ay)*sigma_xtau_z

def hoppingy_t_(site0,site1,p):
	"""
	Hopping energy in 1D system for a, arbitrary hopping p.t_.

	Parameters
	----------
	p.t_ : 	Arbitrary hopping amplitude for which hopping is computed.

	Notes
	-----
	Used e.g. in implementing a pincher. In that case, e.g. t_=p.t_pincher_L or R.
	Since part of the hopping Hamiltonian for which Kwant demands have only the arguments site0, site1, p, t_ needs to be saved in p in order to use this functoin, and cannot be passed as an extra argument to hoppingy_t_.
	"""
	return -p.t_*tau_z - 1j*p.alphahbar/(2.*p.ay)*sigma_xtau_z


def check_asymm_device(p):
	"""
	Check if device is asymmetric (one Al strip) or symmetric (two Al strips).

	Returns
	-------
	p.asymmetric : 	Boolean
					Is used later in naming files and plot titles, and initiates saving the size p.Nxasymmetric of the SC strip and middle retion altogether to the Simplenamespace object p containing the parameters of the system.
	"""

	if len(p.right) == 0:
		p.asymmetric = True	# to be used when printing last blocks of the Hamiltonian
		p.Nxasymmetric = len(p.left) + len(p.middle)
	else:
		p.asymmetric = False
	return p.asymmetric


def print_current_part(part_,i):
	"""
	Prints which part of the system that one is currently attaching to the Builder object sys=kwant.Builder(...).

	Parameters
	----------
	i : 	str.
			Name of the part one is currently building in.
	part_ : 	np.array, list or tuple.
				Contains all the site numbers of the part currently being built.
	"""
	
	if len(part_) == 0:	# if asymmetric device such that a part is simply removed, effectively rendering its length=0
		if part_ != par.right:
			print("%s: %s (site %s to %s)" %(str(misc.round_time(datetime.datetime.now(),round_to=60)),i,part_[0],part_[-1]))
		else:
			print("%s: R abscent" %str(misc.round_time(datetime.datetime.now(),round_to=60)))
	else:
		print("%s: %s (site %s to %s)" %(str(misc.round_time(datetime.datetime.now(),round_to=60)),i,part_[0],part_[-1]))


def make_1D_system(p,ppar):

	"""Function building superconductor for parameters in SimpleNamespace object p.

	Parameters
	----------
	p.pincher : 	If True, make 1D system with an extra site attached with a different tunneling constant, tL and tR on the left and right, respectively. This models the pincher which is affecting this one additional site's hopping amplitude tL/R at the edge of the system.
	p.left : 	Left strip is implemented as the SC 1D system. So p.left = [0], while p.middle and p.right should be empty ([]) when making the 1D system using this function.

	Functions called
	----------------
	onsite_1D : 	onsite different now for 1D system. Not input phase/gamma explicitly such that may specify variables such as Ez at later point in code after finalizing. I think that Kwant has this particular syntax in sys[lat..]=<function_name> without having to input arguments in the function.

	Notes
	-----
	??? (see below)
	"""

	import kwant

	sys = kwant.Builder()

	sys[(lat(x,y) for x in p.left for y in range(p.Ny))] = onsite_1D#(lat(0,0),p) 	#!!!???
	sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy#(lat(0,0),lat(0,1),p)

	if p.pincher == True:
		print(" - Adding pincher")
		"""Attaching extra sites with hoppings p.tL/R:"""
		sys[lat(p.left[0],p.Ny-1)] = onsite_1D_pincher#(lat(p.left[0],p.Ny-1),p)
		sys[lat(p.left[0],0)] = onsite_1D_pincher#(lat(p.left[0],p.Ny-1),p)
		# sys[(p.left,p.)]	#???: How to implement hopping for one specific site, not using neighbors or HoppingKind?
		# sys[0,0] = hoppingy_t_

	"""Attaching leads to this "L" region:"""
	lead = make_lead_onesite_wide_1D(p)

	sys.attach_lead(lead)
	sys.attach_lead(lead.reversed())

	""" 
	Adding SC lead last. In that way, the indexing of the smatrix will work as before when wanting to access the other leads. In order to access from the SC lead, acces the 2nd lead index in e.g. the smatrix.
	"""
	if ppar.SClead == True:
		SClead = make_SC_lead_Ny_wide(p)
		sys.attach_lead(SClead)
		print(" - Attached SC lead")

	return sys

def make_1D_N_2D_SC_system(p,ppar):

	"""

	Notes
	-----
	 - x-coordinate of the 1D semiconducting system is taken to be zero.
	 - x-coordinate of the 2D pincher layer is taken to be at +1. So the pincher layer is 1D, but implemented as a 2D Hamiltonian because it needs to have couplings between it and the middle N.

	"""

	import kwant

	sys = kwant.Builder()

	sys[(lat(x,y) for x in [0] for y in range(p.Ny))] = onsite_1D_semiconductor 
	sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy_1D_semiconductor

	if ppar.pincher_SC_lead == True:
		print(" - Adding SC pincher at SM system - SC lead boundary.")
		"""Attaching extra sites with hoppings p.tx_SC_pincher and p.ty_SC_pincher:"""
		if ppar.SC_pincher_1D == True:
			sys[(lat(x,y) for x in [1] for y in range(p.Ny))] = onsite_1D_superconductor_pincher
			sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy_1D_superconductor_pincher

		else:
			sys[(lat(x,y) for x in [1] for y in range(p.Ny))] = onsite_2D_superconductor_pincher
			sys[kwant.builder.HoppingKind((1,0),lat)] = hoppingx_2D_superconductor_pincher
			sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy_2D_superconductor_pincher

	"""Attaching leads to this "L" region:"""
	lead = make_lead_onesite_wide_1D(p)

	sys.attach_lead(lead)
	sys.attach_lead(lead.reversed())

	""" 
	Adding SC lead last. In that way, the indexing of the smatrix will work as before when wanting to access the other leads. In order to access from the SC lead, acces the 2nd lead index in e.g. the smatrix.
	"""
	if ppar.SClead == True:
		SClead = make_SC_lead_Ny_wide(p)
		sys.attach_lead(SClead)
		print(" - Attached SC lead")

	return sys


def make_1D_NLeft_1D_S_Heff_No_NRight(p,ppar):
	""" 
	Make N-S(Heff) function for calibrating against

	"""
	print("%s: in make_1D_NLeft_1D_S_Heff_No_NRight()" %str(misc.round_time(datetime.datetime.now(),round_to=60)))

	sys = kwant.Builder()

	sys[(lat(x,y) for x in p.left for y in range(p.Ny))] = onsite_1D
	sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy

	if p.pincher == True:
		print(" - Adding pincher")
		"""Attaching extra sites with hoppings p.tL/R:"""
		sys[lat(p.left[0],p.Ny-1)] = onsite_1D_pincher
		sys[lat(p.left[0],0)] = onsite_1D_pincher

	"""Attaching leads to this "L" region:"""
	lead = make_lead_onesite_wide_1D(p)
	sys.attach_lead(lead)

	return sys


def make_1D_Nleft_1D_N_2D_S_2D_SMiddle_No_NRight(p,ppar):
	"""
	Make N-N(S) subject to the calibration of the N-S(Heff) system created in make_1D_NLeft_1D_S_Heff_No_NRight.

	Notes
	-----
	 - x-coordinate of the 1D semiconducting system is taken to be zero.
	 - x-coordinate of the 2D pincher layer is taken to be at +1. So the pincher layer is 1D, but implemented as a 2D Hamiltonian because it needs to have couplings between it and the middle N.

	"""
	print("%s: in make_1D_Nleft_1D_N_2D_S_2D_SMiddle_No_NRight()" %str(misc.round_time(datetime.datetime.now(),round_to=60)))

	import kwant

	sys = kwant.Builder()

	sys[(lat(x,y) for x in [0] for y in range(p.Ny))] = onsite_1D_semiconductor 
	sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy_1D_semiconductor

	if ppar.pincher_SC_lead == True:
		print(" - Adding SC layer (2D Hamiltonian and 1D line of sites) at N system - S lead boundary.")
		"""Attaching extra sites with hoppings p.tx_SC_pincher and p.ty_SC_pincher:"""
		sys[(lat(x,y) for x in [1] for y in range(p.Ny))] = onsite_2D_superconductor_pincher
		sys[kwant.builder.HoppingKind((1,0),lat)] = hoppingx_2D_superconductor_pincher
		sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy_2D_superconductor_pincher

	"""Attaching Left lead to the middle region:"""
	lead = make_lead_onesite_wide_1D(p)
	sys.attach_lead(lead)

	""" 
	Adding SC lead last. In that way, the indexing of the smatrix will work as before when wanting to access the other leads. In order to access from the SC lead, acces the 2nd lead index in e.g. the smatrix.
	"""
	if ppar.SClead == True:
		SClead = make_SC_lead_Ny_wide(p)
		sys.attach_lead(SClead)
		print(" - Attached SC lead")

	return sys


def make_system(p,ppar):

	"""
	Building the different parts by adding onsite/hoppings to sys:
	Since adding lead only to scattering region, I iterate over par.middle first; then I can use just one sys object to attach the leads in the first iteration, and keep the same sys object when continuing on adding the right and left parts.

	Parameters
	----------
	ppar.make_N_wire : 	Boolean.
							If True or 1, a 1D semiconducting wire is constructed with a 2D superconductor connected to it along the whole wire (along the y-direction), with a pincher layer of SC (???given by a 2D hamiltonian but only with one site in width for now ???ok???) at the boundary between the semiconductor and the superconducting lead.

	Notes
	-----
	If 1D pincher is to be added at the edges of the 1D system, this is taken care of in make_1D_system, given that p.pincher==True and a t-value p.t_pincher needs to have been defined.
	"""
	
	print("%s: in make_system()" %str(misc.round_time(datetime.datetime.now(),round_to=60)))

	if ppar.one_D == True:
		
		if ppar.make_N_wire:
			sys = make_1D_N_2D_SC_system(p,ppar)				# make 1D N-N(S)-N system, where middle N is 1D and middle S is 2D.
		elif ppar.make_1D_NLeft_1D_S_Heff_No_NRight:
			sys = make_1D_NLeft_1D_S_Heff_No_NRight(p,ppar) 	# make 1D N-S junction for calibration, where S has the effective Hamiltonian of the proximitized semiconductor
		elif ppar.make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight:	# make 1D N-N(S) function to calibrate, where S is a 2D SC Hamiltonian that should induce the gap in the middle N
			sys = make_1D_Nleft_1D_N_2D_S_2D_SMiddle_No_NRight(p,ppar)
		else:
			sys = make_1D_system(p,ppar)						# make 1D N-S-N system, where middle S is the effective Hamiltonian of the proximitized semiconductor, and there is no middle S lead, meaning that conductances above the gap obrained in this model are not accurate, as modes should in the real physical system have been able to escape through the middle S lead above the gap (where the middle S has available quasiparticle excitations).

	else:
		"""Function building superconductor for parameters in SimpleNamespace object p."""
		sys = kwant.Builder()
		par.asymmetric = check_asymm_device(par)	# True/False, to be used when printing last blocks of Hamiltonian

		for i,part_ in zip(('M','L','R'),(p.middle, p.left, p.right)):
			"""Printing which part is being built at the current loop iteration:"""
			print_current_part(part_,i)
			"""Adding parts to the system:"""
			sys[(lat(x,y) for x in part_ for y in range(p.Ny))] = onsite #  Not input phase/gamma explicitly such that may specify variables such as Ez at later point in code after finalizing. I think that Kwant has this particular syntax in sys[lat..]=<function_name> without having to input arguments in the function.
			sys[kwant.builder.HoppingKind((1,0),lat)] = hoppingx
			sys[kwant.builder.HoppingKind((0,1),lat)] = hoppingy

			"""Attaching leads to the middle region, i.e. the scattering region:"""
			if i == 'M':
				lead = make_lead_onesite_wide(p) #make_lead(p)
				sys.attach_lead(lead)
				sys.attach_lead(lead.reversed())

	return sys


def make_SC_lead_Ny_wide(p):
	"""
	Middle (M) lead - superconducting.
	
	Parameters
	----------
	Parameters are that of superconducting lead, namely


	"""

	sym_SC = kwant.TranslationalSymmetry((1, 0))
	leadSC = kwant.Builder(sym_SC)	# don't need to specify particlehole symmetry explicitly, as it is implicitly implemented in the SC Hamiltonian.

	leadSC[(lat(0, j) for j in range(p.Ny))] = onsite_2D_superconductor_lead
	leadSC[kwant.builder.HoppingKind((0,1),lat)] = hoppingy_2D_superconductor_lead
	leadSC[kwant.builder.HoppingKind((1,0),lat)] = hoppingx_2D_superconductor_lead
	
	return leadSC


def make_lead(p):
	
	"""	Make metallic lead that can be attached to sys (Builder) object with the same conservation_law and particle_hole symmetry as the 2D proximitized topological superconductor being considered.
		Has the width of par.middle, and is translationally invariant in the y-direction."""

	sys_ = kwant.Builder(kwant.TranslationalSymmetry([0,1]),conservation_law=tinyarray.array(np.kron(s_z,I_x)),particle_hole=sigma_ytau_y)	## ???: symmetries - implementing complex conjugation?
	sys_[(lat(x,0) for x in par.middle)] = (2*(p.tx+p.ty) - p.mu)*tau_z
	sys_[kwant.builder.HoppingKind((1,0),lat)] = -p.tx*tau_z
	sys_[kwant.builder.HoppingKind((0,1),lat)] = -p.ty*tau_z

	return sys_


def make_lead_onesite_wide_1D(p):
	sys_ = kwant.Builder(kwant.TranslationalSymmetry([0,1]),conservation_law=sigma_z,particle_hole=tau_y) #conservation_law=sigma_z???
	tx = p.t_N_lead
	muL= tx # par.tx + 2*par.tx*(1-np.cos(np.pi*par.ax/par.LM))	# mu in lead
	middle_xcoordinate = np.int(np.average(p.left))
	sys_[lat(middle_xcoordinate,0)] = np.zeros((4,4))# sys_[(lat(x,0) for x in middleL)] = 0 #(2*(par.tx+par.ty) - muL)*tau_z
	sys_[kwant.builder.HoppingKind((0,1),lat)] = -tx*tau_z
	
	return sys_

