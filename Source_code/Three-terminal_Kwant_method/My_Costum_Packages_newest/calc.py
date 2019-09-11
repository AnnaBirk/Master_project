"""

PACKAGE FOR CALCULATIONS IN KWANT AND FROM RESULTANT QUANTITIES FROM KWANT.

"""
import kwant
from My_Costum_Packages_newest import misc, make, plot
import datetime
import numpy as np


def calc_energies_mu(sys,par,ppar):

	import scipy

	mu_values_1_ = par.mu_values
	if ppar.doubleres == False:
		pass
	else:
		mu_values_2_ = par.mu_values_2
		par.mu_values = mu_values_2_ 	# after done with calculation, re-defining par.mu_values to its original definition, locally saved in mu_values_1_

	if ppar.pincher_SC_lead == True:
		par.energies_mu = []
		par.eigenvectors_mu = []#np.zeros((len(par.mu_values),par.k))
		print(" - sites of N wire only is: ", sites_of_N_wire_only)
		print(" - sites of S wire only is: ", sites_of_S_wire_only)
		for par.mu, i in zip(par.mu_values,range(len(par.mu_values))):

			print(" - mu-value number:", i)
			
			sites_of_N_wire_only = [site for site in range(0,par.Ny)]
			sites_of_S_wire_only = [site for site in range(par.Ny,2*par.Ny)]

			H = sys.hamiltonian_submatrix(args=[par],sparse=True) # !!!???[2],from_sites=sites_of_N_wire_only,to_sites=sites_of_N_wire_only)
			# H = sys.hamiltonian_submatrix(args=[par],sparse=True,from_sites=sites_of_S_wire_only,to_sites=sites_of_S_wire_only)
			H = H.tocsc()
			eigs = scipy.sparse.linalg.eigsh(H,k=par.k,sigma=0)

			sorting_indices = eigs[0].argsort()
			evals = eigs[0][sorting_indices]
			evecs = eigs[1][:,sorting_indices]
			
			par.energies_mu.append(evals)
			par.eigenvectors_mu.append(evecs)
	else:
		par.energies_mu = []
		par.eigenvectors_mu = []#np.zeros((len(par.mu_values),par.k))
		for par.mu, i in zip(par.mu_values,range(len(par.mu_values))):

			print(" - mu-value number:", i)
			 # print(" - Hamiltonian submatrix to be calculated")
			H = sys.hamiltonian_submatrix(args=[par],sparse=True)
			 # print(" - Hamiltonian submatrix calculated")
			H = H.tocsc()
			eigs = scipy.sparse.linalg.eigsh(H,k=par.k,sigma=0)

			sorting_indices = eigs[0].argsort()
			evals = eigs[0][sorting_indices]
			evecs = eigs[1][:,sorting_indices]
			
			par.energies_mu.append(evals)
			par.eigenvectors_mu.append(evecs)

	par.mu_values = mu_values_1_


def calc_G_11_12_S_A(sys_,p,ppar,i,var_name,elapsed_time,biasenergies_asymm=False):

	"""	Computing conductance. Returns G_11.

		Parameters
		----------
		sys_ :  	finalized system with all parts and leads.
		p :  		same parameters as are given to the functions making up the hamiltonian (which is also the parameters given to the kwant.Hamiltonian method).
		energies : 	computed eigenenergies in main().
		biasenergies_symm : 	is True if an odd number of biasenergies where every negative value has a corresponding positive value at the same number of places from the value zero, and vice versa. E.g. biasenergies = -3,-2,-1,0,1,2,3. Is False by default. User needs to ensure that biasenergies are "symmetric" like this before entering biasenergies_symm=True as input.
	
		Notes
		-----
		 - if ppar.Sclead = True, SC lead is attached to the middle region. The lead index for the SC lead is 2, because it is atted after the normal leads. Thus, use index 2 if wanting to access the SC lead values in e.g. the smatrix.
	"""

	print("%.0f'th %s-value, Dt=%.1f" %(i,var_name,elapsed_time))
	G_11 = []
	G_12 = []

	if ppar.oneNLead == True:									# One lead means no nonlocal conductance is to be calculated
		for biasenergy in p.biasenergies:
			smatrix = kwant.smatrix(sys_,biasenergy,args=[p])	#Takes longer for larger system
			
			G_11.append(2 - \
							smatrix.transmission((0,0),(0,0)) + \
							smatrix.transmission((0,1),(0,0)))
	else:
		for biasenergy in p.biasenergies:
			smatrix = kwant.smatrix(sys_,biasenergy,args=[p])	#Takes longer for larger system
																# New in kwant version 1.4: args=[p]-->params=[p] inside kwant.smatrix
			
			G_11.append(2 - \
							smatrix.transmission((0,0),(0,0)) + \
							smatrix.transmission((0,1),(0,0)))

			G_12.append( - smatrix.transmission((0,0),(1,0)) + \
							smatrix.transmission((0,1),(1,0)))

	G_11_S = [(G_11[i] + G_11[-1-i])/2. for i in range(0,len(G_11))]
	G_11_A = [(G_11[i] - G_11[-1-i])/2. for i in range(0,len(G_11))]
	G_12_S = [(G_12[i] + G_12[-1-i])/2. for i in range(0,len(G_12))]
	G_12_A = [(G_12[i] - G_12[-1-i])/2. for i in range(0,len(G_12))]

	return G_11,G_12,G_11_S,G_11_A,G_12_S,G_12_A


def calc_u2mv2_mu(evec_n, mu_values, Ny, filename_, filename):
	"""
	Calculates u and v (Bogoliubon amplitudes) from the eigenvectors, which are in themselves functions of mu/Ez and are given for every site {0,...,Ny} in the 1D system of length Ny.
	 - Definition:
					u(i)^2 = |u_up(i)|^2 + |u_down(i)|^2
					v(i)^2 = |v_up(i)|^2 + |v_down(i)|^2
					where i = site number, up/down = spins and the eigenvectors in the sigma*tau (spin*particle/hole) basis are {|psi(i)>}, given by
					|psi(i)> = (u_up(i), u_down(i), v_down(i), -v_down(i)).
	
	Parameters
	----------
	evec_n : 	{|psi(site i, mu j)>|_n} with shape (par.mu_res,par.Ny*2*2)

	Variables
	---------
	uu, ud, vd, mvu : 	u_up, u_down, v_down, -v_up. Shape: (par.mu_res,par.Ny), not with par.Ny*2*2, since up/down and u/v already specifies the spin/p-h DOFs.
	"""

	print(misc.round_time(datetime.datetime.now(),round_to=60),": in calc_u2mv2_var")

	uu = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	ud = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	mvd = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	vu = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	
	u2_mu_site = np.zeros((len(mu_values),Ny))#,dtype=np.complex_)
	v2_mu_site = np.zeros((len(mu_values),Ny))#,dtype=np.complex_)

	for i in range(len(mu_values)):
		print(" - mu_value number ", i)
		for j in range(Ny):

				[uu[i,j],mvd[i,j],ud[i,j],vu[i,j]] = evec_n[i,4*j:4*(j+1)]

				u2_mu_site[i,j] = np.abs(uu[i,j])**2 + np.abs(ud[i,j])**2
				v2_mu_site[i,j] = np.abs(mvd[i,j])**2 + np.abs(vu[i,j])**2

	u2_mu = np.sum(u2_mu_site,axis=1)
	v2_mu = np.sum(v2_mu_site,axis=1)

	### Also saving the whole u/v2_mu_site file. Then, may access it later, and e.g. plot only for one given site:
	np.save("u2_mu_sites_"+filename_,u2_mu_site)
	np.save("v2_mu_sites_"+filename_,v2_mu_site)
	######################################

	u2mv2_mu = u2_mu-v2_mu

	np.save(filename,u2mv2_mu)

	print(" - Saved datalog file: ", filename)

	return u2mv2_mu, u2_mu, v2_mu, u2_mu_site, v2_mu_site



def calc_u2mv2_var(evec_n, var_values, Ny, filename, savefiles, **kwargs):
	"""
	Calculates u and v (Bogoliubon amplitudes) from the eigenvectors, which are in themselves functions of mu/Ez and are given for every site {0,...,Ny} in the 1D system of length Ny.
	 - Definition:
					u(i)^2 = |u_up(i)|^2 + |u_down(i)|^2
					v(i)^2 = |v_up(i)|^2 + |v_down(i)|^2
					where i = site number, up/down = spins and the eigenvectors in the sigma*tau (spin*particle/hole) basis are {|psi(i)>}, given by
					|psi(i)> = (u_up(i), u_down(i), v_down(i), -v_down(i)).
	
	Parameters
	----------
	evec_n : 		{|psi(site i, mu j)>|_n} with shape (par.mu_res,par.Ny*2*2)

	**kwargs : 		

		write_u_v_var_LR : 	Any value, e.g. True or an int.
							If given, then u and v is also saved, as a function of the variable var, but summed over sites. Also, u and v values are returned for left and right leads. These are defined as the u and v-values at the first (site 0) and last lattice site (site Ny) when e.g. 1D.!!!

	Variables
	---------
	uu, ud, vd, mvu : 	u_up, u_down, v_down, -v_up. Shape: (par.mu_res,par.Ny), not with par.Ny*2*2, since up/down and u/v already specifies the spin/p-h DOFs.
	
	
	Returns
	-------
	if no kwargs given : 	u2mv2_var, u2_var, v2_var
	if one kwarg given : 	u2mv2_var, u2_var, v2_var, u_var, v_var

	"""

	from My_Costum_Packages_newest import misc
	import datetime
	import numpy as np

	print(misc.round_time(datetime.datetime.now(),round_to=60),": in calc_u2mv2_var")

	uu = np.zeros((len(var_values),Ny),dtype=np.complex_)
	ud = np.zeros((len(var_values),Ny),dtype=np.complex_)
	mvd = np.zeros((len(var_values),Ny),dtype=np.complex_)
	vu = np.zeros((len(var_values),Ny),dtype=np.complex_)
	
	u2_var_site = np.zeros((len(var_values),Ny))#,dtype=np.complex_)
	v2_var_site = np.zeros((len(var_values),Ny))#,dtype=np.complex_)

	if len(kwargs) == 1:
		u_var_site = np.zeros((len(var_values),Ny),dtype=np.complex_)
		v_var_site = np.zeros((len(var_values),Ny),dtype=np.complex_)

		for i in range(len(var_values)):
			print(" - var_value number ", i)
			for j in range(Ny):

					[uu[i,j],mvd[i,j],ud[i,j],vu[i,j]] = evec_n[i,4*j:4*(j+1)]

					u2_var_site[i,j] = np.abs(uu[i,j])**2 + np.abs(ud[i,j])**2
					v2_var_site[i,j] = np.abs(mvd[i,j])**2 + np.abs(vu[i,j])**2

					u_var_site[i,j] = uu[i,j] + ud[i,j]		# ??? correct def right?
					v_var_site[i,j] = mvd[i,j] + vu[i,j] 	# ??? same

	elif len(kwargs) == 0:
		for i in range(len(var_values)):
			print(" - var_value number ", i)
			for j in range(Ny):

					[uu[i,j],mvd[i,j],ud[i,j],vu[i,j]] = evec_n[i,4*j:4*(j+1)]

					u2_var_site[i,j] = np.abs(uu[i,j])**2 + np.abs(ud[i,j])**2
					v2_var_site[i,j] = np.abs(mvd[i,j])**2 + np.abs(vu[i,j])**2

	else:
		raise ValueError("kwargs in calc.calc_u2mv2_var is too long or type not understood.")

	u2_var = np.sum(u2_var_site,axis=1)
	v2_var = np.sum(v2_var_site,axis=1)

	if len(kwargs) == 1:
		u_var = np.sum(u_var_site,axis=1)
		v_var = np.sum(v_var_site,axis=1)

	u2mv2_var = u2_var-v2_var

	if savefiles == True:
		np.save("u2mv2_"+filename,u2mv2_var)
		np.save("u2_" + filename,u2_var)
		np.save("v2_" + filename,v2_var)

		print(" - Saved datalog files: ", "u2mv2_"+filename, "\n", "u2_" + filename, "\n", "v2_" + filename)

	if savefiles == True & len(kwargs) == 1:
		np.save("u_"+filename,u_var)
		np.save("v_"+filename,v_var)

		print(" - Saved datalog files: ", "u2mv2_"+filename, "\n", "u2_" + filename, "\n", "v2_" + filename, "\n", "u_" + filename, "\n", "v_" + filename)

	if len(kwargs) == 0:
		return u2mv2_var, u2_var, v2_var

	elif len(kwargs) == 1:
		return u2mv2_var, u2_var, v2_var, u_var, v_var, u_var_site[:,0], u_var_site[:,Ny-1], v_var_site[:,0], v_var_site[:,Ny-1]


def calc_u2mv2_mu_basisredef(evec_n, mu_values, Ny, filename):
	"""

	Redef basis of calc_u2mv2_mu
	"""

	print(misc.round_time(datetime.datetime.now(),round_to=60),": in calc_u2mv2_mu_basisredef")

	uu = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	mvu = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	ud = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	vd = np.zeros((len(mu_values),Ny),dtype=np.complex_)
	
	u2_mu_site = np.zeros((len(mu_values),Ny))
	v2_mu_site = np.zeros((len(mu_values),Ny))

	for i in range(len(mu_values)):
		print(" - mu_value number ", i)
		for j in range(Ny):

				[uu[i,j],mvu[i,j],ud[i,j],vd[i,j]] = evec_n[i,4*j:4*(j+1)]

				u2_mu_site[i,j] = np.abs(uu[i,j])**2 + np.abs(ud[i,j])**2
				v2_mu_site[i,j] = np.abs(mvu[i,j])**2 + np.abs(vd[i,j])**2

	u2_mu = np.sum(u2_mu_site,axis=1)
	v2_mu = np.sum(v2_mu_site,axis=1)

	u2mv2_mu = u2_mu-v2_mu

	np.save(filename,u2mv2_mu)
	print(" - Saved datalog file: ", filename)

	return u2mv2_mu, u2_mu, v2_mu





	E1_mu, Em1_mu, evec1_mu, evecm1_mu = calc.calc_E1_Em1_from_E0_E0prime_mu(E0 = par.energies_mu[:,n[0]], E0prime = par.energies_mu[:,n[1]])


def calc_E1_Em1_from_E0_E0prime_mu(E0,E0prime,evec0,evec0prime,mu_values,N,tol=2):
	"""
	######## DISCLAIMER: ###############

	THIS IS AN OLD FUNCTION. YOU MAY USE AND UPDATE THE FUNCTION alc_E1_Em1_from_E0_E0prime_var INSTEAD.

	####################################

	Calculate electron/hole energies E1 and Em1 from cutting and pasting the energies E0 and E0prime whose values are in
			E0 >= 0, E0prime <= 0
	and whose period is half that of E1 and Em1. E1 and Em1 have max/min values
			max(E1 and Em1) = max(E0), min(E1 and Em1) = min(E0prime).

	This function calculates the derivtive of the energies. When the derivative changes, that's where the cut needs to be done!

	Parameters
	----------
	E0 : 			lowest absolute energy for the solved system, which corresponds to the positive values of the lowest energy pair of modes.

	E0prime : 		lowest absolute energy of the solved system, whose values are negative.

	evec0 : 		eigenvector with the energy E0.
					 - indexing : 	evec0[var_res,N]
					 				where N = 2*2*Ny in spinful system with p-h dofs and is 1D with length Ny. var_res is the number of points in the variable one is considering the energies/evecs a function of, which could be e.g. mu_res or Ez_res.

	evec0prime : 	eigenvector with the energy E0prime.
					 - indexing : 	same as with evec0.

	splineE0 : 					Fitted 'function' of E0. Used to find the derivative of E0.

	splineE0_derivative : 		derivative of fitted 'function' splineE0. Used to find the indices where the sign of the derivative changes, which corresponds to the indices where the energies E1 and Em1 cross as well as the indices where E0 and E0prime have maxima/minima.

	splineE0prime : 			equivalently to above.

	splineE0prime_derivative : 	equivalently to above.

	tol : 						tolerance in indices for which elements in sign_change_indices_E0max and sign_change_indices are considered to be the same. Necessary because they are calculated in different manners, and it thus turns out that when they are supposed to be the same, they vary slightly from each other. Example: for mu_res = 1500, tol=1 would do the job of excluding the points in sign_change_indices_E0max from sign_change_indices. Default value of tol is set to tol = 2.

	sign_change_indices : 		indices in variable (e.g. mu_values) where the derivative of the energy E0 (or equivalently E0prime) changes.

	sign_change_indices_actual_E0max_withintol : 	Indices of sign_change_indices corresponding to the indices in sign_change_indices_E0max. They are calculated by cheching if the indices in sign_change_indices_E0max are inside sign_change_indices_actual_E0max_withintol within the tolerance set by tol.

	m : 											indexes in variable (e.g. mu_values) for which the elements of sign_change_indices that are not equal to the elements of sign_change_indices_actual_E0max_withintol. These are the points where E0and E0prime cross, which is where the cutting and pasting procedure is to be performed.

	Notes
	-----
	 - The sign of the derivative of the energy changes also in between the desired points where it changes. So we need to identify which points are the desired points and which to throw away. We can identify that the points we want are the ones where the energy is approximately zero, while the points to throw away are the ones where the absolute value of the energy is at its maxima.
	 - All elements of the energies E1 and Em1 that are less than m[0] and larger than m[-1] need to be set equal to the appropriate E0 and E0prime. The convention used here is that
	 			 - for indices < m[0] : 	E1 = E0, Em1 = E0prime.
	 Depending on whether len(m) is even or odd, we get that E0 and Em1 should be for indices > m[-1]:
	 			 - indices < m[0] : 			E1_cut = E0prime, Em1_cut = E0.
				 - odd/even indices > m[-1] : 	opposite
	"""

	from scipy.interpolate import UnivariateSpline
	import numpy as np
	import scipy.signal
	from My_Costum_Packages_newest import misc
	import datetime

	print("%s: In calc.calc_E1_Em1_from_E0_E0prime_mu()"%misc.round_time(datetime.datetime.now(),round_to=60))

	splineE0 = UnivariateSpline(mu_values, E0, s=0)
	splineE0prime = UnivariateSpline(mu_values, E0prime, s=0)

	splineE0_derivative = splineE0.derivative(n=1)
	splineE0prime_derivative = splineE0prime.derivative(n=1)

	sign_change_indices = []
	sign_change_indices_E0max = scipy.signal.argrelmax(np.abs(E0))
	sign_change_indices_E0max = np.array(sign_change_indices_E0max[0])

	for i in range(len(splineE0_derivative(mu_values))-1):
		if np.sign(splineE0_derivative(mu_values)[i]) != np.sign(splineE0_derivative(mu_values)[i+1]):
			sign_change_indices.append(i+1)

	"""Getting the actual indices in terms of the indices of sign_change_indices that correspond to the indices in sign_change_indices_E0max:"""
	sign_change_indices_actual_E0max_withintol = []

	for maximum in sign_change_indices_E0max:
		for delta in range(tol):
			# print(" - delta, ", delta)
			if maximum+delta in sign_change_indices:
				# print("hello", maximum)
				sign_change_indices_actual_E0max_withintol.append(maximum+delta)
			elif maximum-delta in sign_change_indices:
				# print("hello", maximum)
				sign_change_indices_actual_E0max_withintol.append(maximum-delta)

	mask = np.isin(sign_change_indices,sign_change_indices_actual_E0max_withintol)

	m = []
	for mask_i,element in zip(mask,sign_change_indices):
		# if True, meaning different elements
		if mask_i == False:
			m.append(element)
	print("m is: ", m)

	if m == []:
		""" If no points where sign of energy changes, then keep the energies as they are """
		E1_cut = E0prime
		Em1_cut = E0
		evec1_cut = evec0prime
		evecm1_cut = evec0

		return E1_cut, Em1_cut, evec1_cut, evecm1_cut

	else:
		E1_cut = np.zeros(len(mu_values))
		Em1_cut = np.zeros(len(mu_values))

		evec1_cut = np.zeros((len(mu_values),N),dtype=np.complex_)
		evecm1_cut = np.zeros((len(mu_values),N),dtype=np.complex_)

		for j in range(len(m)):
			if j % 2 == 0: 	# if j is an even number
				E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
				Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]

				evec1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]
				evecm1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]

			if j % 2 != 0:	# if j is an odd number
				E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
				Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]

				evec1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]
				evecm1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]

		"""
		Assigning values for 
		 - indices < m[0] : 			E1_cut = E0, Em1_cut = E0prime.
		 - odd/even indices > m[-1] : 	

		"""

		E1_cut[0:m[0]] = E0[0:m[0]]
		Em1_cut[0:m[0]] = E0prime[0:m[0]]

		evec1_cut[0:m[0]] = evec0[0:m[0],:]
		evecm1_cut[0:m[0]] = evec0prime[0:m[0],:]

		E1_cut[m[-1]:] = E0[m[-1]:]
		Em1_cut[m[-1]:] = E0prime[m[-1]:]

		evec1_cut[m[-1]:] = evec0[m[-1]:,:]
		evecm1_cut[m[-1]:] = evec0prime[m[-1]:,:]


		return E1_cut, Em1_cut, evec1_cut, evecm1_cut








def calc_E1_Em1_from_E0_E0prime_var(E0,E0prime,evec0,evec0prime,var_values,N,filename,savefiles,tol=2):
	"""
	Calculate electron/hole energies E1 and Em1 from cutting and pasting the energies E0 and E0prime whose values are in
			E0 >= 0, E0prime <= 0 for all values of the variable
	and whose period is half that of E1 and Em1. E1 and Em1 have max/min values
			max(E1 and Em1) = max(E0), min(E1 and Em1) = min(E0prime).

	Calculates the derivtive of the energies. When the derivative changes, that's where the cut needs to be done!

	Parameters
	----------
	E0 : 			lowest energy for the solved system, which in the topological phase corresponds to the positive values of the lowest energy pair of modes.

	E0prime : 		lowest absolute energy of the solved system, whose values are negative.

	evec0 : 		eigenvector with the energy E0.
					 - indexing : 	evec0[var_res,N]
					 				where N = 2*2*Ny in spinful system with p-h dofs and is 1D with length Ny. var_res is the number of points in the variable one is considering the energies/evecs a function of, which could be e.g. mu_res.

	evec0prime : 	eigenvector with the energy E0prime.
					 - indexing : 	same as with evec0.

	splineE0 : 					Fitted 'function' of E0. Used to find the derivative of E0.

	splineE0_derivative : 		derivative of fitted 'function' splineE0. Used to find the indices where the sign of the derivative changes, which corresponds to the indices where the energies E1 and Em1 cross.

	splineE0prime : 			equivalently to above.

	splineE0prime_derivative : 	equivalently to above.

	tol : 						tolerance in indices for which elements in sign_change_indices_E0max and sign_change_indices are considered to be the same. Necessary because they are calculated in different manners, and it thus turns out that when they are supposed to be the same, they vary slightly from each other. Example: for mu_res = 1500, tol=1 would do the job of excluding the points in sign_change_indices_E0max from sign_change_indices. Default value of tol is set to tol = 2.

	sign_change_indices : 		indices in var (e.g. mu_values) where the derivative of the energy E0 (or equivalently E0prime) changes.



	sign_change_indices_actual_E0max_withintol : 	Indices of sign_change_indices corresponding to the indices in sign_change_indices_E0max. They are calculated by cheching if the indices in sign_change_indices_E0max are inside sign_change_indices_actual_E0max_withintol within the tolerance set by tol.

	m : 											indexes in var (e.g. mu_values) for which the elements of sign_change_indices that are not equal to the elements of sign_change_indices_actual_E0max_withintol.

	Notes
	-----
	 - The sign of the derivative of the energy changes also in between the desired points where it changes. So we need to identify which points are the desired points and which to throw away. We can identify that the points we want are the ones where the energy is approximately zero, while the points to throw away are the ones where the absolute value of the energy is at its maxima.
	 - All elements of the energies E1 and Em1 that are less than m[0] and larger than m[-1] need to be set equal to the appropriate E0 and E0prime. The convention used here is that
	 			 - for indices < m[0] : 	E1 = E0prime, Em1 = E0.
	 Depending on whether len(m) is even or odd, we get that E0 and Em1 should be for indices > m[-1]:
	 			 - for indices > m[-1]:
		 			 - for len(m) even : 	E1 = E0	
		 			 - for len(m) odd : 	E1 = E0prime

	"""

	from scipy.interpolate import UnivariateSpline
	import numpy as np
	import scipy.signal
	from My_Costum_Packages_newest import misc
	import datetime

	print("%s: In calc.calc_E1_Em1_from_E0_E0prime_var()"%misc.round_time(datetime.datetime.now(),round_to=60))

	splineE0 = UnivariateSpline(var_values, E0, s=0)
	splineE0prime = UnivariateSpline(var_values, E0prime, s=0)

	splineE0_derivative = splineE0.derivative(n=1)
	splineE0prime_derivative = splineE0prime.derivative(n=1)

	sign_change_indices = []
	sign_change_indices_E0max = scipy.signal.argrelmax(np.abs(E0))
	sign_change_indices_E0max = np.array(sign_change_indices_E0max[0])

	for i in range(len(splineE0_derivative(var_values))-1):
		if np.sign(splineE0_derivative(var_values)[i]) != np.sign(splineE0_derivative(var_values)[i+1]):
			sign_change_indices.append(i+1)

	"""Getting the actual indices in terms of the indices of sign_change_indices that correspond to the indices in sign_change_indices_E0max:"""
	sign_change_indices_actual_E0max_withintol = []

	for maximum in sign_change_indices_E0max:
		for delta in range(tol):
			# print(" - delta, ", delta)
			if maximum+delta in sign_change_indices:
				# print("hello", maximum)
				sign_change_indices_actual_E0max_withintol.append(maximum+delta)
			elif maximum-delta in sign_change_indices:
				# print("hello", maximum)
				sign_change_indices_actual_E0max_withintol.append(maximum-delta)

	mask = np.isin(sign_change_indices,sign_change_indices_actual_E0max_withintol)

	m = []
	for mask_i,element in zip(mask,sign_change_indices):
		# if True, meaning different elements
		if mask_i == False:
			m.append(element)

	E1_cut = np.zeros(len(var_values))
	Em1_cut = np.zeros(len(var_values))

	evec1_cut = np.zeros((len(var_values),N),dtype=np.complex_)
	evecm1_cut = np.zeros((len(var_values),N),dtype=np.complex_)

	# if len(m) > 1:
	# for j in range(len(m)):
	# 	print("j is",j)
	# 	if j % 2 == 0: 	# if j is an even number
	# 		E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
	# 		Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]

	# 		evec1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]
	# 		evecm1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]

	# 	if j % 2 != 0:	# if j is an odd number
	# 		E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
	# 		Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]

	# 		evec1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]
	# 		evecm1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]

	"""
	Assigning values for 
	 - indices < m[0] : 			E1_cut = E0, Em1_cut = E0prime 	| len(m) is odd
	 								opposite 						| len(m) is even
	 - odd/even indices > m[-1] : 	E1_cut = E0prime, Em1_cut = E0 	| len(m) is odd
	 								opposite 						| len(m) is even
	"""
	print(" - m is", m)
	if len(m) % 2 != 0: # odd means first and last assignment of E1_cut and Em1_cut is going to be different
		print(np.shape(E1_cut),np.shape(E0prime),np.shape(evec0prime))
		E1_cut[0:m[0]] = E0prime[0:m[0]]
		Em1_cut[0:m[0]] = E0[0:m[0]]

		evec1_cut[0:m[0]] = evec0prime[0:m[0],:]
		evecm1_cut[0:m[0]] = evec0[0:m[0],:]

		E1_cut[m[-1]:] = E0[m[-1]:]
		Em1_cut[m[-1]:] = E0prime[m[-1]:]

		evec1_cut[m[-1]:] = evec0[m[-1]:,:]
		evecm1_cut[m[-1]:] = evec0prime[m[-1]:,:]
		if len(m) > 1:
			for j in range(0,len(m)-1):
				print(" - j is", j)
				if j % 2 != 0: 	# if after an odd number of crossings
					E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
					Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]

					evec1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]
					evecm1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]
				else: 			# if after an even number of crossings
					E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
					Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]

					evec1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]
					evecm1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]
	else:				# even menas first and last assignment of E1_cut and Em1_cut is going to be the same
		if len(m) > 1:
			E1_cut[0:m[0]] = E0prime[0:m[0]]
			Em1_cut[0:m[0]] = E0[0:m[0]]

			evec1_cut[0:m[0]] = evec0prime[0:m[0],:]
			evecm1_cut[0:m[0]] = evec0[0:m[0],:]

			E1_cut[m[-1]:] = E0prime[m[-1]:]
			Em1_cut[m[-1]:] = E0[m[-1]:]

			evec1_cut[m[-1]:] = evec0prime[m[-1]:,:]
			evecm1_cut[m[-1]:] = evec0[m[-1]:,:]

			for j in range(0,len(m)-1):
				if j % 2 != 0: 	# if after an odd number of crossings
					E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
					Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]

					evec1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]
					evecm1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]
				else: 			# if after an even number of crossings
					E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
					Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]

					evec1_cut[m[j]:m[j+1],:] = evec0[m[j]:m[j+1],:]
					evecm1_cut[m[j]:m[j+1],:] = evec0prime[m[j]:m[j+1],:]

		else:
			E1_cut = E0prime
			Em1_cut = E0

			evec1_cut = evec0prime
			evecm1_cut = evec0
	# elif len(m) == 1:
	# 	print("length of E0prime:", len(E0prime))

	# 	E1_cut[0:m[0]] = E0prime[0:m[0]]
	# 	Em1_cut[0:m[0]] = E0[0:m[0]]

	# 	evec1_cut[0:m[0]] = evec0prime[0:m[0],:]
	# 	evecm1_cut[0:m[0]] = evec0[0:m[0],:]

	# 	E1_cut[m[0]:] = E0prime[m[0]:]
	# 	Em1_cut[m[0]:] = E0[m[0]:]

	# 	evec1_cut[m[0]:] = evec0prime[m[0]:,:]
	# 	evecm1_cut[m[0]:] = evec0[m[0]:,:]

	if savefiles == True:
		np.save("E1_"+filename+".npy",E1_cut)
		np.save("Em1_"+filename+".npy",Em1_cut)
		np.save("evec1_"+filename+".npy",evec1_cut)
		np.save("evecm1_"+filename+".npy",evecm1_cut)

	return E1_cut, Em1_cut, evec1_cut, evecm1_cut




















def calc_E1_Em1_from_E0_E0prime_mu_(E0,E0prime,evec0,evec0prime,mu_values,N,tol=2):
	"""
	Calculate electron/hole energies E1 and Em1 from cutting and pasting the energies E0 and E0prime whose values are in
			E0 >= 0, E0prime <= 0
	and whose period is half that of E1 and Em1. E1 and Em1 have max/min values
			max(E1 and Em1) = max(E0), min(E1 and Em1) = min(E0prime).

	Calculates the derivtive of the energies. When the derivative changes, that's where the cut needs to be done!

	Parameters
	----------
	E0 : 			lowest energy for the solved system, which in the topological phase corresponds to the positive values of the lowest energy pair of modes.

	E0prime : 		lowest absolute energy of the solved system, whose values are negative.

	evec0 : 		eigenvector with the energy E0.
					 - indexing : 	evec0[var_res,N]
					 				where N = 2*2*Ny in spinful system with p-h dofs and is 1D with length Ny. var_res is the number of points in the variable one is considering the energies/evecs a function of, which could be e.g. mu_res.

	evec0prime : 	eigenvector with the energy E0prime.
					 - indexing : 	same as with evec0.

	tol : 			tolerance in indices for which elements in sign_change_indices_E0max and sign_change_indices are considered to be the same. Necessary because they are calculated in different manners, and it thus turns out that when they are supposed to be the same, they vary slightly from each other. Example: for mu_res = 1500, tol=1 would do the job of excluding the points in sign_change_indices_E0max from sign_change_indices. Default value of tol is set to tol = 2.

	sign_change_indices : 	indices in var (e.g. mu_values) where the derivative of the energy E0 (or equivalently E0prime) changes.

	splineE0 : 				Fitted 'function' of E0. Used to find the derivative of E0.

	splineE0_derivative : 	derivative of fitted 'function' splineE0. Used to find the indices where the sign of the derivative changes, which corresponds to the indices where the energies E1 and Em1 cross.

	splineE0prime : 		equivalently to above.

	splineE0prime_derivative : 	equivalently to above.

	m : 					indexes in var (e.g. mu_values) for which the elements of sign_change_indices that are not equal to the elements of sign_change_indices_E0max, within the tolerance tol.

	Notes
	-----
	 - The sign of the derivative of the energy changes also in between the desired points where it changes. So we need to identify which points are the desired points and which to throw away. We can identify that the points we want are the ones where the energy is approximately zero, while the points to throw away are the ones where the absolute value of the energy is at its maxima.
	 - All elements of the energies E1 and Em1 that are less than m[0] and larger than m[-1] need to be set equal to the appropriate E0 and E0prime. The convention used here is that
	 			 - for indices < m[0] : 	E1 = E0, Em1 = E0prime.
	 Depending on whether len(m) is even or odd, we get that E0 and Em1 should be for indices > m[-1]:
	 			 - for indices > m[-1]:
		 			 - for len(m) even : 	E1 = E0prime	
		 			 - for len(m) odd : 	Em1 = E0


	"""

	from scipy.interpolate import UnivariateSpline
	import numpy as np
	import scipy.signal
	from My_Costum_Packages_newest import misc
	import datetime

	print("%s: In calc.calc_E1_Em1_from_E0_E0prime_mu()"%misc.round_time(datetime.datetime.now(),round_to=60))

	splineE0 = UnivariateSpline(mu_values, E0, s=0)
	splineE0prime = UnivariateSpline(mu_values, E0prime, s=0)

	m = np.array(scipy.signal.argrelextrema(np.abs(E0)))[0]
	m_ = np.array(scipy.signal.argrelextrema(np.abs(E0)))[0]


	E1_cut = np.zeros(len(mu_values))
	Em1_cut = np.zeros(len(mu_values))

	evec1_cut = np.zeros((len(mu_values),N),dtype=np.complex_)
	evecm1_cut = np.zeros((len(mu_values),N),dtype=np.complex_)

	for j in range(len(m)-1):
		# print(j,m[j],m[j+1])
		if j % 2 == 0: 	# if j is an even number
			E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
			Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]

		if j % 2 != 0:	# if j is an odd number
			E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
			Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]


	for n in range(len(m_)-1):

		if n % 2 == 0: 	# if n is an even number

			evec1_cut[m_[n]:m_[n+1],:] = evec0[m_[n]:m_[n+1],:]
			evecm1_cut[m_[n]:m_[n+1],:] = evec0prime[m_[n]:m_[n+1],:]

		if n % 2 != 0:	# if n is an odd number

			evec1_cut[m_[n]:m_[n+1],:] = evec0prime[m_[n]:m_[n+1],:]
			evecm1_cut[m_[n]:m_[n+1],:] = evec0[m_[n]:m_[n+1],:]
	"""
	Assigning values for 
	 - indices < m[0] : 			E1_cut = E0, Em1_cut = E0prime.
	 - odd/even indices > m[-1] : 	

	"""

	E1_cut[m[-1]:] = E0[m[-1]:]
	Em1_cut[m[-1]:] = E0prime[m[-1]:]

	evec1_cut[m_[-1]:] = evec0[m_[-1]:,:]
	evecm1_cut[m_[-1]:] = evec0prime[m_[-1]:,:]



	return E1_cut, Em1_cut, evec1_cut, evecm1_cut








def g0_12_asymm(q_L, q_R, omega, En, Gamma_L, Gamma_R, u_L, u_R, v_L, v_R):
	"""
	Calculate asymmetric nonlocal conductance for the single level approximation in the higgin4 note.

	Input Parameters
	----------------
	q_a : 			q_a = Gamma_a(abs(u_a)^2-abs(v_a)^2)

	omega : 		float, array, depending on desired output type.
					The variable in units of energy, meV, that g0_12_asymm is a function of. Is e.g. different from lambda = mu, Ez (e.g.). 
					omega is e.g. biasenergies!

	En : 			Eigenenergy of current eigenstate of consideration.


	Gamma_a : 		Gamma_a = pi*nu*abs(t_a)^2
					where
					t_a = t_pincher in lead number a.

	u_a : 			Bugoliubov u-factor of current eigenstate in lead a.	

	v_a : 			Similarly.


	Calculated Parameters
	---------------------
	q : 			q = sqrt(q_L^2+q_R^2)

	n : 			n = n_L + n_R =def n_0 + n_1

	n_a : 			n_a = Gamma_a(abs(u_a)^2+abs(v_a)^2)

	Notes
	-----
	 - a = "alpha" = L or R in note = first (0th) or second (1st) lead in the current model here.

	"""
	n_L = Gamma_L*(np.abs(u_L)**2+np.abs(v_L)**2)
	n_R = Gamma_R*(np.abs(u_R)**2+np.abs(v_R)**2)
	n = n_L + n_R
	q = np.sqrt(q_L**2+q_R**2)

	g0numerator = 2*En*n_R-32*q_L*Gamma_L*Gamma_R*np.imag(u_L*np.conj(u_R)*np.conj(v_L)*v_R)

	if type(omega) == float:
		g0denominator = (En**2+q**2+2*n_L*n_R-8*Gamma_L*Gamma_R*np.real(u_L*np.conj(u_R)*np.conj(v_L)*v_R)-omega**2)**2 + 4*n**2*omega**2
	elif len(omega) > 1:
		g0denominator = []
		for omegai in omega:
			g0denominator.append((En**2+q**2+2*n_L*n_R-8*Gamma_L*Gamma_R*np.real(u_L*np.conj(u_R)*np.conj(v_L)*v_R)-omegai**2)**2 + 4*n**2*omegai**2)

	g0denominator = np.array(g0denominator)
	# print(np.shape(g0denominator)) is (201, 200)

	g0_12_asym = [4*q_L*omega[i]*(g0numerator/g0denominator[i,:]) for i in range(len(omega))]

	return np.array(g0_12_asym)


def sys_ldos(sys, energy):
    prop_modes = sys.modes(energy=energy)[0]

    ldos = np.zeros((prop_modes.wave_functions.shape[0],), float)
    for i in np.arange(prop_modes.wave_functions.shape[1]):
        ldos += abs(prop_modes.wave_functions[:, i])**2

    return ldos/2.0/np.pi


def lead_ldos(lead, energy):
    prop_modes = lead.modes(energy=energy)[0]

    ldos = np.zeros((prop_modes.wave_functions.shape[0],), float)
    for i in np.arange(prop_modes.wave_functions.shape[1]):
        ldos += abs(prop_modes.wave_functions[:, i])**2

    return ldos/2.0/np.pi


def nu_1D_lead_directly(par, energy):
	nu = []
	lead = make.make_lead_onesite_wide_1D(par) 		# ??? Not the most general way to do this!
	lead = lead.finalized()

	for energyi in energy:
		nu.append(lead_ldos(lead, energyi))
	nu = np.array(nu)
	return nu

def nu_sys_directly(sys_, p, energies):
	"""
	General calculation of density of states for a finite system, sys_, being a finalized kwant.Builder() object, sys_.
	"""
	nu = []
	for energy in energies:
		# nu.append(sys_ldos(sys_, energy))
		nu.append(kwant.ldos(sys_, energy, args=[p]))
	nu = np.array(nu)
	return nu


def calc_g0_12_asymm(sys, par, ppar, lambda_, t_L, t_R, omega, En, evec_n,filename):
	"""
	Loops over lambda variable and calculates g0_12_asymm for the biasenergy (float) or biasenergies (array) given as input omega, using the g0_12_asymm function.

	Parameters
	----------
	lambda_ : 			array.
						Variable that the eigenenergies/eigenstates are calculated as a function of. 
						E.g. mu or Ez.
	
	omega : 			float or array.
						Variable energy that each g0_12_asymm is calculated for / is a function of, for every lambda_.

	par : 				SimpleNamespace object.
						Contains parameters of the system.
						Used to get e.g. biasenergies and t_pincher.

	Calulated Parameters
	--------------------
	un : 				array.
						u-vector containing every u-value for the n'th state under consideration, for every lambda_-value / as a function of lambda_.

	nu : 				array.
						dos !!!???of system as function of lambda_.

	OPTIONALLY CALCULATED:	These are calculated if they are not already given in the par SimpleNamespace class.
	Gamma_L : 			array. 
						Gamma at lead L (1), specified directly in par.
						Could in future implementation make it a function of the density of states nu(lambda_) of the lead, making Gamma_L a function of lambda_, too. 
						Shape: len(lambda_)

	Gamma_R : 			array.
						Similarly.

	Notes
	-----
	 - "n" in En, un and vn symbolizes the index of the eigenstate / eigenenergy under consideration.
	 - ??? below: may want to somehow access the already made lead of the system sys. I haven't figured out how to do this yet.
	 - ??? below: Which of the four orbitals is nu of? sum of all?
	"""

	print(misc.round_time(datetime.datetime.now(),round_to=60),": in calc_g0_12_asymm")

	u2mv2, u2, v2, u, v, un_L, un_R, vn_L, vn_R = calc_u2mv2_var(evec_n=evec_n, var_values=lambda_, Ny=par.Ny, filename=filename, savefiles=True, write_u_v_var_LR=True)

	# !!! for now, only for symmetric pincher:
	par.t_L = par.t_pincher
	par.t_R = par.t_pincher

	if ppar.doubleres == False:
		omega = par.biasenergies 
	else:
		omega = par.biasenergies_2
	###################### UNNECCESARY: NOT USED ANYMORE ####################
	# ###### nu of lead as function of lambda_: #######
	# nu = nu_1D_lead_directly(par, lambda_)
	# plot.plot_nu_nuorb_vs_mu(nu, mu=lambda_, filename=filename)

	# # Summing over the orbitals in nu:
	# nu = np.sum(nu,axis=1)

	# ###### nu of sys as function of lambda_: ######
	# nusys = nu_sys_directly(sys, par, lambda_)

	# plot.plot_nu_nuorb_vs_mu(nusys, mu=lambda_, filename=filename)

	### !!!??? USING SYS DOS, NOT LEAD DOS ###
	# nu = np.sum(nusys,axis=1)
	#########################################################################

	# ??? 2*pi here instead of Karsten's pi convention since nu is divided by 2pi?
	if hasattr(par, 'Gamma_L'):
		pass
	else:
		Error("Gamma_L needs to be specified in par")
		# par.Gamma_L = np.pi*nu*np.abs(t_L)**2
	if hasattr(par, 'Gamma_R'):
		pass
	else:
		Error("Gamma_R needs to be specified in par")
		# par.Gamma_R = np.pi*nu*np.abs(t_R)**2

	par.q_L = par.Gamma_L*(np.abs(un_L)**2-np.abs(vn_L)**2)
	par.q_R = par.Gamma_R*(np.abs(un_R)**2-np.abs(vn_R)**2)


	g0_12_asym = g0_12_asymm(q_L=par.q_L, q_R=par.q_R, omega=omega, En=En, Gamma_L=par.Gamma_L, Gamma_R=par.Gamma_R, u_L=un_L, u_R=un_R, v_L=vn_L, v_R=vn_R)

	np.save("g0_12_asym_"+filename, g0_12_asym)
	print("%s: Saved file %s" %(misc.round_time(datetime.datetime.now(),round_to=60),"g0_12_asym_"+filename+".npy"))

	return g0_12_asym


def get_un_vn(p):
	"""
	Gets all values of u and v for the n'th eigenstate, as a function of the variable that e.g. the eigenenergy is a function of, named lambda_.
	

	Notes
	-----
	 - lambda_ is e.g. mu or Ez, and has the unit of meV.

	"""

	pass