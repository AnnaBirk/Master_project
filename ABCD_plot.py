""" ABCD caluclation and plot as functions of energy, for different barriers. """

def ABCD_of_E_for_Z(ElD, EgD, Z, Delta):
	"""
	Calculating probabilities of Andreev reflection (A), normal reflection (B), transmission with branch-crossing (C), and transmission without branch-crossing (D).

	Parameters 
	----------
	 - ElD : 	energies less than Delta

	 - EgD : 	energies greater than Delta

	 - Z : 		barrier strength

	 - Delta : 	superconducting gap

	Notes
	-----
	 - xik : 	from BCS theory,

	 				xi_k^2 = E_k^2 + Delta_k^2

	 			where E_k is the quasiparticle excitation energy at momentum k and Delta_k is the superconducting gap.

	 - u2 : 	the electronic amplitude is u. 

	 - v2 : 	the hole amplitude is v. Total probability states

	 				|u|^2 + |v|^2 = 1.

	 - gam2 : 	gamma is a constant introduced in the probability amplitudes obtained in the BTK theory,

	 				gamma = u^2 + Z^2 (u^2 - v^2).

	 			We may identify u_k^2 - v_k^2 as the charge of the state with momentum k.

	 - Probabilities for E < Delta : 	is given by one expression. This is to ensure that the expressions for the probabilities are real numbers. Below this energy in the NS system considered here, no single particle transmission is possible into the superconductor and thus 

	 										C = D = 0		for E > Delta,

	 									and the total probability condition reads

	 										A + B = 1.

	 - Probabilities for E > Delta : 	is given by another expression. Above the gap, C and D are in general non-zero.

	"""
	A, B, C, D = np.zeros(len(ElD) + len(EgD)), np.zeros(len(ElD) + len(EgD)), np.zeros(len(ElD) + len(EgD)), np.zeros(len(ElD) + len(EgD))
	xik = np.sqrt(EgD**2 - Delta**2)

	u2 = 0.5 * (1 + xik / EgD)
	v2 = 1 - u2
	gam2 = (u2 + Z**2 * (u2 - v2))**2

	for i in range(len(ElD)):
		### E < D :
		
		A[i] = Delta**2/float(\
						ElD[i]**2 + (Delta**2-ElD[i]**2)*(1+2*Z**2)**2\
					)
		
		B[i] = 1-A[i]
		
		C[i], D[i] = 0, 0

	for j in range(len(EgD)):
		### E > D :
		
		j_ = len(ElD) + j

		A[j_] = u2[j]*v2[j]/float(gam2[j])

		B[j_] = (u2[j] - v2[j])**2 * Z**2 * (1 + Z**2) / float(gam2[j])

		C[j_] = u2[j] * (u2[j] - v2[j]) * (1 + Z**2) / float(gam2[j])

		D[j_] = v2[j] * (u2[j] - v2[j]) * Z**2 / float(gam2[j])

	return A, B, C, D

def plot_ABCD_for_Z(A, B, C, D, energies, size_inches=(2.5,2.), xlim=(0,3), legend=[True, r"$A$", r"$B$", r"$C$", r"$D$"]):

	plt.rc('text', usetex=True)
	font_size = 14
	plt.rcParams.update({'font.size': font_size})  
	fig, ax = plt.subplots(nrows=1, ncols=1)
	fig.set_size_inches(size_inches[0], size_inches[1])

	opacity = 0.9
	ax.plot(energies, A, '--', alpha=opacity)
	ax.plot(energies, B, alpha=opacity)
	ax.plot(energies, C, ':', alpha=opacity)
	ax.plot(energies, D, alpha=opacity)

	""" Cosmetics : """
	ax.set_xlabel(r"$ $")
	ax.set_yticks([0, 0.5, 1.0])
	ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
	ax.set_ylim([0,1.05])
	ax.set_xlim([xlim[0], xlim[1]])
	ax.set_xticks([0,1,2,3])
	ax.set_xticklabels(["0", r"$\Delta$"])
	ax.tick_params(
                            direction="in",
                            which = "major"
                          	)
	if legend[0] == True:
		ax.legend(legend[1:])

	fig.subplots_adjust(left=0.12,bottom = 0.11,right=0.99, top=0.99)

	return fig, ax


""" Main program : """

import numpy as np 

""" Setting particular parameters : """

Delta = 2	# units irrelevant, measuring energy in units of Delta.
N = 10001	# number of points to make energy arrays from
DD = 3		# factor of how many units - 1 of Delta the energies should take for energies above Delta
Z = 3.0		# unitless barrier strength
ElD = np.linspace(0,Delta,N)				# E < Delta array.
EgD = np.linspace(Delta,DD*Delta,DD*N)		# E > Delta array.

""" Calculating probabilities of Andreev reflection, Normal reflection, Transmission with branch-crossing, and transmission without branch-crossing, respectively. """

A, B, C, D = ABCD_of_E_for_Z(ElD, EgD, Z=Z, Delta=Delta)	

import matplotlib.pylab as plt 

""" Plotting A, B, C and D for the given barrier value Z : """
fig, ax = plot_ABCD_for_Z(A, B, C, D, energies = np.concatenate((ElD/Delta, EgD/Delta), axis=None), xlim = (0, DD), legend = [False])
fig.savefig("ABCD_Z=%.1f.pdf" % Z, format="pdf")
plt.show()

