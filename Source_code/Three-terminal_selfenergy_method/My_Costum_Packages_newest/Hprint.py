"""

PACKAGE FOR PRINTING SPECIFIC PARTS OF THE HAMILTONIAN.

"""

def H_nxy_as000110_print(sys,par,n_xystart,printtype):

	""" 
	Function printing particular blocks of the Hamiltonian depending in the input printtype.

	Parameters
	----------
	n_xystart = [to_sites,from_sites] for H00==H_n_xystart.
	printtype = off-diagonal,LMboundary, MRboundary or hopx

	Notes
	-----
	### ???: fix later: in print, use np.round(Hij,#decimals) where #decimals is an integer to get prettier output that is easier to read/compare.
	"""

	print("in Hnxy...")
	if printtype =='off-diagonal':
		"""From this, also the neighboring off-diagonal blocks are printed."""
		H00 = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[n_xystart[0]],from_sites=[n_xystart[1]])
		H00 = H00.tocsc() #'Return a copy of this matrix in Compressed Sparse Column format'
		H00 = H00.todense()
		print("H%.0f,%.0f for Ez=%g:" %(n_xystart[0],n_xystart[1],par.Ez)); print(H00)

		H01 = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[n_xystart[0]],from_sites=[n_xystart[1]+1])
		H01 = H01.tocsc() #'Return a copy of this matrix in Compressed Sparse Column format'
		H01 = H01.todense()
		print("H%.0f,%.0f for Ez=%g:" %(n_xystart[0],n_xystart[1]+1,par.Ez)); print(H01)

		H10 = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[n_xystart[0]+1],from_sites=[n_xystart[1]])
		H10 = H10.tocsc() #'Return a copy of this matrix in Compressed Sparse Column format'
		H10 = H10.todense()
		print("H%.0f,%.0f Ez=%g:" %(n_xystart[0]+1,n_xystart[1],par.Ez)); print(H10)

	elif printtype=="LMboundary":
		"""Accesses hamiltonian for last site of L and first site of M in order to check that Gamma(x) and phi(x) change as they should on the boundary."""
		LMboundary_left_Ny = (par.left[-1] + 1)*(par.Ny) - 1
		LMboundary_leftp1_0 = (par.left[-1] + 1)*(par.Ny)
		HL = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[LMboundary_left_Ny],from_sites=[LMboundary_left_Ny])
		HM = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[LMboundary_leftp1_0],from_sites=[LMboundary_leftp1_0])
		HLyym1 = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[LMboundary_left_Ny],from_sites=[LMboundary_left_Ny-1])
		HL = HL.todense(); HM = HM.todense(); HLyym1 = HLyym1.todense()
		print("H%.0f,%.0f Ez=%g:" %(LMboundary_left_Ny,LMboundary_left_Ny,par.Ez)); print(HL)
		print("H%.0f,%.0f Ez=%g:" %(LMboundary_leftp1_0,LMboundary_leftp1_0,par.Ez)); print(HM)
		print("H%.0f,%.0f Ez=%g:" %(LMboundary_left_Ny,LMboundary_left_Ny-1,par.Ez)); print(HLyym1)

	elif printtype=="MRboundary":
		"""Same as for LM, but now last site of M and first site of R."""
		MRboundary_right_Ny = (par.middle[-1] + 1)*(par.Ny) - 1
		MRboundary_rightp1_0 = (par.middle[-1] + 1)*(par.Ny)
		HM = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[MRboundary_right_Ny],from_sites=[MRboundary_right_Ny])
		HR = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[MRboundary_rightp1_0],from_sites=[MRboundary_rightp1_0])
		HM = HM.todense(); HR = HR.todense()
		print("H%.0f,%.0f Ez=%g:" %(MRboundary_right_Ny,MRboundary_right_Ny,par.Ez)); print(HM)
		print("H%.0f,%.0f Ez=%g:" %(MRboundary_rightp1_0,MRboundary_rightp1_0,par.Ez)); print(HR)

	elif printtype == "hopx":	# not generalized yet for Ny not equal to 160, which is the case for Nx=160.
		H0160 = sys.hamiltonian_submatrix(args=[par],sparse=True,to_sites=[0],from_sites=[160])
		H0160 = H0160.todense()
		print("H%.0f,%.0f Ez=%g:" %(0,160,par.Ez)); print(H0160)

def Hprint(counter_Hprint):

	"""	
	Printing parts of the Hamiltonian the first time main() is called if printtype specified:
		
	printtype
	---
	"off-diagonal": printing H00, H01 and H10 for given par.Ez.
	"LMboundary": printing last site of L and first site of M
	"MRboundary": printing last site of M and first site of R
	"""

	if counter_Hprint == 0:
		par.Ez=par.Ez_values[0]
		H_nxy_as000110_print(sys,par,n_xystart=[0,0],printtype="off-diagonal")
		
