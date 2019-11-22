"""

This script specifies program specific parameters /
what the program will do (ppar), 
physical parameters (par), and physical constants (const).
ppar, par, and const are SimpleNamespace objects. [*]

"""

"""
							PARAMETERS
							----------
const : 		
	Simplenamespace object containing relevant physical constants, 
	such as the speed of light, the electron mass and Planch's reduced constant.

"""

"""
par : 	SimpleNamespace object containing relevant physical parameters.

Notes
-----
 - Parameters special for conductance calculation:
	Beware that the resolution of the biasenergies affects the 
	numerical value of the conductance! 
	If the resolution is too small, one may not resolve the conductance peaks. 
	The max conudctance value saturates such that max conductance = 4 e^2/h 
	for some (high) resolution.
"""

"""
ppar : 			
	Simplenamespace object containing "program specific parameters". 
	Sorted into three categories, (i), (ii), and (iii).
    
    Notes
    -----
        A selection of parameters contained in these are given below. 
        [*] Please see comments thoughout this script for more details.

	(i) The closed system type / geometry :
	--------------------------------------
	- one_D : 	
		Boolean. 
		If True (False), a one(two)-dimensional system will be built 
		in make.make_system().

	- make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight : 		
		Boolean. 
		If defined as instance of ppar and given as True, make.make_system 
		will call 
			make.make_1D_N_2D_SC_system.
		The resultant system builder object is a one-dimensional semiconducting wire 
		in proximity to a 2D superconductor, both in a magnetic field 
		(with the corresponding Zeeman energy, p.Ez), and with Rashba spin-orbit coupling. 
		A superconducitng pincher is also constructed at the barrier between the 
		1D semiconductor and the 2D superconducting lead, given that
			pincher_SC_lead
		is True.

	(ii) The leads being attached to the closed system :
	---------------------------------------------------
	- pincher_SC_lead : 	
		Boolean. 
		If true, a pincher is added in make.make_1D_N_2D_SC_system. Construction of the 
		components of the graph that are at and connected to the interface sites (the
		sites in the superconductor that connect to the semiconductor, and implements the
		pincher), is made from onsite, hoppingx, and hoppingy_2D_superconductor_pincher. 
		They are the components of the Hamiltonian, as specified in make.make_system().

	(iii) How the program is run / what to calculate, read, and ignore :
	--------------------------------------------------------------------
	- generating_g012a : 	
		None (do nothing with g012a) or length-2 or length-3 list/tuple/array.
		If not equal to None, the ...
			... first element is a boolean: 
					True means calculate g012a, 
					False means its datalog file (which should have been generated 
					previously) is to be read and saved in the par.g0_12a variable;
			... second element is a boolean referring to whether or not g0_12a is 
			to be plotted agains biasenergies and the chosen independent variable 
			(which may have been given as a string with ppar.var_name).
			... (optional) third element is an arbitrary value. If the third element 
			is given, g012a is plotted, also with energies and Cooper charges on top of it.
			This representation may be useful in initial calculations, for checking
			that the calculation works. 
	# NAMES #
	 - ppar.var_name : 				str.
									Name of the variable that is varied (e.g. "mu" or "Ez")
	 - ppar.ppar.par_const_name : 	str.
									Name of the variable that is held constant (e.g. "Ez" 
									or "mu" if the ppar.var_name is "mu" or "Ez", respectively.)
	 - ppar.par_const : 			float.
									Numerical value of the constant.

"""


import numpy as np
class SimpleNamespace(object):
    """Contains parameters of the problem"""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def set_params(ppar):
    """Physical constants:"""
	# nm = nanometers
	# meV = milli electron volts
	# c = speed of light
	# J = Joules
	# s = seconds
	# m = meters
	# kg = kilograms
    const = SimpleNamespace(
		# Speed of light :
        c = 2.99792458e17,			# [c]=nm/s
		# Energy corresponding to the mass of an electron :
        m_e = 0.5109989461e9,		# [m_e]=meV (energy units)
                                    # When mass units is used in
                                    # this script (see e.g. par.m_star),
                                    # [mass] = meV/c^2.
		# Reduced Planch's constant : 
        hbar = 6.582119514e-13,  	# [hbar]=meV*s
		# Reduced Planch's constant in SI units : 
        hbar_SI = 1.0547e-34,  		# [hbar_SI]=J*s
		# Speed of light in SI units :
        c_SI = 2.99792458e8,  		# [c_SI]=m/s
		# Mass of an electron in SI units : 
        m_e_SI = 9.10938356e-31  	# [m_e_SI]=kg
    )

	# If the closed system graph is a one-dimensional system : 
    if ppar.one_D == True:
		# The global parameter class : 
        par = SimpleNamespace(
			# No regions ('middle' or 'right') with
			# sites in the x-direction : 
			LM=0, 
			LR=0, 
			# Length of the system in the y-direction :
			LY=1500.,	# [LY]=nm
			# Default value for the chemical potential.
			# May be changed further down in this script : 
			mu=0, 		# [mu]=meV
			# alpha * hbar * c : 
			# [alpha]=nm*meV
			# [hbar]= meV*s
			# [c]=nm/s
			# [alpha*hbar*c]
			# 	=nm*meV*meV*s*nm/s
			# 	=(nm*meV)^2			 
			alphahbar=1.42e-4*const.c*const.hbar,
			Gamma=180e-3, # [Gamma]=meV
			phiL=0, # [phi]=1 (unitless phase)
			phiR=0, 
			k=30	# integer number of lowest energy eigenenergies
					# to numerically diagonalize for the closed
					# system.
			) 
		# Number of sites in the x-direction : 
        par.Nx = 1		# is 1 since 1D system. Unitless.
		# Number of sites in the y-direction : 
        par.Ny = 60		
		# Lattice constant in the y-direction
		# is calculated from the length and the
		# number of sites in said direction : 
        par.ay = float(par.LY)/par.Ny
		# Square lattice implementation : 
        par.ax = par.ay
		# The only region ('left') in the x-direction where
		# sites are occupied : 
        par.LL = par.ax
	# If not 1D, implementing the 2D system described in Hell : 
    else:
        par = SimpleNamespace(
			LM=250., 
			LL=1.e3, 
			LR=1.e3, 
			LY=4.e3, 
			mu=0, 
			alphahbar=1.42e-4*const.c*const.hbar, 
			Gamma=180e-3, 	# Induced gap from self-energy due to
							# the proximity effect (see Hell)
			phiL=0, 
			phiR=0, 
			k=30
			)
        par.Nx = Nx 		# IMHERE <<--- should be a problem. 
        # Number of sites in y-direction, determined
        # by the aspect ratio, Nx/Ny = 1.625 in the
        # article by Hell et. al.
        par.Ny = int(round(float(par.Nx/1.625)))  # rounding off to nearest integer
        # Length of middle + left + right regions in the
        # x-direction, divided by number of sites in x-
        # direction, is the lattice constant in the x-
        # direction : 
        par.ax = float(par.LM+par.LL+par.LR)/par.Nx
        # Length along the wire in the y-direction, 
        # divided by the number of sites in said direction : 
        par.ay = float(par.LY)/par.Ny

    # The induced p-wave gap in the left and the right regions :
    par.GammaL = par.Gamma
    par.GammaR = par.Gamma 	                    # GammaL/R are used in make.gamma_phi_fn(site,p)
    # Effective mass of the semiconductor : 
    par.m_star = 0.023*const.m_e/(const.c)**2   # [m_star]=meV/c^2
    
    """ If the semiconducting wire Hamiltonian is not the effective, 
        proximitized wire : 
    IF two-terminal, NS system : 
       - Left lead: one-dimensional, semiconductor
       - Middle region (scattering region): one-dimensional semiconductor
       - Middle lead: two-dimensional superconductor
       - Right lead: not present
    OR two-terminal, NS system with a superconducting layer inside
    the scattering region, along the whole length of the semiconductor
    in the scattering region : 
       - Left lead: one-dimensional, semiconductor
       - Middle region (scattering region): one-dimensional semiconductor,
                                           connected to a two-dimensional
                                           superconductor along the whole length
                                           of the semiconductor
       - Middle lead: two-dimensional superconductor
       - Right lead: not present
    """
    if ppar.make_1DNLeft_1DN_2DSMiddle_No_NRight \
        or ppar.make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight:   # booleans, 
                                                            # specified in script containing 
                                                            # the main program
        # The effective mass in the semiconductor is typically ~ 20 times
        # the effective mass in the superconductor : 
        par.m_N = 20.*par.m_star
    
    """ If three-terminal system, for energies below the effective gap :
       - Left lead : semiconducting, one-dimensional
       - Middle (scattering region) : effective, rashba wire with applied field
                                       parallel to the wire
       - Right lead : semiconducting, one-dimensional
    Note
    ----
    This models energies below the gap. For energies above the gap, a
    superconducting lead should be attached along the length of the wire.
    In that way, the current measured in the normal leads is correctly
    correlated with the current going through the superconducting lead.
    (Above the gap, both Cooper pairs and supercurrent travels in the
    superconducting lead.
    Below the gap, supercurrent does not travel in said lead.)
    """
    if ppar.make_1D_Heff_LR_leads:
        # Effective mass is m_star, no distinction betwenn N and S mass 
        # because now the effective Hamiltonian is modelled in the scattering
        # region : 
        par.m_N = par.m_star
    par.tx = const.hbar**2/(2*par.m_N*par.ax**2)    # [t]=meV(=(meV*s)^2/(meV/(nm/s)^2*nm^2))
    par.ty = const.hbar**2/(2*par.m_N*par.ay**2)
    par.hbar = const.hbar


    """ 
    REAL SPACE DIMENSIONS 
    
    Notes
    -----
     - Setting symmetric or asymmetric device. 
     - If asymmetric: par.right=[]. This is obtained when defining par.LR = 0.
     - Lattice site numbering starts at 0 and ends at N-1, if there are N sites.
     """
    # Encoding all lattice sites, in arrays, for the left, middle and right regions:
    # Defined such that NL=NR, by making NM even/odd depending on Nx:
    par.left    = np.arange(0,  round(float(par.LL)/par.ax))
    par.middle  = np.arange(    round(float(par.LL)/par.ax),
                                round(float(par.LL+par.LM)/par.ax))
    par.right   = np.arange(    round(float(par.LL+par.LM)/par.ax),
                                round(float(par.LL+par.LM+par.LR)/par.ax))


    """ 
    DOMAINS AND RESOLUTIONS OF INDEPENDENT VARIABLES AND PARAMETERS 

    User action required
    --------------------
     - Comment (1) or (2), depending on which of the two variables
        (chemical potential, par.mu, or field, par.Ez) is constant, and which is
        varied.
    """
    # Domain for biasenergies : 
    par.bias_maxmin = 0.25                          # [meV]
    par.biasenergies = np.linspace(-par.bias_maxmin,
                                   par.bias_maxmin, par.biasEres) 	

    """ (1) Varying chemical potential, keeping field constant : """
    par.biasEres = 1001 # resolution in biasenergies
    par.Ezres = 1       # resolution in field
    par.mures = 1001    # resolution in chem. pot.
    # Ez is the constant, mu is the variable. 
    # In master thesis, V_Z is the independent variable (Ez = 2*V_Z).
    par.Ez = 0.8                                        # [Ez] = meV
    par.mu_values = np.linspace(-0.6, 0.6, par.mures)   # [mu] = meV
    # par.mu needs to be defined as the first value that the Hamiltonian
    # will be calculated for : 
    par.mu = par.mu_values[0]
    # An array for Ez_values needs to be defined for later :
    par.Ez_values = np.array([par.Ez])
    ppar.var_name, ppar.par_const_name, ppar.par_const = "mu", "Ez", par.Ez
    """ End of (1) """

    """ (2) Varying field, keeping chemical potential constant : """
    # par.Ezres = 1001
    # par.mures = 1
    # par.mu = 0.	
    # par.Ez_values = np.linspace(0,1,par.Ezres) 	# [Ez] = meV
    # par.Ez = par.Ez_values[0]
    # par.mu_values = np.array([par.mu])
    # ppar.var_name,ppar.par_const_name,ppar.par_const = "Ez","mu",par.mu
    """ End of (2) """

    """ 'DOUBLERES' : 															 
        used when wanting to load conductance files with different resolutions.
    """
    ppar.doubleres = True
    if ppar.doubleres == True:
        """(2) Vs. Ez : """
        # par.mures_2 = 1
        # par.Ezres_2 = 200
        # par.biasEres_2 = 201
        # par.mu_2 = 0.
        # par.mu_values_2 = np.array([par.mu_2])
        # par.Ez_values_2 = np.linspace(0,10,par.Ezres_2)
        # par.biasenergies_2 = np.linspace(-0.5,0.5,par.biasEres_2)
        """ End of (2) """
        """(1) Vs. mu : """
        par.mures_2 = 1001  # 251#IMHERE##1001
        par.Ezres_2 = 1
        par.biasEres_2 = 1001  # 251#IMHERE##1001
        par.Ez_2 = 0.8  # 1.6#IMHERE##0.8 			# didived by 2 gives V_Z in thesis!!!
        par.mu_values_2 = np.linspace(-.6, .6, par.mures_2)
        par.Ez_values_2 = np.array([par.Ez_2])
        par.biasenergies_2 = np.linspace(-par.bias_maxmin,
                                         par.bias_maxmin, par.biasEres_2)
        """ End of (1) """


    """ 
    PINCHER
    
    If par.pincher = True, the end sites in the scattering region
    have a different onsite potential from sites in the bulk. The
    tight-binding parameter at these end sites is
        par.t_pincher.
    In order for this site to act as a potential barrier, the
    coefficient par.t_p_coeff = par.t_pincher/par.ty needs to be
    larger than 1.
    """
    par.pincher = True
    par.t_p_coeff = 3
    par.t_pincher = par.t_p_coeff*par.ty

    """
    If the user wants to use par.t_pincher as the tight-binding
    parameter of
         - only the last site
         - only the hopping between the last and the next-to last
            sites, or
         - both of the two above,
    set
        ppar.pinchertype = "onsite", "hopping", or "onsite and hopping",
    respectively.
	"""
    ppar.pinchertype = "onsite"


    """
    SINGLE ANDREEV LEVEL MODEL
     - decide a value for the end 'weights', 
            pi * (lead DOS) * |t_alpha|^2
        for every value of the variable : 
    """
    par.Gamma_L = np.pi*np.abs(par.ty)**2*np.ones(par.mures)  # = 3.37 meV
    par.Gamma_R = np.pi*np.abs(par.ty)**2*np.ones(par.mures)



    """
	FOR Calculating/loading 
	 			E1, Em1, evec1, evecm1,
				u, v, u^2 - v^2 for the "n'th pair" of eigenvectors:

	Parameters
	----------
	n : 		list containing indices of two eigenstates that 
                have same absolute value of their energy for every 
                value of the variable.
    
    Notes
    -----
     - For conductance probe calculations, you want n to contain
        the indices for the two lowest energy modes. For instance,
        if par.k = 30 eigenstates are diagonalized, the states are
        sorted such that the 14'th and the 15'th elements contains
        the two lowest energy solutions. 
	"""
    # given par.k = 30 :
    ppar.n = [14,15]        # ppar category (iii)


    """ Tight-binding parameter for the semiconducting lead(s) : 
    Note : Band-width needs to be large compared to biasenergies considered. 
    ----
    """
    par.t_N_lead = 1.  # meV


    """
    If model with superconducting (SC) lead attached along the whole scattering region : 
        Need to specify some separate parameters for the superconducting lead.
    - Note : 
        Matching potentials : 
        2tx_SC + 2ty_SC + mu_SC = 0
        mu_N = 0
    """
    if ppar.make_1DNLeft_1DN_2DS_2DSMiddle_No_NRight \
        or ppar.make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight \
            or ppar.make_1D_NLeft_1D_S_Heff_No_NRight \
                or ppar.make_1DNLeft_1DN_2DSMiddle_No_NRight:

        # ppar.pincher_SC_lead = True

        # Tight-binding parameter of semiconductor, 
        # separate from that of the SC lead:
        par.t_N = par.ty

        par.t_SC_lead_coeff = 0.0001
        # Tight-binding parameter of SC lead, x-direction : 
        par.tx_SC_lead = par.t_SC_lead_coeff*par.ty
        # Tight-binding parameter of SC lead, y-direction : 
        par.ty_SC_lead = par.t_SC_lead_coeff*par.tx
        # Chemical potential in superconductor : 
        #     - Note : 
        #       The user may consider matching effective onsite
        #       potential of SC lead to that of the scattering region. 
        #           2tx_SC + 2ty_SC + mu_SC = 0
        #           mu_N = 0
        par.mu_SC_lead = -2.0*(par.tx_SC_lead+par.ty_SC_lead)
        # SC gap in SC lead : 
        par.Delta_SC_lead = 0.180 

        """ Defining Ez of S and S lead.
        Differs from that of N because effective g-factor is different. """ 
        # Ratio of g-factors, g(of SC)/g(of semiconductor) : 	
        par.g_S_over_g_N = -1./5
        # Effective field inside superconducting region(s) : 
        par.Ez_2D_S_lead = par.g_S_over_g_N*par.Ez
        par.Ez_2D_S = par.g_S_over_g_N*par.Ez

        # If 1D semiconductor connected to 2D SC lead, there is no layer
        # inside the scattering region along the wire, implementing a pincher.
        if ppar.make_1DNLeft_1DN_2DSMiddle_No_NRight == False:
            # Thus, setting pincher parameters to be the same as parameters in the SC lead :
            par.ty_SC_pincher_coeff = 1.
            par.ty_SC_pincher = par.ty_SC_pincher_coeff*par.ty
            par.tx_SC_pincher = par.ty_SC_pincher_coeff*par.tx
            par.mu_SC_pincher = par.mu_SC_lead
            par.Delta_pincher = par.Delta_SC_lead

    # If not two normal leads, i.e. not a three-terminal device, 
    # let ppar 'know' that only the local differential conductance
    # should be calculated : 
    if ppar.make_1D_NLeft_1D_N_2D_S_2D_SMiddle_No_NRight \
        or ppar.make_1D_NLeft_1D_S_Heff_No_NRight:
        # to be used in calc to only calculate G_11, and not G_12 as well : 
        ppar.oneNLead = True
    else:
        ppar.oneNLead = False


    ###


    return const, par, ppar
