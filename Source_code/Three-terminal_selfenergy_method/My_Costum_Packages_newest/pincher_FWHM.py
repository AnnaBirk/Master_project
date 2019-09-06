"""

PACKAGE FOR CALCULATING FWHM OF A FUNCTION f(x) WITH A PEAK AROUND x=0

"""


def FWHM(x,f):
	"""
	Finds FWHM of the function f by interpolating f'-max(f')/2 (where f' is f where its minimum value is set to zero) such that at its half max, this function crosses the y-axis. Thus, the difference in this interpolated function's roots is the FWHM.

	Parameters
	----------
	f : 	Function that has a symmetric peak at x=0 and which is only given for the precise domain in which the peak resides. Its minimum value is assumed to be at f[0] and f[-1]. Due to assuming f(x) is symmetric, f[0] and f[-1] should be equal within a small tolerance.
	x : 	Domain of the function f. Dimension of x and f has to be equal.

	Variables
	---------
	fmin : 		minimum of f(x), calculated as the average of f[0] and f[-1] in order to get a better precision for the FWHM.
	spline: 	used on f shifted such that fmin is at zero AND shifted such that at half of the original (f-fmin) value, it crosses the x-axis. Hence, the roots of spline gives the points one calculates the FWHM from.
	r1, r2 : 	defined assuming the (zero-bias) peak in f is symmetric.

	Notes
	-----
	 - s=0 means UnivariateSpline interpolates through all data points, {(x,f(x))}.
	 - The function does not pick out a certain peak if the function f has several peaks. Thus, it is important that the user gives a function f(x) that only has one peak.


	 ??? 	PROBLEM: fmin = find_1st_fmin(f[int(len(f)/2.):-1])
	 		SOLUTION: the peak of f(x)=conductance in my experiment is equal to the effective superconducting gap parameter, Gamma in par (parameter class). Thus, give already nice function f(x) for x in [-Gamma,Gamma] into the current function.
	"""
	from scipy.interpolate import UnivariateSpline

	fmin = (f[0]+f[-1])/2.
	if fmin != 0:
		print(" - fmin not equal to zero, calculated as: %g, where f[0] = %g, f[-1] = %g" %(fmin,f[0],f[-1]))

	spline = UnivariateSpline(x, (f-fmin)-np.max(f-fmin)/2. , s=0)
	roots = spline.roots()
	roots = np.abs(roots)
	min_root = np.min(roots)
	r2, r1 = min_root, -min_root

	return(r2-r1)


def calc_FWHM(x,x_domain,f):
	"""
	Uses FWHM() to get FWHM of a function f inside a domain in the x-values, x_domain, for a given t_pincher and writes t_pincher, FWHM to a file with the given filename. 

	Parameters
	----------
	x : 		list/np.array/tuple
	x_domain : 	list/array.
				start and stop value of x that the peak is within.
	f : 		list, np.array, tuple
				function for which FWHM is calculated.
	
	Variables calculated
	--------------------
	i0/1 : 	first/last index corresponding to x_domain.
	"""

	i0 = index_of_value(x,x_domain[0])
	i1 = index_of_value(x,x_domain[1])

	return FWHM(x[i0:i1],f[i0:i1])


def calc_write_FWHM_vs_t_pincher_values(par,G_values,biasenergies,biasenergies_domain,filename_w,pincher_coeff,pincher_coeffres):
	"""
	Write a file with FWHM values for a range of t_pincher values. The t_pincher values are not written themselves for ease of loading the data, but they are specified in the filename such that the user may manually generate the t_p_values when e.g. plotting the FWHM's.

	Parameters
	----------
	G : 			list, np.array, tuple.
					Conductance values vs. e.g. biasenergies. The specific conductance and variable it is a function of should be secifiedin filename_w.
	biasenergies : 
	biasenergies_domain : 	list/array.
							start and stop value of x that the peak is within.	
	filename_w : 	str.
					filename for FWHM data file.
	pincher_coeff : 	list, np.array, tuple.
					First and last values of t_pincher.
	pincher_coeffres : 		float.
					Step size between values inside pincher_coeff, assumed to be linearly spaced.
	"""
	
	par.FWHM = []

	for t_p,G in zip(pincher_coeff*par.ty,G_values):
		par.FWHM.append(calc_FWHM(biasenergies,biasenergies_domain,G))
	np.save(filename_w,FWHM)

