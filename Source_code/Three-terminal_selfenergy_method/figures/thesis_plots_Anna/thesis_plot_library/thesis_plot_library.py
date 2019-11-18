# thesis_plot_library

import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import SymLogNorm
from scipy.signal import argrelextrema

def calc_E1_Em1_from_E0_E0prime_var(E0,E0prime,var_values,tol=2,use_derivative=True,**kwargs):
    """
    Calculate electron/hole energies E1 and Em1 from cutting and pasting the energies E0 and E0prime whose values are in
            E0 >= 0, E0prime <= 0 for all values of the variable
    and whose period is half that of E1 and Em1. E1 and Em1 have max/min values
            max(E1 and Em1) = max(E0), min(E1 and Em1) = min(E0prime).

    Calculates the derivtive of the energies. When the derivative changes, that's where the cut needs to be done!

    Parameters
    ----------
    E0 :            lowest energy for the solved system, which in the topological phase corresponds to the positive values of the lowest energy pair of modes.

    E0prime :       lowest absolute energy of the solved system, whose values are negative.

    splineE0 :                  Fitted 'function' of E0. Used to find the derivative of E0.

    splineE0_derivative :       derivative of fitted 'function' splineE0. Used to find the indices where the sign of the derivative changes, which corresponds to the indices where the energies E1 and Em1 cross.

    splineE0prime :             equivalently to above.

    splineE0prime_derivative :  equivalently to above.

    tol :                       tolerance in indices for which elements in sign_change_indices_E0max and sign_change_indices are considered to be the same. Necessary because they are calculated in different manners, and it thus turns out that when they are supposed to be the same, they vary slightly from each other. Example: for mu_res = 1500, tol=1 would do the job of excluding the points in sign_change_indices_E0max from sign_change_indices. Default value of tol is set to tol = 2.

    sign_change_indices :       indices in var (e.g. mu_values) where the derivative of the energy E0 (or equivalently E0prime) changes.

    sign_change_indices_actual_E0max_withintol :    Indices of sign_change_indices corresponding to the indices in sign_change_indices_E0max. They are calculated by cheching if the indices in sign_change_indices_E0max are inside sign_change_indices_actual_E0max_withintol within the tolerance set by tol.

    m :                                             indexes in var (e.g. mu_values) for which the elements of sign_change_indices that are not equal to the elements of sign_change_indices_actual_E0max_withintol.

    **kwargs:
                 - evecs :                              list of two arrays.

                                                        evecs = [evec0,evec0prime]
                                                        needs to have same length/shape as energies.
                 
                                                        an array (e.g. u^2-v^2 or eigenfunction) that also needs to be cut at the same points where the enrgies are cut. 


    Notes
    -----
     - The sign of the derivative of the energy changes also in between the desired points where it changes. So we need to identify which points are the desired points and which to throw away. We can identify that the points we want are the ones where the energy is approximately zero, while the points to throw away are the ones where the absolute value of the energy is at its maxima.
     - All elements of the energies E1 and Em1 that are less than m[0] and larger than m[-1] need to be set equal to the appropriate E0 and E0prime. The convention used here is that
                 - for indices < m[0] :     E1 = E0prime, Em1 = E0.
     Depending on whether len(m) is even or odd, we get that E0 and Em1 should be for indices > m[-1]:
                 - for indices > m[-1]:
                     - for len(m) even :    E1 = E0 
                     - for len(m) odd :     E1 = E0prime

    """

    from scipy.interpolate import UnivariateSpline
    import scipy.signal

    print(" - In calc.calc_E1_Em1_from_E0_E0prime_var()")

    splineE0 = UnivariateSpline(var_values, E0, s=0)
    splineE0prime = UnivariateSpline(var_values, E0prime, s=0)

    splineE0_derivative = splineE0.derivative(n=1)
    splineE0prime_derivative = splineE0prime.derivative(n=1)

    sign_change_indices = []
    sign_change_indices_E0max = scipy.signal.argrelmax(np.abs(E0))
    sign_change_indices_E0max = np.array(sign_change_indices_E0max[0])



    if use_derivative == False:
        for i in range(len(splineE0(var_values))-1):
            if np.sign(splineE0(var_values)[i]) != np.sign(splineE0(var_values)[i+1]):
                sign_change_indices.append(i+1)

    else:
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

    if "evecs" in kwargs:
        [evec0, evec0prime] = kwargs["evecs"]

        evec1_cut = np.zeros(np.shape(evec0),dtype=np.complex_)
        evecm1_cut = np.zeros(np.shape(evec0),dtype=np.complex_)
    else:
        evec1_cut = np.zeros(len(var_values),dtype=np.complex_)
        evecm1_cut = np.zeros(len(var_values),dtype=np.complex_)
        [evec0, evec0prime] = [evec1_cut,evecm1_cut] #dummys not to be used

    if len(m) % 2 != 0: # odd means first and last assignment of E1_cut and Em1_cut is going to be different
        E1_cut[0:m[0]] = E0prime[0:m[0]]
        Em1_cut[0:m[0]] = E0[0:m[0]]

        evec1_cut[0:m[0]] = evec0prime[0:m[0]]
        evecm1_cut[0:m[0]] = evec0[0:m[0]]

        E1_cut[m[-1]:] = E0[m[-1]:]
        Em1_cut[m[-1]:] = E0prime[m[-1]:]
        evec1_cut[m[-1]:] = evec0[m[-1]:]
        evecm1_cut[m[-1]:] = evec0prime[m[-1]:]

        if len(m) > 1:
            for j in range(0,len(m)-1):
                # print(" - j is", j)
                if j % 2 != 0:  # if after an odd number of crossings
                    E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
                    Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
                    evec1_cut[m[j]:m[j+1]] = evec0prime[m[j]:m[j+1]]
                    evecm1_cut[m[j]:m[j+1]] = evec0[m[j]:m[j+1]]

                else:           # if after an even number of crossings
                    E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
                    Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
                    evec1_cut[m[j]:m[j+1]] = evec0[m[j]:m[j+1]]
                    evecm1_cut[m[j]:m[j+1]] = evec0prime[m[j]:m[j+1]]

    else:               # even menas first and last assignment of E1_cut and Em1_cut is going to be the same
        if len(m) > 1:
            E1_cut[0:m[0]] = E0prime[0:m[0]]
            Em1_cut[0:m[0]] = E0[0:m[0]]
            evec1_cut[0:m[0]] = evec0prime[0:m[0]]
            evecm1_cut[0:m[0]] = evec0[0:m[0]]

            E1_cut[m[-1]:] = E0prime[m[-1]:]
            Em1_cut[m[-1]:] = E0[m[-1]:]
            evec1_cut[m[-1]:] = evec0prime[m[-1]:]
            evecm1_cut[m[-1]:] = evec0[m[-1]:]

            for j in range(0,len(m)-1):
                if j % 2 != 0:  # if after an odd number of crossings
                    E1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
                    Em1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
                    evec1_cut[m[j]:m[j+1]] = evec0prime[m[j]:m[j+1]]
                    evecm1_cut[m[j]:m[j+1]] = evec0[m[j]:m[j+1]]

                else:           # if after an even number of crossings
                    E1_cut[m[j]:m[j+1]] = E0[m[j]:m[j+1]]
                    Em1_cut[m[j]:m[j+1]] = E0prime[m[j]:m[j+1]]
                    evec1_cut[m[j]:m[j+1]] = evec0[m[j]:m[j+1]]
                    evecm1_cut[m[j]:m[j+1]] = evec0prime[m[j]:m[j+1]]


        else:
            E1_cut = E0prime
            Em1_cut = E0
            evec1_cut = evec0prime
            evecm1_cut = evec0

    if "evecs" in kwargs:
        return E1_cut, Em1_cut, evec1_cut, evecm1_cut
    else:
        return E1_cut, Em1_cut


def calc_GS_GA_from_G0_Gm0(G0, Gm0):
    """
    Calculates the symmetric and antisymmetric conductances from the conductances for positive biasenergies, G0, and the condutcance for negative biasenergies, Gm0. 

    Parameters
    ----------
     - G0 and Gm0 :     arrays.
                        They are single-valued conductances for each var-value.

    Returns :           GS, GA 
    -------
    """

    GS = np.zeros(len(G0))
    GA = np.zeros(len(G0))

    for i in range(len(G0)):    # for every var-value index
        GA[i] = (G0[i] - Gm0[i])/2.
        GS[i] = (G0[i] + Gm0[i])/2.

    return GS, GA

def set_share_axes(axs, target=None, sharex=False, sharey=False):
    if target is None:
        target = axs.flat[0]
    # Manage share using grouper objects
    for ax in axs.flat:
        if sharex:
            target._shared_x_axes.join(target, ax)
        if sharey:
            target._shared_y_axes.join(target, ax)
    # Turn off x tick labels and offset text for all but the bottom row
    if sharex and axs.ndim > 1:
        for ax in axs[:-1,:].flat:
            ax.xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
            ax.xaxis.offsetText.set_visible(False)
    # Turn off y tick labels and offset text for all but the left most column
    if sharey and axs.ndim > 1:
        for ax in axs[:,1:].flat:
            ax.yaxis.set_tick_params(which='both', labelleft=False, labelright=False)
            ax.yaxis.offsetText.set_visible(False)

def plot_pcolormesh_4subplots_wzoom_a(  x1, y1, z1, x1_zoom, y1_zoom, z1_zoom,
                                        x2, y2, z2, x2_zoom, y2_zoom, z2_zoom, 
                                        colormap_log_2_value,
                                        colormap_limits=[False,0,0], colormap_norm=[False,0], 
                                        colormap_name=None, colormap_log=[False,0], 
                                        colormap_ticks = [False,0], font_size=10):
    plt.rc('text', usetex=True)
    font_size = 8
    plt.rcParams.update({'font.size': font_size})    

    fig, ax = plt.subplots(nrows=2, ncols=2, gridspec_kw = {'height_ratios':[3.5, 2]})
    set_share_axes(ax[0:1,0:2],sharey=True)

    if colormap_limits[0]==True:
        z_min, z_max = colormap_limits[1], colormap_limits[2]
    else:
        z_max = z1.max()
        z_min = z1.min()
    # matching limits for 2nd variable:
    z_min2, z_max2 = z_min, z_max

    if colormap_log[0] == True:
        CF1 = ax[0,0].pcolormesh(x1, y1, z1, vmin=z_min, vmax=z_max, cmap=colormap_name, norm=SymLogNorm(colormap_log[1],vmin=z_min, vmax=z_max), rasterized=True)

        CF1_zoom = ax[1,0].pcolormesh(x1_zoom, y1_zoom, z1_zoom, vmin=z_min, vmax=z_max, cmap=colormap_name, norm=SymLogNorm(colormap_log_2_value,vmin=z_min, vmax=z_max), rasterized=True)

        CF2 = ax[0,1].pcolormesh(x2, y2, z2, vmin=z_min2, vmax=z_max2, cmap=colormap_name, norm=SymLogNorm(colormap_log[1],vmin=z_min2,vmax=z_max2), rasterized=True)
    
        CF2_zoom = ax[1,1].pcolormesh(x2_zoom, y2_zoom, z2_zoom, vmin=z_min2, vmax=z_max2, cmap=colormap_name, norm=SymLogNorm(colormap_log_2_value,vmin=z_min2,vmax=z_max2), rasterized=True)
    else:
        CF1 = ax[0,0].pcolormesh(x1, y1, z1, vmin=z_min, vmax=z_max, cmap=colormap_name, norm=MidpointNormalize(midpoint=colormap_norm[1]), rasterized=True)

        CF1_zoom = ax[1,0].pcolormesh(x1_zoom, y1_zoom, z1_zoom, vmin=z_min, vmax=z_max, cmap=colormap_name, norm=MidpointNormalize(midpoint=colormap_norm[1]), rasterized=True)

        CF2 = ax[0,1].pcolormesh(x2, y2, z2, vmin=z_min2, vmax=z_max2, cmap=colormap_name, norm=MidpointNormalize(midpoint=colormap_norm[1]), rasterized=True)
    
        CF2_zoom = ax[1,1].pcolormesh(x2_zoom, y2_zoom, z2_zoom, vmin=z_min2, vmax=z_max2, cmap=colormap_name, norm=MidpointNormalize(midpoint=colormap_norm[1]), rasterized=True)

    """ COLORBAR: One for upper panel, one for lower """
    cbar_len = 0.61 #0.77#0.654
    cbar_width = 0.015
    cbar_x1 = 0.25#0.17
    cbar_y1 = 0.95
    ticks1 = [-0.1,-0.01,-0.001,0,0.01,0.001,0.1]
    ticks2 = [-0.1,-0.01,0,0.01,0.1]
    shifty = 0.02   # shift in lower panel to make space for xlabels in upper panels

    cbar1_ax = fig.add_axes([cbar_x1, cbar_y1, cbar_len, cbar_width])
    ticks1_ = [-0.1,-0.001,0,0.001,0.1]

    cbar1 = fig.colorbar(CF1, cax=cbar1_ax, orientation="horizontal", ticks=ticks1_)

    cbar2_ax = fig.add_axes([cbar_x1, cbar_y1-0.54\
                                                    -0.8*shifty, cbar_len, cbar_width])
    cbar2 = fig.colorbar(CF1_zoom, cax=cbar2_ax, orientation="horizontal", ticks=ticks2)
    
    cbar1_ax.tick_params(direction='in')
    cbar1_ax.text(-0.23,-2,r"$g_{LR}^{0,asym}$",fontsize=font_size+2)
    cbar2_ax.tick_params(direction='in')
    cbar2_ax.text(-0.23,-2,r"$g_{LR}^{0,asym}$",fontsize=font_size+2)
    cbar1.update_ticks()

    """ AXES AND ADJUSTMENTS: """
    ax[0,0].set_ylabel(r'$V\ [\mu eV]$', fontsize=font_size+2, va='center')
    ax[0,0].set_xlabel((r"$\mu\ [\mu eV]$"), fontsize=font_size+2)
    ax[0,1].set_xlabel(r"$V_Z\ [\mu eV]$", fontsize=font_size+2)
    ax[1,0].set_ylabel(r'$V\ [\mu eV]$', fontsize=font_size+2, va='center')
    ax[1,0].set_xlabel((r"$\mu\ [\mu eV]$"), fontsize=font_size+2)
    ax[1,1].set_xlabel(r"$V_Z\ [\mu eV]$", fontsize=font_size+2)
    ax[0,0].set_yticks([-200,0,200]) # may 'decouple' from data/cause problems later. beware.
    ax[0,1].set_yticks([-200,0,200])

    for i in [0,1]:
        for j in [0,1]:
            ax[i,j].tick_params(direction='in')

    # fig.set_size_inches(3.39, 3.1)  # letter-size 
    fig.set_size_inches(4, 4)
    fig.subplots_adjust(wspace=0.23, hspace=0.55+\
                                                    8*shifty, bottom = 0.13, left = 0.14, right=0.97, top=0.90)
    return fig, ax, cbar1


def calc_QLR(psi_components):
    wf_E0_L = psi_components[:,:,0,0]   # print(E0_index,E0prime_index) ---> 0,1
    wf_E0_R = psi_components[:,:,-1,0]

    wf_Em0_L = psi_components[:,:,0,1]   # print(E0_index,E0prime_index) ---> 0,1
    wf_Em0_R = psi_components[:,:,-1,1]

    ### for E0:
    # print(np.shape(wf_E0_R))
    u_up_E0_L = wf_E0_L[:,0]
    u_down_E0_L = wf_E0_L[:,1]
    v_up_E0_L = wf_E0_L[:,2]
    v_down_E0_L = -wf_E0_L[:,3]

    u_up_E0_R = wf_E0_R[:,0]
    u_down_E0_R = wf_E0_R[:,1]
    v_up_E0_R = wf_E0_R[:,2]
    v_down_E0_R = -wf_E0_R[:,3]

    u2_up_E0_L = np.multiply(np.conj(u_up_E0_L), u_up_E0_L)
    u2_down_E0_L = np.multiply(np.conj(u_down_E0_L), u_down_E0_L)
    v2_up_E0_L = np.multiply(np.conj(v_up_E0_L), v_up_E0_L)
    v2_down_E0_L = np.multiply(np.conj(v_down_E0_L), v_down_E0_L)

    u2_up_E0_R = np.multiply(np.conj(u_up_E0_R), u_up_E0_R)
    u2_down_E0_R = np.multiply(np.conj(u_down_E0_R), u_down_E0_R)
    v2_up_E0_R = np.multiply(np.conj(v_up_E0_R), v_up_E0_R)
    v2_down_E0_R = np.multiply(np.conj(v_down_E0_R), v_down_E0_R)

    q_E0_R = u2_up_E0_R + u2_down_E0_R - v2_up_E0_R - v2_down_E0_R

    n_E0_R = u2_up_E0_R + u2_down_E0_R + v2_up_E0_R + v2_down_E0_R

    q_E0_L = u2_up_E0_L + u2_down_E0_L - v2_up_E0_L - v2_down_E0_L

    n_E0_L = u2_up_E0_L + u2_down_E0_L + v2_up_E0_L + v2_down_E0_L

    Q_E0_R = q_E0_R / n_E0_R  # checked by plotting: has correct values! That is, not the extra factor of 2.

    Q_E0_L = q_E0_L / n_E0_L

    ### for Em0:
    u_up_Em0_L = wf_Em0_L[:,0]
    u_down_Em0_L = wf_Em0_L[:,1]
    v_up_Em0_L = wf_Em0_L[:,2]
    v_down_Em0_L = -wf_Em0_L[:,3]

    u_up_Em0_R = wf_Em0_R[:,0]
    u_down_Em0_R = wf_Em0_R[:,1]
    v_up_Em0_R = wf_Em0_R[:,2]
    v_down_Em0_R = -wf_Em0_R[:,3]

    u2_up_Em0_L = np.multiply(np.conj(u_up_Em0_L), u_up_Em0_L)
    u2_down_Em0_L = np.multiply(np.conj(u_down_Em0_L), u_down_Em0_L)
    v2_up_Em0_L = np.multiply(np.conj(v_up_Em0_L), v_up_Em0_L)
    v2_down_Em0_L = np.multiply(np.conj(v_down_Em0_L), v_down_Em0_L)

    u2_up_Em0_R = np.multiply(np.conj(u_up_Em0_R), u_up_Em0_R)
    u2_down_Em0_R = np.multiply(np.conj(u_down_Em0_R), u_down_Em0_R)
    v2_up_Em0_R = np.multiply(np.conj(v_up_Em0_R), v_up_Em0_R)
    v2_down_Em0_R = np.multiply(np.conj(v_down_Em0_R), v_down_Em0_R)

    q_Em0_R = u2_up_Em0_R + u2_down_Em0_R - v2_up_Em0_R - v2_down_Em0_R

    n_Em0_R = u2_up_Em0_R + u2_down_Em0_R + v2_up_Em0_R + v2_down_Em0_R

    q_Em0_L = u2_up_Em0_L + u2_down_Em0_L - v2_up_Em0_L - v2_down_Em0_L

    n_Em0_L = u2_up_Em0_L + u2_down_Em0_L + v2_up_Em0_L + v2_down_Em0_L

    Q_Em0_R = q_Em0_R / n_Em0_R  # checked by plotting: has correct values! That is, not the extra factor of 2.

    Q_Em0_L = q_Em0_L / n_Em0_L

    return  q_E0_R, n_E0_R, q_E0_L, n_E0_L, Q_E0_R, Q_E0_L,\
            q_Em0_R, n_Em0_R, q_Em0_L, n_Em0_L, Q_Em0_R, Q_Em0_L 


def calc_Q_dE_mu_b_phase_shift(mu, b, QR_mu, dEdmu, QR_b, dEdb):
    import scipy.signal 
    QR_mu_zeros = scipy.signal.argrelmin(np.abs(QR_mu))
    QR_b_zeros = scipy.signal.argrelmin(np.abs(QR_b))
    
    dEdmu_zeros = scipy.signal.argrelmin(np.abs(-dEdmu))
    dEdb_zeros = scipy.signal.argrelmin(np.abs(-dEdb))
    
    QR_mu_maxmin = scipy.signal.argrelmax(np.abs(QR_mu))
    QR_b_maxmin = scipy.signal.argrelmax(np.abs(QR_b))

    dEdmu_maxmin = scipy.signal.argrelmax(np.abs(-dEdmu))
    dEdb_maxmin = scipy.signal.argrelmax(np.abs(-dEdb))

    QR_mu_zeros_and_extrema = mu[np.sort(np.concatenate((QR_mu_zeros, QR_mu_maxmin), axis=1)[0])]
    dEdmu_zeros_and_extrema = mu[np.sort(np.concatenate((dEdmu_zeros, dEdmu_maxmin), axis=1)[0])]

    QR_b_zeros_and_extrema = b[np.sort(np.concatenate((QR_b_zeros, QR_b_maxmin), axis=1)[0])]
    dEdb_zeros_and_extrema = b[np.sort(np.concatenate((dEdb_zeros, dEdb_maxmin), axis=1)[0])]

    # plt.figure(1)
    # plt.plot(QR_b_zeros_and_extrema,'.-')
    # plt.plot(dEdb_zeros_and_extrema,'.-')
    # plt.grid("on",which="both")
    # plt.show()

    # phase_shifts_mu = np.abs(dEdmu_zeros_and_extrema - QR_mu_zeros_and_extrema)
    # ValueError: operands could not be broadcast together with shapes (17,) (15,) 


    phase_shifts_mu = np.abs(dEdmu_zeros_and_extrema - QR_mu_zeros_and_extrema)
    phase_shifts_b = np.abs(dEdb_zeros_and_extrema - QR_b_zeros_and_extrema)

    return  QR_mu_zeros_and_extrema, dEdmu_zeros_and_extrema, phase_shifts_mu,\
            QR_b_zeros_and_extrema, dEdb_zeros_and_extrema, phase_shifts_b

def calc_dEdvar(mu, b, E0_mu, E0_b):
    from scipy.interpolate import UnivariateSpline
    smoothing_factor = 0.02
    splineE0_mu = UnivariateSpline(mu, E0_mu, s=smoothing_factor)
    der_splineE0_mu = splineE0_mu.derivative()
    dEdmu = der_splineE0_mu(mu)

    splineE0_b = UnivariateSpline(b, E0_b, s=smoothing_factor)
    der_splineE0_b = splineE0_b.derivative()
    dEdb = der_splineE0_b(b)

    return dEdmu, dEdb

def plot_E0_gS_gA_and_QR_dE0_gSovergA_4subplots(    mu,E0_mu,gS_mu,gA_mu,QR_mu,
                                                    b,E0_b,gS_b,gA_b,QR_b):

    plt.rc('text', usetex=True)
    font_size = 8
    plt.rcParams.update({'font.size': font_size})
    fig, ax = plt.subplots(nrows=2, ncols=2)
    # set_share_axes(ax[0:1,0:2],sharey=True)
    set_share_axes(ax[1:2,0:2],sharey=True)


    """ PLOTTING ENERGIES"""
    color="gray"
    ax[0,0].plot(mu,E0_mu,'--',color=color,zorder=0)
    ax[0,1].plot(b,E0_b,'--',color=color,zorder=0)

    """ ENERGY AXIS SETTINGS """
    # ax[0,0].set_yticks([0], minor=True)
    ax[0,0].get_yaxis().set_visible(False)
    legend_00 = ax[0,0].legend([r"$E_0\ [a.u.]$"],frameon=False,loc="upper center")
    ax[0,1].get_yaxis().set_visible(False)
    ax[0,1].legend([r"$E_0\ [a.u.]$"],frameon=False,loc="upper center")

    range_scale_E = 1.3
    E0_mu_absmax = np.max(np.abs(E0_mu))
    ax[0,0].set_ylim(-E0_mu_absmax*range_scale_E,E0_mu_absmax*range_scale_E)
    E0_b_absmax = np.max(np.abs(E0_b))
    ax[0,1].set_ylim(-E0_b_absmax*range_scale_E,E0_b_absmax*range_scale_E)        

    """ 00 AXIS """
    color = "C0"
    ax_00_gA = ax[0,0].twinx()
    ax_00_gA.plot(mu,gA_mu*100,color=color,zorder=5)

    # ax_00_gA.spines['left'].set_color(color)
    ax_00_gA.set_ylabel(r"$g_{LR}^{0,asym}\ [10^{-2}],\ $",color=color)
    ax_00_gA.ticklabel_format(axis="y",style="scientific", scilimits=(4,-1),color=color)
    ax_00_gA.set_yticks([-1,0,1])#,which = "major")
    ax_00_gA.tick_params(
                            direction="in",
                            which = "major")
    ax_00_gA.yaxis.set_label_coords(-0.25,0.2)


    color = "C1"
    ax_00_gS = ax_00_gA.twinx()
    ax_00_gS.plot(mu,gS_mu,color=color,zorder=10)
    ax_00_gS.set_ylabel(r"$g_{LR}^{0,sym}\ [10^{-3}]$",color=color)
    ax_00_gS.yaxis.set_label_coords(-0.25,0.8
)
    ax_00_gS.tick_params(
                            axis="y",
                            which="both",
                            bottom=False,
                            top = False,
                            right = False,
                            left = False,
                            labelright = False
                            )

    """ 01 AXIS """
    color = "C1"
    ax_01_gS = ax[0,1].twinx()

    ax_01_gS.plot(b,gS_b,color=color)
    # ax_01_gS.set_ylabel(r"$g_{LR}^{0,sym}\ [10^{-3}]$",color=color))
    ax_01_gS.tick_params(
                            direction="in",
                            which = "major")
    # ax_01_gS.yaxis.set_label_coords(1.15,0.8)

    color = "C0"
    ax_01_gA = ax_01_gS.twinx()
    ax_01_gA.plot(b,gA_b*100,color=color)

    # ax_01_gA.set_ylabel(r"$g_{LR}^{0,asym}\ [10^{-2}], $",color=color)
    ax_01_gA.ticklabel_format(axis="y",style="scientific", scilimits=(4,-1),color=color)
    ax_01_gA.tick_params(
                            direction="in",
                            which = "major")
    # ax_01_gA.yaxis.set_label_coords(1.15,0.2)
    
    ax_01_gS.tick_params(
                            axis="y",
                            which="both",
                            bottom=False,
                            top = False,
                            right = False,
                            left = False,
                            labelright = False,
                            labelleft = False
                            )
    ax_01_gA.tick_params(
                        axis="y",
                        which="both",
                        bottom=False,
                        top = False,
                        right = True,
                        left = False,
                        labelright = True,
                        labelleft=False
                        )
    ax_01_gA.set_yticks([-5,0,5])


    """ AXIS 00 CONDUCTANCE AXES SETTINGS """
    range_scale = 1.1
    gS_mu_absmax = np.max(np.abs(gS_mu))
    gA_mu_absmax = np.max(np.abs(gA_mu))
    ax_00_gS.set_ylim(-gA_mu_absmax*0.1*range_scale,gA_mu_absmax*0.1*range_scale)
    ax_00_gA.set_ylim([-gA_mu_absmax*100*range_scale,gA_mu_absmax*100*range_scale])

    """ AXIS 01 CONDUCTANCE AXES SETTINGS """
    gS_b_absmax = np.max(np.abs(gS_b))
    gA_b_absmax = np.max(np.abs(gA_b))
    ax_01_gA.set_ylim(-gA_b_absmax*100*range_scale,gA_b_absmax*100*range_scale)
    ax_01_gS.set_ylim(-gA_b_absmax*0.1*range_scale,gA_b_absmax*0.1*range_scale)

    """ AXIS 10 CHARGE QR """
    sci_factor = 100
    color = 'k'
    # ax[1,0] = ax[1,0]#.twinx()
    ax[1,0].plot(mu,QR_mu*sci_factor,color=color)

    color = 'k'
    ax_QR_b = ax[1,1].twinx()   # Fixes axes scaling problem! / axis 11 and 10 charge being forced to be same for some reason
    # removing ax[1,1] ticks has to be done since this axis not directly used for plotting
    ax[1,1].tick_params(
                        axis="y",
                        which="both",
                        bottom=False,
                        top = False,
                        right = False,
                        left = False,
                        labelright = False,
                        labelleft=False)

    ax_QR_b.plot(b[:-1],QR_b[:-1]*sci_factor,color=color)

    """ AXIS 10 CHARGE AXES SETTINGS """
    QR_mu_absmax = np.max(np.abs(QR_mu))
    QR_b_absmax = np.max(np.abs(QR_b))
    
    
    ax_QR_b.set_ylim(-QR_b_absmax*range_scale*sci_factor,QR_b_absmax*range_scale*sci_factor)
    ax[1,0].set_ylim(-QR_mu_absmax*range_scale*sci_factor,QR_mu_absmax*range_scale*sci_factor)

    ax[1,0].tick_params(
                        axis="y",
                        which="both",
                        bottom=False,
                        top = False,
                        right = False,
                        left = True,
                        labelright = False,
                        labelleft=True,
                        direction="in"
                        )

    ax[1,0].set_ylabel(r"$Q_R\ [10^{-2}],\ $",color=color)
    ax[1,0].yaxis.set_label_coords(-0.155,0.15)


    """ AXIS 10 gS/gA """
    gS_over_A_mu = gS_mu/gA_mu
    color = "C2"
    alpha=0.5
    ax_10_gSA = ax[1,0].twinx()
    ax_10_gSA.plot(mu,gS_over_A_mu*sci_factor,color=color,alpha=alpha,linewidth=0.8)

    """ AXIS 11 gS/gA """
    gS_over_A_b = gS_b/gA_b
    color = "C2"
    ax_11_gSA = ax_QR_b.twinx()
    ax_11_gSA.plot(b[:-1],gS_over_A_b[:-1]*sci_factor,color=color,alpha=alpha,linewidth=0.8)

    """ AXIS 10 gS/gA AXIS SETTINGS """
    ax_10_gSA.set_ylim(-QR_mu_absmax*range_scale*sci_factor,QR_mu_absmax*range_scale*sci_factor)
    # ax_10_gSA.get_yaxis().set_visible(False)
    ax_10_gSA.tick_params( axis="y",
                            which="both",
                            left = False,
                            right = False,
                            bottom=False, 
                            top = False,
                            direction="in",
                            labelright=False,
                            labelleft=False)

    ax_10_gSA.set_ylabel(r"$\frac{g_{LR}^{0,sym}}{g_{LR}^{0,asym}}\ [10^{-2}]$",color=color)
    ax_10_gSA.yaxis.set_label_coords(-0.27,0.7)

    # """ AXIS 11 gS/gA AXIS SETTINGS """
    ax_11_gSA.set_ylim(-QR_b_absmax*range_scale*sci_factor,QR_b_absmax*range_scale*sci_factor)
    # ax_11_gSA.get_yaxis().set_visible(False)    
    ax_11_gSA.tick_params( axis="y",
                            which="both",
                            left = False,
                            right = False,
                            bottom=False, 
                            top = False,
                            direction="in",
                            labelright=False,
                            labelleft=False)



    """ GETTING dEdvar """
    dEdmu, dEdb = calc_dEdvar(mu,b,E0_mu,E0_b)

    """ AXIS 10 dEdmu """
    color = "gray"
    ax_dEdmu = ax[1,0].twinx()
    ax_dEdmu.plot(mu,-dEdmu,'--',color=color)

    """ AXIS 11 dEdb """
    ax_dEdb = ax_QR_b.twinx()
    ax_dEdb.plot(b,-dEdb,'--',color=color)

    """ GETTING SHIFT IN D ENERGY FROM CHARGE WHEN CHARGE HAS EXTREMA AND ZERO CROSSINGS"""
    QR_mu_zeros_and_extrema, dEdmu_zeros_and_extrema, phase_shifts_mu,\
    QR_b_zeros_and_extrema, dEdb_zeros_and_extrema, phase_shifts_b = calc_Q_dE_mu_b_phase_shift(mu, b, QR_mu, dEdmu, QR_b, dEdb)
    
    # ax_dEdmu.plot((QR_mu_zeros_and_extrema+dEdmu_zeros_and_extrema)/2., phase_shifts_mu/1000,'.-')
    # ax_dEdb.plot((QR_b_zeros_and_extrema+dEdb_zeros_and_extrema)/2., phase_shifts_b/1000.,'.-')

    """ AXIS 10 11 AXIS ADJUSTMENTS """
    dEdmu_absmax = np.max(np.abs(dEdmu))
    ax_dEdmu.set_ylim(-dEdmu_absmax*range_scale_E,dEdmu_absmax*range_scale_E)

    range_scale_E = 1.2
    dEdb_absmax = np.max(np.abs(dEdb))
    ax_dEdb.set_ylim(-dEdb_absmax*range_scale_E,dEdb_absmax*range_scale_E)

    ax_dEdmu.get_yaxis().set_visible(False)
    ax_dEdmu.legend([r"$-dE_0/d\mu\ [a.u.]$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.1,edgecolor="white")
    ax_dEdb.get_yaxis().set_visible(False)
    ax_dEdb.legend([r"$-dE_0/dV_Z\ [a.u.]$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.1,edgecolor="white")

    # linestyle = (0, (5, 1)) # densely dashed
    ax[0,0].plot(mu,np.zeros(len(mu)),linestyle=':',linewidth=0.5,color='k')
    ax[0,1].plot(b,np.zeros(len(b)),linestyle=':',linewidth=0.5,color='k')
    ax_dEdmu.plot(mu,np.zeros(len(mu)),linestyle=':',linewidth=0.5,color='k')
    ax_dEdb.plot(b,np.zeros(len(b)),linestyle=':',linewidth=0.5,color='k')

    """ AXES AND ADJUSTMENTS: """
    # ax[0,0].set_ylabel(r'$V\ [\mu eV]$', fontsize=font_size+2, va='center')

    ax[1,0].set_xlabel((r"$\mu\ [\mu eV]$"), fontsize=font_size+2)
    ax[1,1].set_xlabel(r"$V_Z\ [\mu eV]$", fontsize=font_size+2)

    for i in [0,1]:
        for j in [0,1]:
            ax[i,j].tick_params(direction='in')

    # fig.set_size_inches(3.39, 3.1)
    fig.set_size_inches(4., 4.)
    # fig.subplots_adjust(wspace=0.23, hspace=0.55, bottom = 0.13, left = 0.14, right=0.97)
    fig.subplots_adjust(left=0.13,wspace=0.1,right=0.93)

    ax_QR_b.tick_params(   axis="y",
                        which="both",
                        left = False,
                        right = True,
                        bottom=False, 
                        top = False,
                        direction="in",
                        labelright=True,
                        labelleft=False) # needs to be after cond. plot in order to update these tick params properly

    """ ANNOTATING SUBPLOTS """
    # xy=(0.02,0.908)
    xy=(0.018,0.89)

    # xy=(0.03,0.8)

    t_00 = ax[0,0].annotate("(a)", xy=xy, xycoords='axes fraction',fontsize=font_size+2)
    t_01 = ax[0,1].annotate("(b)", xy=xy, xycoords='axes fraction',fontsize=font_size+2)
    t_10 = ax_dEdmu.annotate("(c)", xy=xy, xycoords='axes fraction',fontsize=font_size+2)
    t_11 = ax_dEdb.annotate("(d)", xy=xy, xycoords='axes fraction',fontsize=font_size+2)

    t_10.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))
    t_11.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))

    t_00.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))
    t_01.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))


    data = [phase_shifts_mu, phase_shifts_b, QR_mu_zeros_and_extrema, dEdmu_zeros_and_extrema, QR_b_zeros_and_extrema, dEdb_zeros_and_extrema]
    return fig, ax, data 





def autocorr(x):
    result = np.correlate(x, x, mode='valid') # also have tried "same" --> local correlation
    # return result[result.size/2:]
    # return result[len(result)/2:]
    # print(result)
    # return np.average(result)
    # print(np.shape(result))
    return result


def P_l_function(L, rho_k_array):
    """
    Parameters
    ----------
     - L :              number of sites

     - rho_k_array :    |psi(site)|^2 = |u up(site)|^2 + |u down(site)|^2 + |v up(site)|^2 + |v down(site)|^2

    """

    # print("hello0", np.shape(rho_k_array))
    P_l = np.zeros(np.shape(rho_k_array))

    for l in range(L):
        P_l[:,l] = np.sum( rho_k_array[:,:l+1], axis=1)

    return P_l     


def A_L(L, P_l_array,axis=1):
    """
    Computes A_L which quantifies localization of wavefunction.

    Parameters
    ----------
     - L :          number of sites in system

     - P_l_array :  P_l for every l in 1, 2, ..., L, where

                    P_l = sum_k=1^l |psi(k)|^2

                    for a given energy.
    """
    return 1/float(L) * np.sum( np.exp( 2*np.pi * 1j * P_l_array ) , axis=axis)


import matplotlib
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    import matplotlib.ticker
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

def plot_wf2_vs_var_and_site_AN_Dtop(mu, b, sites_mu, wf2_mu, wf2_b, wf2_cbar_label, midzero=False, Dtop=False, Dtop2=[False], wf2_ave=False,lognorm=[False],x_cbarlabel=-1.,cbar_x1=0.35,cbar_len = 0.5):
    wf2_E0_mu_sum = wf2_mu 
    wf2_E0_b_sum = wf2_b 
    sites_mu_arr = np.arange(0,sites_mu)

    plt.rc('text', usetex=True)
    font_size = 8
    plt.rcParams.update({'font.size': font_size})
    fig, ax = plt.subplots(nrows=2, ncols=2, gridspec_kw = {'height_ratios':[3.3, 2]})
    fig.set_size_inches(4., 3.3)
    fig.subplots_adjust(top=0.84,wspace=0.07)
    set_share_axes(ax[1:2,0:2],sharey=True)
    set_share_axes(ax[0:1,0:2],sharey=True)



    # ax[0,0].pcolormesh(x1, y1, z1, vmin=z_min, vmax=z_max, cmap=colormap_name, norm=SymLogNorm(colormap_log[1],vmin=z_min, vmax=z_max), rasterized=True)


    # mesh_01 = ax[0,1].pcolormesh(b,sites,np.transpose(wf2_E0_b_sum),cmap="seismic")
    if midzero == True:
        wf2_mu_absmax = np.max(np.abs(wf2_E0_mu_sum))
        wf2_b_absmax = np.max(np.abs(wf2_E0_b_sum))

        if lognorm[0] == True:
            wf2_E0_mu_absmax = np.max(np.abs(wf2_E0_mu_sum))

            mesh_00 = ax[0,0].pcolormesh(mu,sites_mu_arr,np.transpose( wf2_E0_mu_sum),norm=SymLogNorm(lognorm[1],vmin=-wf2_E0_mu_absmax,vmax=wf2_E0_mu_absmax),cmap="seismic",rasterized=True)
            mesh_01 = ax[0,1].pcolormesh(b,sites_mu_arr,np.transpose(wf2_E0_b_sum),norm=SymLogNorm(lognorm[1],vmin=-wf2_E0_mu_absmax,vmax=wf2_E0_mu_absmax),cmap="seismic",rasterized=True)  
        else:
            mesh_00 = ax[0,0].pcolormesh(mu,sites_mu_arr,np.transpose( wf2_E0_mu_sum),cmap="seismic",rasterized=True)
            mesh_01 = ax[0,1].pcolormesh(b,sites_mu_arr,np.transpose(wf2_E0_b_sum),cmap="seismic",rasterized=True)
            mesh_00.norm = matplotlib.colors.Normalize(vmin=-wf2_mu_absmax,vmax=wf2_mu_absmax)
            mesh_01.norm = matplotlib.colors.Normalize(vmin=-wf2_b_absmax,vmax=wf2_b_absmax)
    else:
        mesh_00 = ax[0,0].pcolormesh(mu,sites_mu_arr,np.transpose( wf2_E0_mu_sum),cmap="seismic",rasterized=True)
        mesh_01 = ax[0,1].pcolormesh(b,sites_mu_arr,np.transpose(wf2_E0_b_sum),cmap="seismic",rasterized=True)
    """ Adding colorbar """
    cbar_width = 0.015
    cbar_y1 = 0.95#0.95

    cbar_00_ax = fig.add_axes([cbar_x1, cbar_y1, cbar_len, cbar_width])
    if lognorm[0] == True:
        cbar_00 = fig.colorbar(mesh_00, cax=cbar_00_ax, orientation="horizontal")#, format=OOMFormatter(-2, mathText=False))
    else:
        cbar_00 = fig.colorbar(mesh_00, cax=cbar_00_ax, orientation="horizontal", format=OOMFormatter(-2, mathText=False))

    cbar_00_ax.annotate(wf2_cbar_label,xy=(x_cbarlabel,-2), xycoords='axes fraction',fontsize=font_size+2)

    # cbar_00_ax.tick_params(
    #                         which="both",
    #                         direction="out",
    #                         axis="both"         )

    ax[0,0].set_ylabel(r"$site$")

    L = sites_mu
    rho_k_array_mu = wf2_E0_mu_sum 
    rho_k_array_b = wf2_E0_b_sum

    P_l_array_mu = P_l_function(L, rho_k_array_mu)
    P_l_array_b = P_l_function(L, rho_k_array_b)

    A_L_array_mu = np.abs(A_L(L,P_l_array_mu))
    A_L_array_b = np.abs(A_L(L,P_l_array_b))

    ax[1,0].plot(mu,A_L_array_mu,color='k')
    ax[1,1].plot(b,A_L_array_b,color='k')

    idx, maxAL = find_where_and_nearest(A_L_array_b, np.max(A_L_array_b))
    print(" max at b, A_L:", b[idx], maxAL)

    idx_mu, minAL_mu = find_where_and_nearest(A_L_array_mu[:500], np.min(A_L_array_mu[:500]))
    idx_mu_max, maxAL_mu = find_where_and_nearest(A_L_array_mu[:600], np.max(A_L_array_mu[:600]))

    idx_mu2 = argrelextrema(A_L_array_mu, np.less)[0][1]

    idx_mu_argrelmin2, argrelminAL_mu2 = find_where_and_nearest(A_L_array_mu, A_L_array_mu[idx_mu2])

    print(" min 1 at mu, A_L:", idx_mu, mu[idx_mu], minAL_mu)
    print(" max at mu, A_L:", idx_mu_max, mu[idx_mu_max], maxAL_mu)
    print(" min 2 at mu, A_L:", idx_mu_argrelmin2, mu[idx_mu_argrelmin2], argrelminAL_mu2)

    ax[1,0].set_xlabel(r"$\mu\ [\mu eV]$")
    ax[1,1].set_xlabel(r"$V_Z\ [\mu eV]$")

    ax[1,0].set_ylim(0,1.1)
    ax[1,1].set_ylim(0,1.1)


    if wf2_ave == True:
        """ PLOTTING SUM VALUE OF WF2 QUANTITY ALONG REAL SPACE AXIS, AS A FUNCTION OF THE VARIABLE """
        wf2_mu_ave =np.sum(wf2_E0_mu_sum,axis=1)#/len(wf2_E0_mu_sum[0,:])
        wf2_b_ave =np.sum(wf2_E0_b_sum,axis=1)#/len(wf2_E0_b_sum[0,:])

        ax_spin_10_legend = ax[1,0].twinx()
        ax_spin_11_legend = ax[1,1].twinx()

        ax_spin_10_legend.set_zorder(6)
        ax_spin_11_legend.set_zorder(6)

        ax_spin_10_legend.plot(mu,np.abs(wf2_mu_ave),zorder=1,color="C1",rasterized=True)
        ax_spin_11_legend.plot(b,np.abs(wf2_b_ave),zorder=1,color="C1",rasterized=True)

        ax[1,0].set_ylim(0,1.1)
        ax[1,1].set_ylim(0,1.1)

        ax_spin_10_legend.set_ylim(0,1.1)
        ax_spin_11_legend.set_ylim(0,1.1)

        ax[0,0].set_xlim(mu[0],600)

        ax[1,0].set_ylabel(r"$\left|A_N^{E_0}\right|\ [\ ]$")

        ax_spin_10_legend.legend([r"$|\sum_i{q_{i,\uparrow}^{E_0}-q_{i,\downarrow}^{E_0}}|\ [m^{-1}]$"],frameon=True,loc="lower center",borderpad=0.1,handletextpad=0.3,edgecolor="white")
        ax_spin_11_legend.legend([r"$|\sum_i{q_{i,\uparrow}^{E_0}-q_{i,\downarrow}^{E_0}}|\ [m^{-1}]$"],frameon=True,loc="lower center",borderpad=0.1,handletextpad=0.3,edgecolor="white")

        ax_spin_10_legend.get_yaxis().set_visible(False)
        ax_spin_11_legend.get_yaxis().set_visible(False)



        

    if Dtop == True:
        """ ALSO PLOTTING TOPOLOGICAL GAP FOR INFINITE WIRE (OREG 2010) """
        b_constant = 400     # mu eV
        mu_constant = 0      # mu eV
        Delta = 180 # mu eV
        Delta_top_mu = np.abs(b_constant - np.sqrt( Delta**2 + np.multiply(mu,mu) ))
        Delta_top_b = np.abs(b - np.sqrt( Delta**2 + mu_constant**2))

        color = "gray"
        ax[1,0].plot(mu, Delta_top_mu*1e-2,color=color, alpha=0.6)
        ax[1,1].plot(b, Delta_top_b*1e-2, color=color, alpha=0.6)

        ax[1,0].set_ylabel(r"$\left|A_N^{E_0}\right|\ [\ ],\ \Delta_{top}\ [10^{2}\, \mu eV]$")

        ax[1,0].set_ylim(0,2.5)
        ax[1,1].set_ylim(0,2.5)

    if Dtop2[0] == True:
        """ ALSO PLOTTING TOPOLOGICAL GAP found numerically (OREG 2010) """
        color = "gray"

        legend_10_ax = ax[1,0].twinx()
        legend_11_ax = ax[1,1].twinx()
        legend_10_ax.set_ylim(0,2.5)
        legend_11_ax.set_ylim(0,2.5)

        legend_10_ax.plot(Dtop2[1], Dtop2[2]*1e-2,color=color, alpha=0.6)
        legend_11_ax.plot(Dtop2[3], Dtop2[4]*1e-2, color=color, alpha=0.6)

        ax[1,0].set_ylabel(r"$\left|A_N^{E_0}\right|\ [\ ]$")

        ax[1,0].set_ylim(0,2.5)
        ax[1,1].set_ylim(0,2.5)


        legend_10_ax.get_yaxis().set_visible(False)
        legend_10_ax.legend([r"$ \Delta_{top}\ [10^2\,\mu eV]$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.3,edgecolor="white")

        legend_11_ax.get_yaxis().set_visible(False)
        legend_11_ax.legend([r"$ \Delta_{top}\ [10^2\,\mu eV]$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.3,edgecolor="white")


    """ ANNOTATING SUBPLOTS with (a)-(d)"""
    xy=(0.018,0.9)
    xy2 = (0.018,0.83)

    t_00 = ax[0,0].annotate("(a)", xy=xy, xycoords='axes fraction',fontsize=font_size+2)
    t_01 = ax[0,1].annotate("(b)", xy=xy, xycoords='axes fraction',fontsize=font_size+2)
    t_10 = ax[1,0].annotate("(c)", xy=xy2, xycoords='axes fraction',fontsize=font_size+2, zorder=51)
    t_11 = ax[1,1].annotate("(d)", xy=xy2, xycoords='axes fraction',fontsize=font_size+2)

    t_10.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))
    t_11.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))
    t_00.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))
    t_01.set_bbox(dict(facecolor="white", alpha=1, edgecolor="white",pad=0))

    fig.texts.append(ax[1,0].texts.pop())

    return fig, ax, cbar_00_ax, [P_l_array_mu, P_l_array_b, A_L_array_mu, A_L_array_b]

def find_where_and_nearest(array, value):
    """
    Returns index and array[index] where value is closest to an array element.
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def calc_wf_data(wf_mu, wf_b):

    sites_mu = len(wf_mu[0,0,:,0])
    sites_mu_arr = np.arange(0,sites_mu)

    wf_E0_mu = wf_mu[:,:,:,0]
    wf_E0_b = wf_b[:,:,:,0]

    wf_E0_mu_conj = np.conj(wf_E0_mu)
    wf_E0_b_conj = np.conj(wf_E0_b)

    wf_E0_mu_abs_squared = np.real(np.multiply(wf_E0_mu_conj,wf_E0_mu))
    wf_E0_b_abs_squared = np.real(np.multiply(wf_E0_b_conj,wf_E0_b))

    wf2_E0_mu_sum = np.sum(wf_E0_mu_abs_squared,axis=1)   # = |psi(site)|^2 = |u up(site)|^2 + |u down(site)|^2 + |v up(site)|^2 + |v down(site)|^2 (SUMS OVER E-H,SPIN)

    wf2_E0_b_sum = np.sum(wf_E0_b_abs_squared,axis=1)

    """ SPIN POLARIZATION """
    # = |u up(site)|^2 - |u down(site)|^2 - |v up(site)|^2 + |v down(site)|^2 
    u2_up_mu, u2_down_mu, v2_up_mu, v2_down_mu = wf_E0_mu_abs_squared[:,0,:], wf_E0_mu_abs_squared[:,1,:], wf_E0_mu_abs_squared[:,2,:], wf_E0_mu_abs_squared[:,3,:]

    u2_up_b, u2_down_b, v2_up_b, v2_down_b = wf_E0_b_abs_squared[:,0,:], wf_E0_b_abs_squared[:,1,:], wf_E0_b_abs_squared[:,2,:], wf_E0_b_abs_squared[:,3,:]

    spinpol_E0_mu = u2_up_mu - u2_down_mu - v2_up_mu + v2_down_mu
    spinpol_E0_b = u2_up_b - u2_down_b - v2_up_b + v2_down_b



    return wf2_E0_mu_sum, wf2_E0_b_sum, spinpol_E0_mu, spinpol_E0_b, u2_up_mu, u2_down_mu, v2_up_mu, v2_down_mu, u2_up_b, u2_down_b, v2_up_b, v2_down_b


def calc_plot_wf2_vs_var_and_site(mu, wf_mu, b, wf_b, Jeroen=[False,0], phase_data=[False,0], Q_realspace=[False,0], Dtop2=[False], level_cross_mu=[False]):

    wf2_E0_mu_sum, wf2_E0_b_sum, spinpol_E0_mu, spinpol_E0_b, u2_up_mu, u2_down_mu, v2_up_mu, v2_down_mu, u2_up_b, u2_down_b, v2_up_b, v2_down_b = calc_wf_data(wf_mu, wf_b)

    q_E0_R_mu, n_E0_R_mu, q_E0_L_mu, n_E0_L_mu, Q_E0_R_mu, Q_E0_L_mu,\
    q_Em0_R_mu, n_Em0_R_mu, q_Em0_L_mu, n_Em0_L_mu, Q_Em0_R_mu, Q_Em0_L_mu = calc_QLR(wf_mu[:])

    q_E0_R_b, n_E0_R_b, q_E0_L_b, n_E0_L_b, Q_E0_R_b, Q_E0_L_b,\
    q_Em0_R_b, n_Em0_R_b, q_Em0_L_b, n_Em0_L_b, Q_Em0_R_b, Q_Em0_L_b = calc_QLR(wf_b[:])



    """ NORMALIZING     / CHECK IF NORMALIZED : TRUE """
    # wf2_E0_mu_norm_const_var_i = np.sum(wf2_E0_mu_sum, axis = 1)
    # wf2_E0_b_norm_const_var_i = np.sum(wf2_E0_b_sum, axis=1)  

    # print(wf2_E0_b_norm_const_var_i) # returns 1 for all var values, so no need to normalize again as done here below:
    # for i in range(len(wf2_E0_mu_norm_const_var_i)):
    #     wf2_E0_mu_sum[i,:] = wf2_E0_mu_sum[i,:] / wf2_E0_mu_norm_const_var_i[i]
    #     wf2_E0_b_sum[i,:] = wf2_E0_b_sum[i,:] / wf2_E0_b_norm_const_var_i[i]

    ###### PLOTTING #######
    print("-------SPIN---------")
    sites_mu = len(wf_mu[0,0,:,0])
    if level_cross_mu[0] == True:
        level_cross_mu_index, level_cross_mu_array = find_where_and_nearest(mu, level_cross_mu[1])
        fig_spin, ax_spin, cbar_00_spin_ax, spin_P_l_A_L_mu_b = plot_wf2_vs_var_and_site_AN_Dtop(mu[:level_cross_mu_index], b, sites_mu, spinpol_E0_mu[:level_cross_mu_index], spinpol_E0_b, r"$q_{i,\uparrow}^{E_0}(\lambda)-q_{i,\downarrow}^{E_0}(\lambda)\ [m^{-1}]$", midzero=True, wf2_ave=True, x_cbarlabel=-0.7 , cbar_x1 = 0.4)#,lognorm=[True,0.005])

    else:
        fig_spin, ax_spin, cbar_00_spin_ax = plot_wf2_vs_var_and_site_AN_Dtop(mu, b, sites_mu, spinpol_E0_mu, spinpol_E0_b, r"$q_{\uparrow}-q_{\downarrow}\ [m^{-1}]$", midzero=True, wf2_ave=True, x_cbarlabel=-0.36 )#,lognorm=[True,0.005])

    ax_spin[1,0].axvspan(-600, -350, facecolor='b', alpha=0.1)
    ax_spin[1,0].axvspan(350, 600, facecolor='b', alpha=0.1)
    ax_spin[1,0].set_xlim(-600,600)

    ax_spin[1,1].axvspan(0, 200, facecolor='b', alpha=0.1)
    ax_spin[1,1].set_xlim(0,800)




    if level_cross_mu[0] == True:
        for axis in [ax_spin[0,0], ax_spin[1,0]]:
            axis.axvspan(level_cross_mu[1],mu[-1], facecolor='gray', alpha=0.97, zorder=5,hatch="//")


    if Q_realspace[0] == True:
        Q_E0_mu = u2_up_mu + u2_down_mu - (v2_up_mu + v2_down_mu)

        Q_E0_b = u2_up_b + u2_down_b - (v2_up_b + v2_down_b)

        print("-------Q---------")
        fig_Q, ax_Q, cbar_00_Q_ax, Q_P_l_A_L_mu_b = plot_wf2_vs_var_and_site_AN_Dtop( mu, b, sites_mu, Q_E0_mu, Q_E0_b, r"$Q_i^{E_0}(\lambda)\ [\ ]$", midzero=True, Dtop=False,lognorm=[True,0.001], x_cbarlabel=-0.4)

        ax_Q[1,0].axvspan(-600, -350, facecolor='b', alpha=0.1)
        ax_Q[1,0].axvspan(350, 600, facecolor='b', alpha=0.1)
        ax_Q[1,0].set_xlim(-600,600)

        ax_Q[1,1].axvspan(0, 200, facecolor='b', alpha=0.1)
        ax_Q[1,1].set_xlim(0,800)

        Q_E0_mu_global = np.sum(Q_E0_mu,axis=1)
        # ax_QR_mu = ax_Q[1,0].twinx()
        # # ax_QR_mu.plot(mu,Q_E0_mu[:,-1]*5000.) # Q_R
        # ax_QR_mu.plot(mu,Q_E0_mu_global,color="C1")
        # Q_mu_absmax = np.max(np.abs(Q_E0_mu_global))
        # ax_QR_mu.set_ylim(-Q_mu_absmax,Q_mu_absmax)
        # ax_QR_mu.legend([r"$\sum_i Q_i$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.1,edgecolor="white")
        ax_Q_dE_10 = ax_Q[1,0].twinx()
        dQ_mu = 1-np.abs(Q_E0_R_mu + Q_realspace[1])
        ax_Q_dE_10.plot(mu, dQ_mu, color="gray", alpha=0.6,zorder=5)
        ax_Q_dE_10.set_ylim(0,1.1)

        ax_Q_dE_11 = ax_Q[1,1].twinx()
        ax_Q_dE_11.plot(b, 1-np.abs(Q_E0_R_b + Q_realspace[2]), color="gray", alpha=0.6)
        ax_Q_dE_11.set_ylim(0,1.1)



        idx_mu1_dQ = argrelextrema(dQ_mu, np.less)[0][0]

        idx_mu_argrelmin1_dQ, argrelmindQ_mu1 = find_where_and_nearest(dQ_mu, dQ_mu[idx_mu1_dQ])

        print(" min 1 for 1-abs(QR-integrated(Q)) at mu, dQ..:", idx_mu_argrelmin1_dQ, mu[idx_mu_argrelmin1_dQ], argrelmindQ_mu1)


        Q_E0_b_global = np.sum(Q_E0_b,axis=1)
        # ax_QR_b = ax_Q[1,1].twinx()
        # # ax_QR_b.plot(b,Q_E0_b[:,-1]*5000.)  # QR
        # ax_QR_b.plot(b,Q_E0_b_global,color="C1")   # global Q
        # Q_b_absmax = np.max(np.abs(Q_E0_b_global))
        # ax_QR_b.set_ylim(-Q_b_absmax,Q_b_absmax)
        # ax_QR_b.legend([r"$\sum_i Q_i$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.1,edgecolor="white")

        # ax_QR_b.plot(b,-Q_realspace[2], '--', color="gray")

        # ax_Q[1,1].plot(b, np.abs(Q_E0_b_global + Q_realspace[2]), color="red", alpha=0.5)

        """ LABELLING """
        ax_Q[1,0].set_ylabel(r"$\left|A_N^{E_0}\right|\, [\ ]$")
        ax_Q_dE_10.legend([r"$1-(Q_R + dE/d\mu)\, [\ ]$"],frameon=True,loc="lower center",borderpad=0.1,handletextpad=0.3,edgecolor="white")
        ax_Q_dE_11.legend([r"$1-(Q_R + dE/dV_Z)\, [\ ]$"],frameon=True,loc="lower center",borderpad=0.1,handletextpad=0.3,edgecolor="white")

        ax_Q_dE_10.get_yaxis().set_visible(False)
        ax_Q_dE_11.get_yaxis().set_visible(False)


        
        if level_cross_mu[0] == True:
            for axis in [ax_Q[1,0], ax_Q[0,0]]:
                axis.axvspan(level_cross_mu[1],mu[-1], facecolor=level_cross_mu[2], alpha=0.97, zorder=50,hatch="//")




        # QR_mu_i, dEdmu_i, ph_mu_i,\
        # QR_b_i, dEdb_i, ph_b_i = calc_Q_dE_mu_b_phase_shift(mu[mu_domain:mu_domain_2], b[b_domain:b_domain_2], np.sum(Q_E0_mu,axis=1)[mu_domain:mu_domain_2], dEdmu[mu_domain:mu_domain_2], np.sum(Q_E0_b,axis=1)[b_domain:b_domain_2], dEdb[b_domain:b_domain_2])

        # ax_dEdmu = ax_Q[1,0].twinx()
        # ax_dEdmu.plot(mu,-Q_realspace[1], '--', color="gray")
        
        # ax_ph_mu_global = ax_Q[1,0].twinx()
        # color="C2"

        # ax_ph_mu_global.plot((QR_mu_i+dEdmu_i)/2.,ph_mu_i,'.-',color=color,alpha=0.6)
        # ph_mu_global_absmax = np.max(np.abs(ph_mu_i))*2
        # ax_ph_mu_global.set_ylim(0,ph_mu_global_absmax)


    print("-------PROB. DENSITY---------")
    fig, ax, cbar_00_ax, wf2_P_l_A_L_mu_b = plot_wf2_vs_var_and_site_AN_Dtop(mu, b, sites_mu, wf2_E0_mu_sum, wf2_E0_b_sum, r"$|\psi_i^{E_0}(\lambda)|^2\ [m^{-1}]$", Dtop=False,Dtop2=Dtop2,x_cbarlabel=-0.52)

    # fig_, ax_ = plt.subplots()
    # ax_.plot(Dtop2[2][:])
    # plt.show()
    # adsf
    """ when Dtop=False :"""
    if Dtop2[0] == False:
        ax[1,0].set_ylim(0,1.1)
        ax[1,1].set_ylim(0,1.1)
        ax[1,0].set_ylabel(r"$\left|A_N^{E_0}\right|\ [\ ]$")
    """"""
    


    if  phase_data[0] == True:
        QdE_ave_mu_i = phase_data[1]
        ph_mu_i = phase_data[2]

        QdE_ave_b_i = phase_data[3]
        ph_b_i = phase_data[4]

        ax_ph_mu = ax[1,0].twinx()
        color="C2"
        ax_ph_mu.plot(QdE_ave_mu_i,ph_mu_i,'.-',color=color, alpha=0.6)
        ax_ph_mu.legend([r"$\phi_{\frac{dE_0}{d\mu}, Q_R} \ [\mu eV]$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.1,edgecolor="white")

        ph_mu_absmax = np.max(np.abs(ph_mu_i))*2
        ax_ph_mu.set_ylim(0,ph_mu_absmax)

        ax_ph_b = ax[1,1].twinx()
        color="C2"
        ax_ph_b.plot(QdE_ave_b_i,ph_b_i,'.-',color=color, alpha=0.6)
        ax_ph_b.legend([r"$\phi_{\frac{dE_0}{d V_Z}, Q_R} \ [\mu eV]$"],frameon=True,loc="upper center",borderpad=0.1,handletextpad=0.1,edgecolor="white")

        ph_b_absmax = np.max(np.abs(ph_b_i))*2
        ax_ph_b.set_ylim(0,ph_b_absmax)





    """ COLORING OUT NON-TOPOLOGICAL REGION """
    ax[1,0].axvspan(-600, -357.2114, facecolor='b', alpha=0.1)
    ax[1,0].axvspan(357.2114, 600, facecolor='b', alpha=0.1)
    ax[1,0].set_xlim(-600,600)

    ax[1,1].axvspan(0, 180, facecolor='b', alpha=0.1)
    ax[1,1].set_xlim(0,800)


    if level_cross_mu[0] == True:
        for axis in [ax[0,0], ax[1,0]]:
            axis.axvspan(level_cross_mu[1],mu[-1], facecolor=level_cross_mu[2], alpha=0.97, zorder=50,hatch="//")


    plt.show()

    return fig, ax, fig_spin, fig_Q, wf2_P_l_A_L_mu_b, spin_P_l_A_L_mu_b, Q_P_l_A_L_mu_b









def plot_psi2_cut_cohlength_fit(wf_mu, wf_b): 

    ### OLD FUNCTION - 
    ########### EXTRACTING APPROXIMATE COHERENCE LENGTH IN MIDDLE OF TOP REGION ########################
    # print(np.shape(wf_mu))

    # wf_mu_cut = wf_mu[500,:,:,0]



    wf2_E0_mu_sum, wf2_E0_b_sum, spinpol_E0_mu, spinpol_E0_b, u2_up_mu, u2_down_mu, v2_up_mu, v2_down_mu, u2_up_b, u2_down_b, v2_up_b, v2_down_b = calc_wf_data(wf_mu, wf_b)

    # print(np.shape(wf2_E0_mu_sum)) (1000, 200)
    plt.figure(1)
    plt.plot(np.log(wf2_E0_mu_sum[500,:100]))

    x,y = np.arange(17,25),np.log(wf2_E0_mu_sum[500,17:25])
    x,y = np.arange(18,24),np.log(wf2_E0_mu_sum[500,18:24])

    # print(np.shape(x), np.shape(y)) (10,) (10,)
    B, logA = np.polyfit(x,y,1)
    # print(np.sum(wf2_E0_mu_sum[500,:])) --> 0.9999999999999998 approx 1; OK 

    """ SECOND FIT TRIAL: UPPER ESTIMATE OF COHERENCE LENGTH """
    max_idx, psi2_max = find_where_and_nearest(wf2_E0_mu_sum[500,:100],np.max(wf2_E0_mu_sum[500,:100]))

    x_max_mid, wf2_E0_mu_sum_max_mid = np.arange(max_idx,100), wf2_E0_mu_sum[500,max_idx:100]

    B2, logA2 = np.polyfit(x_max_mid, np.log(wf2_E0_mu_sum_max_mid), 1)

    print(B2, logA2)

    print(B,logA)
    plt.hold("on")
    plt.plot(x,logA+B*x)
    plt.plot(x_max_mid,logA2+B2*x_max_mid)

    """ FIGURE COSMETICS : """

    font_size = 10
    plt.rcParams.update({'font.size': font_size})

    plt.rc('text', usetex=True)

    fig, ax = plt.subplots()
    fig.set_size_inches(3.39, 3)

    plt.plot(wf2_E0_mu_sum[500,:100],color='k')

    """ PLOTTING FITTING : """
    x_fulldomain = np.arange(0,100)
    color="C0"
    ax_fit = ax.twinx()
    ax_fit.plot(x_fulldomain, np.exp(logA)*np.exp(B*x_fulldomain),'--', color=color)

    color = "C1"
    # ax_fit2 = ax.twinx()
    ax_fit.plot(x_fulldomain[max_idx:], np.array(np.exp(logA2)*np.exp(B2*x_fulldomain))[max_idx:], '--', color=color)

    ############# For own understanding:
    # xi_OregL = 30. # sites
    # ax_fit.plot(x_fulldomain, np.exp(logA)/12.*np.exp(-x_fulldomain/xi_OregL),'--', color='gray')
    #############

    wf2_mu_max = np.max(wf2_E0_mu_sum)*1.1
    wf2_mu_min = np.min(wf2_E0_mu_sum)

    ax.set_ylim([wf2_mu_min,wf2_mu_max])
    ax_fit.set_ylim([wf2_mu_min,wf2_mu_max])
    # ax_fit2.set_ylim([wf2_mu_min,wf2_mu_max])

    ax_fit.tick_params(
                        axis="y",
                        which="both",
                        bottom=False,
                        top = False,
                        right = False,
                        left = False,
                        labelright = False,
                        labelleft = False
                        )

    ax_fit.legend([r"$A e^{-x/\xi_A}$", r"$B e^{-x/\xi_B}$"],frameon=True,bbox_to_anchor=(0.4, 0.9),borderpad=0.1,handletextpad=0.1,edgecolor="white")
    leg = ax.legend([r"$A\approx %.1f\,m^{-1}$"%(np.exp(logA)) + "\n" + r"$\xi_A\approx%.0f\,sites$" %(-1/float(B)) + "\n" + \
                r"$B\approx %.2f\,m^{-1}$"%(np.exp(logA2)) + "\n" + r"$\xi_B\approx%.0f\,sites$" %(-1/float(B2))]\
                ,frameon=True,bbox_to_anchor=(0.4, 0.65),borderpad=0.1,handletextpad=0.1,edgecolor="white")
    leg.legendHandles[0].set_visible(False)

    ax.grid("on",axis="both", color='lightgray', linestyle='--', linewidth=1)
    plt.rc('text', usetex=True)

    ax.ticklabel_format(axis="y",style="scientific", scilimits=(4,-2))
    ax.set_xlabel(r"$site$")
    ax.set_ylabel(r"$|\psi_i(\mu=0\,\mu eV)|^2\ [m^{-1}]$")
    fig.tight_layout()
    plt.show()






def plot_psi2_cut_cohlength(wf_mu, wf_b, coherence_length_th, logA_th): 

    """ PLOTTING WAVEFUNCTION CUT WITH ALSO EXPONENTIAL FROM THEORETICAL COH LENGTH GIVEN AS INPUT """

    wf2_E0_mu_sum, wf2_E0_b_sum, spinpol_E0_mu, spinpol_E0_b, u2_up_mu, u2_down_mu, v2_up_mu, v2_down_mu, u2_up_b, u2_down_b, v2_up_b, v2_down_b = calc_wf_data(wf_mu, wf_b)
    x = np.arange(0,100)

    plt.figure(1)
    plt.plot(x, np.log(wf2_E0_mu_sum[500,:100]))

    B_th = -1./(coherence_length_th)
    plt.plot(x,logA_th+B_th*x)
    plt.show()


    """ FIGURE COSMETICS : """
    font_size = 10
    plt.rcParams.update({'font.size': font_size})

    plt.rc('text', usetex=True)

    fig, ax = plt.subplots()
    fig.set_size_inches(3.39, 3)

    plt.plot(wf2_E0_mu_sum[500,:100],color='k')

    """ PLOTTING FITTING : """
    color="C1"
    ax_fit = ax.twinx()
    ax_fit.plot(x, np.exp(logA_th)*np.exp(B_th*x),'--', color=color)

    wf2_mu_max = np.max(wf2_E0_mu_sum)*1.1
    wf2_mu_min = np.min(wf2_E0_mu_sum)

    ax.set_ylim([wf2_mu_min,wf2_mu_max])
    ax_fit.set_ylim([wf2_mu_min,wf2_mu_max])

    ax_fit.tick_params(
                        axis="y",
                        which="both",
                        bottom=False,
                        top = False,
                        right = False,
                        left = False,
                        labelright = False,
                        labelleft = False
                        )

    ax_fit.legend([r"$A e^{-x/\xi}$"],frameon=True,bbox_to_anchor=(0.4, 0.9),borderpad=0.1,handletextpad=0.1,edgecolor="white")
    leg = ax.legend([r"$A\approx %.1f\,m^{-1}$"%(np.exp(logA_th)) + "\n" + r"$\xi\approx%.0f\,sites$" %(-1/float(B_th)) ]\
                #+ "\n" + \
                #r"$B\approx %.2f\,m^{-1}$"%(np.exp(logA2)) + "\n" + r"$\xi_B\approx%.0f\,sites$" %(-1/float(B2))]\
                ,frameon=True,bbox_to_anchor=(0.4, 0.65),borderpad=0.1,handletextpad=0.1,edgecolor="white")
    leg.legendHandles[0].set_visible(False)

    ax.grid("on",axis="both", color='lightgray', linestyle='--', linewidth=1)
    plt.rc('text', usetex=True)

    ax.ticklabel_format(axis="y",style="scientific", scilimits=(4,-2))
    ax.set_xlabel(r"$site$")
    ax.set_ylabel(r"$|\psi_i(\mu=0\,\mu eV)|^2\ [m^{-1}]$")
    fig.tight_layout()
    plt.show()





def plot_Pl_and_cut_wf_spinpol_Q(wf2, spinpol, Q, Pl_wf2, Pl_spinpol, Pl_Q, cut_idx):
    fig_cuts, ax_cuts = plt.subplots(nrows=1,ncols=2)
    ax_cuts0_Pl = ax_cuts[0].twinx()

    import scipy.fftpack
    N = len(wf2[cut_idx,:])
    dt = 1.
    x = np.linspace(1,N*dt,N)
    y1 = wf2[cut_idx,:]
    Pl1 = Pl_wf2[cut_idx,:]
    y2 = spinpol[cut_idx,:]
    Pl2 = Pl_spinpol[cut_idx,:]
    y3 = Q[cut_idx,:]
    Pl3 = Pl_Q[cut_idx,:]

    y1f = scipy.fftpack.fft(y1)
    y2f = scipy.fftpack.fft(y2)
    y3f = scipy.fftpack.fft(y3)
    xf = np.linspace(0,1./(2*dt), N/2.)

    color="C0"
    ax_cuts[0].plot(x,y1, color=color)
    ax_cuts0_Pl.plot(x,Pl1,'k--', color=color)

    color="C1"
    ax_cuts[0].plot(x,y2, color=color)
    ax_cuts0_Pl.plot(x,Pl2,'k--', color=color)

    color="C2"
    ax_cuts[0].plot(x,y3, color=color)
    ax_cuts0_Pl.plot(x,Pl3,'k--', color=color)

    ax_cuts[0].legend(["wf2","spinpol","Q"])
    ax_cuts0_Pl.legend(["Pl wf2", "Pl spinpol", "Pl Q"])

    ax_cuts[1].plot(xf,2./N*np.abs(y1f[:N//2]))
    ax_cuts[1].plot(xf,2./N*np.abs(y2f[:N//2]))
    ax_cuts[1].plot(xf,2./N*np.abs(y3f[:N//2]))
    ax_cuts[1].set_title("Fourier transform")

    fig_cuts.subplots_adjust(wspace=0.33)

    plt.show()
