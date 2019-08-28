
import numpy as np
import datetime
from My_Costum_Packages_newest import misc

import matplotlib.pyplot as plt 
import matplotlib.colors as colors
plt.rc('text', usetex=True)
plt.rcParams.update({'font.size': 18})


def plot_E_vs_mu(par):
	"""
	Old plotting function. Prettier plot from pplot_E_vs_var.
	"""
	f = plt.figure()
	plt.plot(par.mu_values,par.energies_mu)
	plt.ylabel("$Energies\, [meV]$",fontsize=16)
	plt.xlabel("$\mu\, [meV]$",fontsize=16)
	plt.savefig("Envsmu_" + par.filename + ".pdf")
	plt.hold("off")
	print("%s: saved figure %s" %(round_time(datetime.datetime.now(),round_to=60),("Envsmu_" + par.filename + ".pdf")))


def plot_E_vs_Ez(par,xlabel="$E_z\, [meV]$",ylabel="$E_n\, [meV]$",fontsize=16):
	
	"""
	Plots energies vs. Ez.

	Parameters
	----------
	par.asymmetric : boolean.
		If True, asymmetric par.energies is plotted. The value of Nx for the asymmetric device is appended in the figure title, and "_A" is added to the figure name. If False, symmetric par.energies is plotted, with no alteration to the figure title and filename.
	
	Notes
	-----
	Old plotting function. Prettier plot from pplot_E_vs_var.
	"""
	
	figtitle = str.format("$N_x$:%g,$N_y$:%g,$\mu$:%g,$\Gamma$:%g,$\phi_L$:%g,$\phi_R$:%g,$a_x$:%.0f,$a_y$:%.0f" %(par.Nx,par.Ny,par.mu,par.Gamma,par.phiL,par.phiR,par.ax,par.ay))
	filename_ = "EvsEz_" + par.filename 

	f = plt.figure()
	plt.plot(par.Ez_values,par.energies)
	plt.xlabel(xlabel,fontsize=fontsize)
	plt.ylabel(ylabel,fontsize=fontsize)
	
	if par.one_D == True:
		plt.title(figtitle,fontsize=fontsize)
	elif par.asymmetric == True:
		plt.title(figtitle + ";$N_x^A$:%g" %par.Nxasymmetric,fontsize=fontsize)
		filename += "_A"
	else:
		plt.title(figtitle,fontsize=16)

	"""Saving plot of energies vs. Ez:"""
	f.savefig(filename_ + ".pdf")
	f.hold("off")
	print("%s: saved figure %s" %(round_time(datetime.datetime.now(),round_to=60),(filename_ + ".pdf")))



def plot_map(dependent_var,x,y,filename,figtitle,xlabel,ylabel,cm):
	"""
	Plots general map of data of dependent_var vs independent_vars[0]=x_values, ---[1]=y_values, with filename, figtitle, labels and colormap cm.

	Parameters
	---------
		par: parameters. Take values of biasenergies used in calculation to plot against.
		independent_vars : 	list of two sets of independent variable values, to be plotted on the x and y axis.
		cm : 	colormap. E.g. "seismic" = red and blue, "PuBu" = pink/white to blue.
	"""
	fig,ax = plt.subplots()
	X,Y = np.meshgrid(x,y)
	dependent_var = np.transpose(dependent_var)

	simplecontour = ax.contourf(X,Y,dependent_var)
	pcm = ax.pcolormesh(X,Y,dependent_var,cmap = cm)
	fig.colorbar(pcm)

	ax.set_title(figtitle,fontsize=16)
	ax.set_xlabel(xlabel,fontsize=16)
	ax.set_ylabel(ylabel,fontsize=16)

	fig.set_size_inches(4, 3)
	plt.tight_layout()

	plt.savefig(filename)	
	print("%s: Saved figure"%misc.round_time(datetime.datetime.now(),round_to=60), filename)	


def plot_G_ij(par,scaling,G_ij,filename,figtitle,xlabel,ylabel):

	"""
		par: parameters. Take values to plot against.
	"""
	fig,ax = plt.subplots()

	G_ij=np.transpose(G_ij)

	simplecontour = ax.contourf(par.Ez_values,par.biasenergies,G_ij)


	pcm = ax.pcolormesh(par.Ez_values,par.biasenergies,G_ij,\
						norm = colors.SymLogNorm(\
						linthresh=0.0001,linscale=0.0001,vmin=G_ij.min(),vmax=G_ij.max()),cmap="seismic") #vmin=G_11.min(),vmax=G_11.max()))

	fig.colorbar(pcm,extend="both")
	# pcm = ax.pcolormesh(par.Ez_values,par.biasenergies,G_ij,\
	# 					norm = colors.PowerNorm(gamma=1./4.),cmap="seismic") #vmin=G_11.min(),vmax=G_11.max()))
	# fig.colorbar(pcm,extend="max")

	ax.set_title(figtitle,fontsize=16)
	ax.set_xlabel(xlabel,fontsize=16)
	ax.set_ylabel(ylabel,fontsize=16)
	colorbar = fig.colorbar(simplecontour)
	plt.tight_layout()
	
	plt.savefig(filename + "SymLogNorm.pdf")	#"_powerlaw_gamma_1over4.pdf") #
	print("%s: Saved figure"%round_time(datetime.datetime.now(),round_to=60), filename + "SymLogNorm.pdf")	#"_powerlaw_gamma_1over4.pdf") # 




def plot_G_ij_var(par,scaling,G_ij,var_values_str,filename,figtitle,xlabel,ylabel,cm,biasenergies_2=[],Ez_values_2=[],u2mv2_factor=1,**kwargs):#var_values_=[], E1_var=[], Em1_var=[], u2mv2_E1_var=[], u2mv2_Em1_var=[]):

	"""
	Parameters
	---------
		par: 		parameters. Take values of biasenergies used in calculation to plot against.
	
		scaling: 	"SymLogNorm", "PowerNorm" or "Linear".
	
		G_ij : 		array.
					conductance of choice.

		var_values_str: 	string name of variable_values to be plotted against, other than the biasenergies

		filename : 			str.

		figtitle : 			str.
							Conductance plotted.

		xlabel, ylabel : 	str.
							
		cm : 				str.
							colormap. E.g. "seismic", "PuBu"

		u2mv2_factor : 		int, float.
							Factor to scale the u2mv2-values such that their max value equals the max value of the energy of the corresponding state.
		**kwargs : 	
			legend : 			list, tuple.

			var_values_ : 		array.
								Variable values for which energies and u2mv2 are plotted against. Distinct from "par."var_values_str, enabeling one to plot with e.g. a higher resolution in the energies/u2mv2 than in the conductances, which may prove useful because valculation time of the conductances is highly limited by the chosen resolution in var_values_str.

			E1_var : 			array, list, tuple.
								Energy E_1, the lowest absolute value of energy. Partner to E_m1. Function of var_values.

			Em1_var : 			same.

			u2mv2_E1_var : 		array, list, tuple.
								u^2-v^2 for the E1 energy branch.
			
			biasenergies_2 : 	array.
								Other biasenergies for which the conductance is plotted against. Used when conductance resolution in bias is different from that of the energies.

			Ez_values_2 : 		Similarly.


	Example usage
	-------------
	 -> Plotting only conductance:
	 	--------------------------
	 	######## vs Ez:
		plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_11_Ez,var_values_str="Ez_values",filename="G_11_vsEbias1mu_seismic_"+par.filename,figtitle="$G_{11}\ [e^2/h]$",xlabel="$E_z\ [meV]$",ylabel="$E_{bias,1}$",cm="seismic")
		############# vs mu:
		# plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_11_mu,var_values_str="mu_values",filename="G_11_vsEbias1mu_seismic_"+par.filename,figtitle="$G_{11}\ [e^2/h]$",xlabel="$\mu\ [meV]$",ylabel="$E_{bias,1}$",cm="seismic")
	
	 -> Plotting both conductance, E_(m)1 and u2mv2_(m)1, with resolution 1 for Energies and Cooper charges and resolution 2 for conductance: 
	 	----------------------------------
		legend = ["$E_1\, [meV]$","$E_{-1}\, [meV]$","$(u^2_1-v^2_1)$","$(u^2_{-1}-v^2_{-1})$"]

		if ppar.doubleres == True:
			plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_12_A_Ez,var_values_str="Ez_values",filename="G_12_A_vsEbias1mu_seismic_"+ppar.filename,figtitle="$G_{12}^A\ [e^2/h]$",xlabel="$E_z\ [meV]$",ylabel="$E_{bias,1}$",cm="seismic",biasenergies_2=par.biasenergies_2,Ez_values_2=par.Ez_values_2,u2mv2_factor=np.abs(np.max(par.E1_Ez))/np.abs(np.max(par.u2mv2_E1_Ez)),var_values_=par.Ez_values, E1_var=par.E1_Ez, Em1_var=par.Em1_Ez, u2mv2_E1_var=par.u2mv2_E1_Ez, u2mv2_Em1_var=par.u2mv2_Em1_Ez,legend=legend)
		else:
			plot.plot_G_ij_var(par,scaling="Linear",G_ij=par.G_12_A_Ez,var_values_str="Ez_values",filename="G_12_A_vsEbias1mu_seismic_"+ppar.filename,figtitle="$G_{12}^A\ [e^2/h]$",xlabel="$E_z\ [meV]$",ylabel="$E_{bias,1}$",cm="seismic",u2mv2_factor=np.abs(np.max(par.E1_Ez))/np.abs(np.max(par.u2mv2_E1_Ez)),var_values_=par.Ez_values, E1_var=par.E1_Ez, Em1_var=par.Em1_Ez, u2mv2_E1_var=par.u2mv2_E1_Ez, u2mv2_Em1_var=par.u2mv2_Em1_Ez,legend=legend)

		import matplotlib.pylab as plt
		plt.show()


	 -> Plotting for analytical model:
	 	------------------------------
	 	legend = ["$E_1\, [meV]$","$E_{-1}\, [meV]$","$(u^2_1-v^2_1)$","$(u^2_{-1}-v^2_{-1})$"]

		plot.plot_G_ij_var(par,scaling="Linear",G_ij=np.transpose(g0_12_a),var_values_str="mu_values",filename="g012a_nus_vsEbimu_seis_"+ppar.filename,figtitle="$g^0_{12,a}\ [e^2/h]$",xlabel="$\mu\ [meV]$",ylabel="$E_{bias}$",cm="seismic",u2mv2_factor=np.abs(np.max(par.E1_mu))/np.abs(np.max(par.u2mv2_E1_mu)),var_values_=par.mu_values, E1_var=par.E1_mu, Em1_var=par.Em1_mu, u2mv2_E1_var=par.u2mv2_E1_mu, u2mv2_Em1_var=par.u2mv2_Em1_mu,legend=legend)

		import matplotlib.pylab as plt
		plt.show()
	"""


	fig,ax = plt.subplots()

	fig.set_size_inches(8, 4.3)

	G_ij=np.transpose(G_ij)

	cm = plt.cm.get_cmap(cm)
	cm.set_under("w")

	G_abs_max = np.max(np.abs(G_ij))

	if len(biasenergies_2) == 0:
		exec("global simplecontour; simplecontour = ax.contourf(par." + var_values_str + ",par.biasenergies,G_ij)")

	elif len(biasenergies_2) != 0:
		exec("global simplecontour; simplecontour = ax.contourf(par." + var_values_str + "_2,par.biasenergies_2,G_ij)")
	
	if scaling == "SymLogNorm":
		exec('global pcm; pcm = ax.pcolormesh(par.%s,par.biasenergies,G_ij,\
						norm = colors.SymLogNorm(\
						linthresh=0.0001,linscale=0.0001,vmin=G_ij.min(),vmax=G_ij.max()),cmap=cm)' %(var_values_str))

	elif scaling == "PowerNorm":
		exec('global pcm; pcm = ax.pcolormesh(par.%s,par.biasenergies,G_ij,\
							norm = colors.PowerNorm(gamma=1./4.),cmap=cm)' %(var_values_str))

	elif scaling == "Linear":
		if len(biasenergies_2) == 0:
			exec('global pcm; pcm = ax.pcolormesh(par.%s,par.biasenergies,G_ij,cmap=cm,vmin=-1*G_abs_max, vmax=G_abs_max)' %(var_values_str))
		elif len(biasenergies_2) != 0:
			exec('global pcm; pcm = ax.pcolormesh(par.%s_2,par.biasenergies_2,G_ij,cmap=cm,vmin=-1*G_abs_max, vmax=G_abs_max)' %(var_values_str))

		plt.gcf().subplots_adjust(bottom=0.15)
		plt.gcf().subplots_adjust(left=0.2)
	
	
	if len(kwargs) == 0:
		cb = fig.colorbar(pcm)#,extend='both')
		cb.formatter.set_powerlimits((-3, 4))
		cb.update_ticks()

	if len(kwargs)!=0: #len(var_values_) != 0 & len(E1_var) != 0 & len(Em1_var) != 0 & len(u2mv2_E1_var) != 0 & len(u2mv2_Em1_var) != 0 :
		print(" - adding energies and u2mv2 to conductance plot.")
		u2mv2_factor = np.round(u2mv2_factor,2)
		
		my_blue = plt.cm.get_cmap("PuBu",lut=100)(85)
		my_red = plt.cm.get_cmap("OrRd")(200)

		ax.plot(kwargs["var_values_"], kwargs["E1_var"], color=my_blue, linewidth=1.5, label=kwargs["legend"][0])
		ax.plot(kwargs["var_values_"], kwargs["Em1_var"], "--", color=my_red, linewidth=1.5,label=kwargs["legend"][1])
		ax.plot(kwargs["var_values_"], kwargs["u2mv2_E1_var"]*u2mv2_factor, color="black", linewidth=1.5,label=kwargs["legend"][2]+"$\cdot %s$"%u2mv2_factor)
		ax.plot(kwargs["var_values_"], kwargs["u2mv2_Em1_var"]*u2mv2_factor, "--", color="black", linewidth=1.5,label=kwargs["legend"][3]+"$\cdot %s$"%u2mv2_factor)
		
		
		box = ax.get_position()
		# making separate axis for colorbar
		cbaxes = fig.add_axes(ax)#[box.x0+box.width,box.y0+box.height*0.25,box.width*0.05,box.height*0.75]) 

		# cbaxes.set_position([box.x0+box.width,box.y0+box.height*0.25,box.width*0.05,box.height*0.75]) 
	
		cb = fig.colorbar(pcm, orientation='vertical',ax=cbaxes, shrink=0.75)
		# # Shrink current axis's height by 10% on the bottom

		ax.set_position([box.x0, box.y0 + box.height*0.25,
		                 box.width*0.75, box.height * 0.75])


		# Put a legend below current axis
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
		          fancybox=True, ncol=2,fontsize=12)


	ax.set_title(figtitle,fontsize=16)
	ax.set_xlabel(xlabel,fontsize=16)
	ax.set_ylabel(ylabel,fontsize=16)
	ax.ticklabel_format(style='sci',scilimits=(-3,4))

	plt.savefig(filename + "_" + scaling + ".pdf")	
	print("%s: Saved figure"%misc.round_time(datetime.datetime.now(),round_to=60), filename + scaling + ".pdf")	


def plot_G_ij_var_muconst(par,G_ij,index_of_varconst,filename,figtitle,xlabel,ylabel):
	"""
	Plots conductance G_ij against par.biasenergies for a constant value of the second variable it was calculated/mapped for.
	
	Parameters
	----------
	G_ij : 	np.array.
			Map of conductances against values of an arbitrary variable, call it 'var', and biasenergies. The indexes of its first axis indexes var, while on the second axis, they index biasenergies.

	"""
	fig,ax = plt.subplots()

	plt.plot(par.biasenergies,G_ij[index_of_varconst,:])

	ax.set_title(figtitle,fontsize=16)
	ax.set_xlabel(xlabel,fontsize=16)
	ax.set_ylabel(ylabel,fontsize=16)
	ax.ticklabel_format(style='sci',scilimits=(-3,4))

	plt.savefig(filename + "_" + ".pdf")	
	print("%s: Saved figure"%round_time(datetime.datetime.now(),round_to=60), filename + ".pdf")	

# OLD:
# def plot_G_1_G_2_varconst(par,G_1,G_2,index_of_varconst,filename,figtitle,xlabel,ylabel,legend,**kwargs):
# 	"""

# 	Parameters
# 	----------
# 	**kwargs : 	optionally also setting a third conductance, G_3
# 	"""
# 	fig,ax = plt.subplots()

# 	plt.plot(par.biasenergies,G_1[index_of_varconst,:])
# 	plt.hold("on")
# 	plt.plot(par.biasenergies,G_2[index_of_varconst,:])
# 	if len(kwargs) == 1:
# 		plt.hold("on")
# 		plt.plot(par.biasenergies,kwargs["G_3"][index_of_varconst,:])
# 	if len(kwargs) == 2:
# 		plt.hold("on")
# 		plt.plot(par.biasenergies,kwargs["G_4"])
# 	plt.hold("off")

# 	ax.set_title(figtitle,fontsize=16)
# 	ax.set_xlabel(xlabel,fontsize=16)
# 	ax.set_ylabel(ylabel,fontsize=16)
# 	print(len(legend))
# 	ax.legend(legend,loc=5)
# 	plt.grid("on")

# 	plt.savefig(filename + ".pdf")	

# 	print("%s: Saved figure"%misc.round_time(datetime.datetime.now(),round_to=60), filename + ".pdf")

def plot_G_1_G_2_varconst(par,G_1,G_2,index_of_varconst,filename,figtitle,xlabel,ylabel,legend,ylim,**kwargs):
	"""

	Parameters
	----------
	**kwargs : 	optionally also setting a third conductance, G_3
	"""
	fig,ax = plt.subplots()

	ax.plot(par.biasenergies,G_1[index_of_varconst,:])
	ax.hold("on")
	ax.plot(par.biasenergies,G_2[index_of_varconst,:])

	if len(kwargs) == 1:
		print("len kwargs = 1")
		plt.hold("on")
		ax.plot(par.biasenergies,kwargs["G_3"][index_of_varconst,:])
	if len(kwargs) == 2:
		print("len kwargs = 2")
		plt.hold("on")
		ax.plot(par.biasenergies,kwargs["G_3"][index_of_varconst,:])
		ax.hold("on")
		ax.plot(par.biasenergies,kwargs["G_4"][index_of_varconst,:])
	if len(kwargs) == 3:
		print("len kwargs = 3")
		plt.hold("on")
		ax.plot(par.biasenergies,kwargs["G_3"][index_of_varconst,:])
		ax.hold("on")
		ax.plot(par.biasenergies,kwargs["G_4"][index_of_varconst,:])
		ax.hold("on")
		ax.plot(par.biasenergies,kwargs["G_5"][index_of_varconst,:])

	plt.hold("off")

	ax.set_title(figtitle,fontsize=16)
	ax.set_xlabel(xlabel,fontsize=16)
	ax.set_ylabel(ylabel,fontsize=16)
	ax.set_ylim(ylim)
	ax.ticklabel_format(style='sci',scilimits=(-3,4))

	ax.legend(legend,loc=4,fontsize=16)
	plt.grid("on")

	plt.savefig(filename + ".pdf")	

	print("%s: Saved figure"%misc.round_time(datetime.datetime.now(),round_to=60), filename + ".pdf")




def pplot_E_vs_var(par,energies,var_values,figfilename,xlim,ylim,title,xlabel,ylabel,add_zoom,ylim_zoom,yticks,file_extension,u2mv2,u2_mu,v2_mu):
	"""
	Plots energy spectrum versus a variable whose values are given by energies and var_values, respectively.

	Parameters
	----------
	energies : 	array. Shape: (len(var_values),par.k)
	u2mv2 :  	array.
				If you do not want to plot u2mv2, provide is as an empty list [].

	"""
	if add_zoom==False:
		fig,ax = plt.subplots()
	else:
		fig,ax = plt.subplots(nrows=2,sharex='all',gridspec_kw = {'height_ratios':[3, 1]})

	fig.set_size_inches(4.5, 4)
	ax[0].set_title(title)#,fontsize=18)
	ax[0].set_ylabel(ylabel)#,fontsize=16)
	ax[0].set_xlim([xlim[0],xlim[1]])
	ax[0].set_ylim([ylim[0],ylim[1]])
	colors = [plt.cm.PuBu(x) for x in np.linspace(0,1,100)]
	ax[0].plot(var_values,energies,color=colors[-15],linewidth=1.5)
	ax[0].ticklabel_format(style='sci',scilimits=(-3,4))

	
	if add_zoom != False:
		ax[1].set_xlabel(xlabel)#,fontsize=16)
		ax[1].set_xlim([xlim[0],xlim[1]])
		ax[1].set_ylim([ylim_zoom[0],ylim_zoom[1]])
		ax[1].set_yticks(yticks)
		ax[1].ticklabel_format(style='sci',scilimits=(-3,4))
		ax[1].grid("on")
		if len(u2mv2) < 1:
			ax[1].plot(var_values,energies,color=colors[-15],linewidth=1.5)
		# elif u2mv2 is given as input:
		elif len(u2mv2) > 1:
			ax[1].plot(var_values,energies,color=colors[-15],linewidth=1.5)
			ax[1].plot(var_values,u2mv2*0.01,"black")
			ax[1].plot(var_values,u2_mu*0.01,"--")
			ax[1].plot(var_values,v2_mu*0.01,"--")

			# ax[1].scatter(var_values,u2mv2*0.0004,cmap="seismic")
		else:
			Inputerror(" - u2mv2 length cannot be accessed in order to determine whether it is to be plotted.")
	plt.tight_layout()
	ax[0].grid("on")

	
	fig.savefig(figfilename + file_extension,bbox_inches='tight',pad_inches = 0.1)

def pplot_E_u2mv2_1_m1_vs_var(par,E1,Em1,u2mv2_E1,u2mv2_Em1,u2_E1,v2_E1,u2_Em1,v2_Em1,var_values,figfilename,xlim,ylim,title,xlabel,ylabel,file_extension,legend,one_ax,yticks,u2mv2_factor):

	"""
	Plots energy spectrum versus a variable whose values are given by energies and var_values, respectively.

	Parameters
	----------
	E1 : 			array. 
					Shape: (len(var_values),par.k)
	
	u2mv2_E1 :  	array.
					If you do not want to plot u2mv2, provide is as an empty list [].
					(Similarly for the rest of the input arrays one may plot here)
	
	var_values : 	array.
					Values on x-axis.

	xlim : 			array.
					Lowest and highest value on the x-axis, specifying x-axis domain.

	ylim : 			array.
					same, but for y-axis.

	file_extension : 	str.
						e.g. ".pdf", ".eps", ".png", ".jpg"

	legend : 			list.
						Elements are strings for the legend.

	one_ax : 			True/False.
						If True, plot all data in one subplot. If False, plot E_1 and u2mv2_E1 in one subplot and the m1 data in another subplot.

	yticks : 			list.
						Ticks to be displayed on y-axis.

	u2mv2_factor : 		factor to multiply u2mv2 with in order for its values to be comparable to the energies.


	Notes
	-----
	 - until now, u2mv2_factor was larger for longer device (4000 nm) than for shorter device (800 nm).


	Example usage
	-------------
		 - Example: Plotting u^2 - v^2 vs.mu:
		import matplotlib.pyplot as plt
		from My_Costum_Packages_newest import plot

		#		 - For LY = 4000 nm:
		plot.pplot_E_u2mv2_1_m1_vs_var(par,par.E1_mu,par.Em1_mu,u2mv2_E1=par.u2mv2_E1_mu, u2mv2_Em1=par.u2mv2_Em1_mu, u2_E1=par.u2_E1_mu, v2_E1=par.v2_E1_mu, u2_Em1=par.u2_Em1_mu, v2_Em1=par.v2_Em1_mu, var_values=par.mu_values,figfilename="E_u2mv2_1_m1_vsmu_"+par.filename,xlim=[0,1],ylim=[-8e-4,8e-4],title="",xlabel="$\mu\ [meV]$",ylabel="$meV$",file_extension=".pdf",legend=["$E_1$","$(u^2_1-v^2_1)\cdot 0.01$","$E_{-1}$","$(u^2_{-1}-v^2_{-1})\cdot 0.01$"],one_ax=False,yticks=[-8e-4,0,8e-4],u2mv2_factor = 0.01)

		#		 - For LY = 400 nm:
		plot.pplot_E_u2mv2_1_m1_vs_var(par,par.E1_mu,par.Em1_mu,u2mv2_E1=par.u2mv2_E1_mu, u2mv2_Em1=par.u2mv2_Em1_mu, u2_E1=par.u2_E1_mu, v2_E1=par.v2_E1_mu, u2_Em1=par.u2_Em1_mu, v2_Em1=par.v2_Em1_mu, var_values=par.mu_values,figfilename="E_u2mv2_1_m1_vsmu_"+par.filename,xlim=[0,1],ylim=[-0.3,0.3],title="",xlabel="$\mu\ [meV]$",ylabel="$meV$",file_extension=".pdf",legend=["$E_1$","$(u^2_1-v^2_1)\cdot 0.1$","$E_{-1}$","$(u^2_{-1}-v^2_{-1})\cdot 0.1$"],one_ax=False,yticks=[-0.3,0,0.3],u2mv2_factor = 0.1)
		

		 - Example: Plotting u^2 - v^2 vs.Ez:
		import matplotlib.pyplot as plt
		from My_Costum_Packages_newest import plot

		plot.pplot_E_u2mv2_1_m1_vs_var(par,par.E1_Ez,par.Em1_Ez,u2mv2_E1=par.u2mv2_E1_Ez, u2mv2_Em1=par.u2mv2_Em1_Ez, u2_E1=par.u2_E1_Ez, v2_E1=par.v2_E1_Ez, u2_Em1=par.u2_Em1_Ez, v2_Em1=par.v2_Em1_Ez, var_values=par.Ez_values,figfilename="E_u2mv2_1_m1_vsEz_"+par.filename,xlim=[par.Ez_values[0],par.Ez_values[-1]],ylim=[-0.4,0.4],title="",xlabel="$Ez\, [meV]$",ylabel="$meV$",file_extension=".pdf",legend=["$E_1$","$(u^2_1-v^2_1)\cdot 0.1$","$E_{-1}$","$(u^2_{-1}-v^2_{-1})\cdot 0.1$"],one_ax=False,yticks=[-0.4,0,0.4],u2mv2_factor = 0.1)


	"""



	# PuBu_colors = [plt.cm.PuBu(x) for x in np.linspace(-1,1,100)]
	# my_blue = PuBu_colors[-15]
	my_blue = plt.cm.get_cmap("PuBu",lut=100)(85)
	# my_blue = my_blue(0)
	my_red = plt.cm.get_cmap("OrRd")(200)	# Dark red around argument=500 and up
	my_orange = plt.cm.get_cmap("autumn")(0.4)
	# my_dark_red = plt.cm.get_cmap("OrRd")(1000)


	if one_ax == True:
		fig,ax = plt.subplots()#gridspec_kw = {'height_ratios':[3, 1]}

		ax.plot(var_values,E1,color=my_blue,linewidth=1.5)
		ax.plot(var_values,Em1,color=my_red,linewidth=1.5)
		ax.plot(var_values,u2mv2_E1*u2mv2_factor,color=black,linewidth=1.5)
		ax.plot(var_values,u2mv2_Em1*u2mv2_factor,color=black,linewidth=1.5)

		ax.ticklabel_format(style='sci',scilimits=(-3,4))

		ax.set_title(title)#,fontsize=18)
		ax.set_ylabel(ylabel)#,fontsize=16)
		ax.set_xlabel(xlabel)
		ax.set_xlim([xlim[0],xlim[1]])
		ax.set_ylim([ylim[0],ylim[1]])

	elif one_ax == False:
		fig,ax = plt.subplots(nrows=2,sharex='all',sharey='all')

		ax[0].plot(var_values,E1,color=my_blue,linewidth=1.5)
		ax[1].plot(var_values,Em1,color=my_red,linewidth=1.5)
		ax[0].plot(var_values,u2mv2_E1*u2mv2_factor,color="black",linewidth=1.5)
		ax[1].plot(var_values,u2mv2_Em1*u2mv2_factor,color="black",linewidth=1.5)

		ax[0].ticklabel_format(style='sci',scilimits=(-3,4))
		ax[1].ticklabel_format(style='sci',scilimits=(-3,4))



		ax[0].set_title(title)#,fontsize=18)

		ax[0].yaxis.set_label_coords(-0.12,-0.23)

		ax[0].set_ylabel(ylabel)#,fontsize=16)
		ax[1].set_xlabel(xlabel)
		ax[0].set_xlim([xlim[0],xlim[1]])
		ax[0].set_ylim([ylim[0],ylim[1]])
		ax[1].set_xlim([xlim[0],xlim[1]])
		ax[1].set_ylim([ylim[0],ylim[1]])

		ax[0].set_yticks(yticks)
		ax[1].set_yticks(yticks)

		ax[0].grid(which="both")
		ax[1].grid(which="both")	

		ax[0].legend(legend[0:2],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=16)#(handles=[line1], loc=1)
		ax[1].legend(legend[2:4],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=16)

	fig.set_size_inches(4.5, 3)
	plt.subplots_adjust(hspace=0.4)

	# plt.tight_layout()

	
	fig.savefig(figfilename + file_extension,bbox_inches='tight',pad_inches = 0.1)




def plot_u2mv2(par):
	plt.hold("off")
	plt.plot(par.mu_values,np.real(par.u2mv2_mu))
	plt.hold("on")
	plt.plot(par.mu_values,np.imag(par.u2mv2_mu))
	# plt.hold("on")
	# plt.plot(par.mu_values[peaks],np.real(par.u2mv2_mu[peaks]),"x")

	plt.xlabel("$\mu\ [meV]$")
	plt.ylabel("$u^2 - v^2$")
	plt.title("$E_z=1.6\,meV$")
	plt.savefig("u2mv2_vs_mu_" + par.filename + ".pdf")
	# plt.hold("on")
	# plt.plot(par.mu_values,np.imag(par.u2mv2_mu))
	print(3)
	plt.show()



def plot_nu_nuorb_vs_mu(nu, mu, filename):
	"""
	Notes
	-----
	 - assumes each orbital has same nu.
	"""

	fig,ax = plt.subplots()
	fig.set_size_inches(4,3)
	for i in range(len(nu[0,:])):
		plt.plot(mu,nu[:,i],"--",linewidth=1.5)
	# ??? Which of the four orbitals is nu of? sum of all?
	nu = np.sum(nu,axis=1)
	plt.grid("on")
	plt.plot(mu,nu,linewidth=1.5)
	# plt.plot(omega,nu)
	plt.xlabel("$\mu\ [meV]$")
	plt.ylabel(r"$\nu\ [meV^{-1} nm^{-3}]$")
	# plt.legend(["Single orbital","Orbitals summed"],fontsize=12)
	plt.tight_layout()
	plt.savefig("nu_vs_mu_"+filename+".pdf")