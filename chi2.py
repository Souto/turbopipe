#!

'''chi2.py: To obtain the best abundance fit

	'''

__author__ = 'Diogo Souto'
__email__  = 'souto.at.on.br'

 
###################################
#       Importing modules         #
###################################

import numpy as np
import matplotlib.pyplot as plt
import os
import glob as glob
from scipy.stats import chisquare as chi
import sys

###################################
#       Auxiliary Methods         #
###################################


def lim(lamb,lim1,lim2):

	# The function lim() return each index for the file *lim.txt
	return  np.where((lamb>lim1) & (lamb<lim2))[0]
	''' Method for selecting and fiting region with chi-squared
	
	Parameters
	----------
	lamb: numpy.ndarray
		array with wavelenghs to be used in chi2 and plots
	
	lim1: float
		inferior limit for chi2 method 
		
	lim2: float
		superior limit for chi2 method
		
	Returns
	-------
	An array with the index for each region selected
	
	'''



def chi2fun(flux1,flux2,lamb,lim1,lim2):
	''' Method for selecting and fiting region with chi-squared
	
	Parameters
	----------
	flux1: numpy.ndarray
		array with fluxes for interpoleted observed spectra 
		
	flux2: numpy.ndarray
		array for synthetic spectra 
		
	lamb: numpy.ndarray
		array with wavelenghs to run chi2 method 
	
	lim1: float
		inferior limit for chi2 method 
		
	lim2: float
		superior limit for chi2 method f
		
	Returns
	-------
	O_chi : numpy.ndarray
		value of chi2 and p value
	
	
	Comments: It's possible look at the each region by each chi2 analyse in each abundance ran
		  using the command at the terminal: python chi2.py plot
		  The word plot will allow to show this.
	'''


	# Seleting wavelenths for chi squared adjust
	aux = lim(lamb,lim1,lim2)
	# Adjusting and renormalyzing the spectra
	O_chi = chi(norm(flux1[aux]),norm(flux2[1][aux]),ddof=9,axis=0)

	# Plot each region for Abundance
	if 'plot' in sys.argv:
		plt.scatter(aux,norm(flux1[aux]),marker='o')
		plt.scatter(aux,norm(flux2[1][aux]),marker='^',color='g')
		plt.legend(("Observed","synthetic"),loc="lower left")
		plt.plot([min(aux),max(aux)],[1,1],'--k')		
		plt.plot([min(aux),max(aux)],[1,1],'--k')
		plt.show()

	# Return the chi2 coeficient and the limits of each region
	
	return O_chi[0]

##################################################################################################


def norm(norm1):
	'''Method for renormalize the synthetic and/or the observed spectra in a
	   region to be used in the chi2 analyses.
	   
	   A simple renormalization consist in take the maximum value of the array and divide
	   all array by it. However, we assume here that 15 % of the higher value of the flux 
	   (above the continuum) should be an error and we neglect these data.
	   
	   In other words, the spectra were renormalized using the higher value at 85 % of all data
	   
	   
	Parameters
	----------
	flux1: numpy.ndarray
		Array with a vector to be normalized
		
		
	Returns
	-------
	norm : numpy.ndarray
		Array with renormalized spectra
	
	'''
	norm = norm1/(np.sort(norm1)[len(norm1) - int(0.15*len(norm1))])

	return norm




###################################
#         Ploting Figures         #
###################################


	
	
###################################
#          Run script 	          #
###################################
def run():	
	'''Run the main module for chi-square analyses.
	   That routine takes all synthetic spectra that were made in
	   auto.py and uses chi2 method to find, store and plot the best
	   fit for each abundance in each region analysed in the spectra.
	   
	'''
	# Changing path do open files
	os.chdir('/home/diogo/turbopipe/trash/')
	
	# Open observed spectra to be used in chi2 method
	kp138 = np.loadtxt('starkp138_newtest.txt',usecols=[0,1],skiprows=1).T

	# Defining common frame and interpolating stuff
	xvals_kp138 = np.linspace(15000,17000,40001)
	yinterp_kp138 = np.interp(xvals_kp138,kp138[0],kp138[1])
	
		
	# Searching for synthetic spectra
	syn_Oxygen = np.sort(glob.glob('convol_teste_trash_Oxygen_*.txt'))

	# Open all synthetic spectra to Oxygen analysis
	syn_O = np.array([(np.loadtxt(i,usecols=[0,1]).T) for i in syn_Oxygen])


	# Open the file with the wavelenght limits to be used in chi2 method
	O_lim_inf,O_lim_sup = np.loadtxt('Oxygen_lim.dat',usecols=[0,1],dtype='float').T

	
	# Run the function chi2fun and store each value of chi2 for each region analysed
	O_chi = np.array([[chi2fun(yinterp_kp138,syn_O[k],xvals_kp138,i,j) for i,j in zip(O_lim_inf, O_lim_sup)] for k in xrange(len(syn_O))]).T

	# Find the minimum chi2 for the several abundances tested
	best_O_chi = ([min(O_chi[i]) for i in xrange(len(O_lim_inf))])


	# Find the respective abundance used in the best fit (minimum chi2)
	ab_O_chi = np.ravel([np.where(O_chi[i]==min(O_chi[i]))[0] for i in xrange(len(O_lim_inf))])

	# Save all the abundances used to determine best chi2
	Oxygen_values = ([str(i).split('_')[-1].split('.t')[0] for i in syn_Oxygen])	
	Oxygen_values = np.float64(Oxygen_values)
	
	# Save the abundances value correspondent to lower chi2 for all 
	# regions analysed in 'Oxygen_lim.dat'
	Oxygen_abundances = ([str(i).split('_')[-1].split('.t')[0] for i in syn_Oxygen[ab_O_chi]])
	Oxygen_abundances = np.float64(Oxygen_abundances)
	Oxygen = [np.mean(Oxygen_abundances),np.std(Oxygen_abundances)]


	print "Making and Saving figures"

	'''
	Defining the subplot type. I choose to plot the figure in 3 window
	1 - The observed spectra overploted with best abundances and min and max abundance analysed
	2 - The residual difference between the observed and best synthesis
	3 - Abundance over chi2 that shows the minimum chi2 obtained corralated with the best abundance 
	    the sigma of the chi2 are shown as dashed line.

	Returns
	-------
	All figures related to Oxygen abundance
	'''
	
	# Makes one plot for region
	for i in xrange(len(O_lim_inf)):

		# Defining the format of the figure
		plt.figure(i)
		plt.subplots(nrows=2, ncols=2,figsize=(15,10))
		plt.subplots_adjust(left  = 0.07,right = 0.95,bottom = 0.07,top = 0.95,wspace = 0.15,hspace = 0.15)

		# plot window 1 (i,j) = (0,[0-1])
		plt.subplot2grid((2, 2), (0, 0),colspan=2)

		# Stellar ID on the title of the figure
		plt.title('Star Name',fontsize=8)
	
		# Selecting region of lim function for the synthetic and observed spectra
		region_aux = lim(xvals_kp138,O_lim_inf[i]-6,O_lim_sup[i]+6)
		region_aux_obs = lim(kp138[0],O_lim_inf[i]-6,O_lim_sup[i]+6)		

		# Plot the synthetis. Best fit (blue line); Max abundance analysed (red dashed); Min abundance analysed (black dashed)
		plt.plot(syn_O[ab_O_chi][i][0][region_aux],norm(syn_O[ab_O_chi][i][1][region_aux]),'-b',alpha=0.8)
		plt.plot(syn_O[0][0][region_aux],norm(syn_O[0][1][region_aux]),'--',color='black',alpha=0.7)
		plt.plot(syn_O[len(syn_O) -1][0][region_aux],norm(syn_O[len(syn_O) -1][1][region_aux]),'--r',alpha=0.7)
		
		# Plot the observed spectra in green dots
		plt.scatter(kp138[0][region_aux_obs],norm(kp138[1][region_aux_obs]),color='g',marker='o',s=15)
	
		# Define box stuffs
		plt.legend(('Best fit with A(O) =  %.2f $\pm$ %.4f dex' % (Oxygen_abundances[i], best_O_chi[i]),"Minimun of A(O) = " +str(min(Oxygen_values)) +" dex analysed","Maximum of A(O) = " +str(max(Oxygen_values)) +" dex analysed","Observed Spectra of NameStar"),fontsize=10, loc='lower left', scatterpoints=1, markerscale=1,   frameon=False)		
		plt.xlim(O_lim_inf[i]-5,O_lim_sup[i]+5)
		plt.ylim(0.5,1.1)
		plt.ticklabel_format(useOffset=False)	
		plt.xlabel(r'Wavelength ($\AA{}$)',fontsize=12);plt.xticks(fontsize=8)
		plt.ylabel(r'Intensity',fontsize=12);plt.yticks(fontsize=8)


		# plot window 2 (i,j) = (1,0)
		plt.subplot2grid((2, 2), (1,0))
		
		# Plot the residual - difference between best fit and observed spectra -.
		plt.plot(syn_O[ab_O_chi][i][0][region_aux],norm(syn_O[ab_O_chi][i][1][region_aux]) - norm(yinterp_kp138[region_aux]),'-b',alpha=0.8)
		
		# Plot a dashed line indicating the residual of 0.
		plt.plot([O_lim_inf[i]-5,O_lim_sup[i]+5],[0,0],'--k',linewidth=1.5)

		# Define box stuffs
		plt.xlim(O_lim_inf[i]-5,O_lim_sup[i]+5)
		plt.ylim(-0.2,0.2)
		plt.ticklabel_format(useOffset=False)	
		plt.xlabel(r'Wavelength ($\AA{}$)',fontsize=12);plt.xticks(fontsize=8)
		plt.ylabel(r'Residual',fontsize=12);plt.yticks(fontsize=8)

	
		# plot window 3 (i,j) = (1,1)
		plt.subplot2grid((2, 2), (1,1))
		
		# Plot the abundances of element over chi2 results
		plt.plot(Oxygen_values,O_chi[i],alpha=0.8)
		
		# Makes a line to shown the minimun chi2 obtained and the dashed lines are the standart deviation of chi2 method
		plt.plot([Oxygen_abundances[i],Oxygen_abundances[i]],[best_O_chi[i]+np.std(O_chi[i]),max(O_chi[i])-np.std(O_chi[i])],'-k')
		plt.plot([Oxygen_abundances[i]-np.std(O_chi[i]),Oxygen_abundances[i]-np.std(O_chi[i])],[best_O_chi[i]+np.std(O_chi[i]),max(O_chi[i])-np.std(O_chi[i])],'--k')
		plt.plot([Oxygen_abundances[i]+np.std(O_chi[i]),Oxygen_abundances[i]+np.std(O_chi[i])],[best_O_chi[i]+np.std(O_chi[i]),max(O_chi[i])-np.std(O_chi[i])],'--k')	
		# Define box stuffs
		plt.xlim(min(Oxygen_values)-.03,max(Oxygen_values)+.03)
		plt.ticklabel_format(useOffset=False)
		plt.ylabel(r'$\chi^{2}$ values',fontsize=12);plt.xticks(fontsize=8)
		plt.xlabel(r'Oxygen abundance dex',fontsize=12);plt.yticks(fontsize=8)
#		plt.show()
		plt.savefig('new_figure'+str(i) +'.jpg',dpi=200)
		plt.clf()
	










	# End of the program
	print "Chi2.py finished well"
	
	
if __name__=='__main__':
	run()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
