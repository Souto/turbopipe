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
	
	'''
	# Seleting wavelenths for chi squared adjust
	aux = np.where((lamb>lim1) & (lamb<lim2))[0]
	# Adjusting 
	O_chi = chi(flux1[aux],flux2[1][aux],ddof=9,axis=0)
	
	# plot chi2 method
	#if '-plot' in sys.argv:
	#	plt.scatter(min(aux)-20,min(aux)+20,color='c')
	#	plt.scatter(lim1,syn_O[1],color='r')
	#	plt.show()

	# Return the chi2 coeficient
	return O_chi

###################################
#          Run script 	          #
###################################
def run():	
	'''run module
	'''
	# Changing path do open files
	os.chdir('/home/diogo/turbopipe/trash/')

	# Open observed spectra to be used in chi2 method
	kp138 = np.loadtxt('starkp138_newtest.txt',usecols=[0,1],skiprows=1).T

	# Defining common frame and interpolating stuff
	xvals_kp138 = np.linspace(15000,17000,20001)
	yinterp_kp138 = np.interp(xvals_kp138,kp138[0],kp138[1])
	
		
	# Searching for synthetic spectra
	syn_Oxygen = np.sort(glob.glob('convol_teste_trash_Oxygen_*.txt'))

	# Open all synthetic spectra to Oxygen analysis
	syn_O = np.array([(np.loadtxt(i,usecols=[0,1]).T) for i in syn_Oxygen])


	# Open the file with the wavelenght limits to be used in chi2 method
	O_lim_inf,O_lim_sup = np.loadtxt('Oxygen_lim.dat',usecols=[0,1],dtype='float').T

	
	# Run the function chi2fun and save each value of chi2 for each region analysed
	O_chi = np.array([np.array([chi2fun(yinterp_kp138,syn_O[k],xvals_kp138,i,j) for i,j in zip(O_lim_inf, O_lim_sup)]).T[0] for k in xrange(len(syn_O))]).T


	# Find the minimum chi2 for the several abundances tested
	best_O_chi = ([min(O_chi[i]) for i in xrange(len(O_lim_inf))])


	# Find the respective abundance used in the best fit (minimum chi2)
	ab_O_chi = np.array([np.where(O_chi[i]==min(O_chi[i])) for i in xrange(len(O_lim_inf))])

	# Save all the abundances used to determine best chi2
	Oxygen_values = ([str(i).split('_')[-1].split('.t')[0] for i in syn_Oxygen])	
	
	
	# Save the abundances value correspondent to lower chi2 for all 
	# regions analysed in 'Oxygen_lim.dat'
	Oxygen_abundances = ([str(i).split('_')[-1].split('.t')[0] for i in syn_Oxygen[ab_O_chi]])
	Oxygen_abundances = np.float64(Oxygen_abundances)
	Oxygen = [np.mean(Oxygen_abundances),np.std(Oxygen_abundances)]

	###	ploting figures	###
	
	print "Making and Saving figures"
	for i in xrange(len(O_lim_inf)):
		plt.figure(i)
		plt.title('Star Name',fontsize=8)
		plt.plot((syn_O[ab_O_chi])[i][0][0][0],(syn_O[ab_O_chi])[i][0][0][1],'-b',alpha=0.8)
		plt.plot(syn_O[0][0],syn_O[0][1],'--',color='black',alpha=0.7)
		plt.plot(syn_O[len(syn_O) -1][0],syn_O[len(syn_O) -1][1],'--r',alpha=0.7)
		plt.scatter(kp138[0],kp138[1],color='g',marker='o',s=15)
		plt.legend(('Best fit with A(O) = ' +str(Oxygen_abundances[i]) +' $\pm$ '+str(round(best_O_chi[i], 3)) +' dex',"Minimun of A(O) = " +str(min(Oxygen_values)) +" dex analysed","Max of A(O) = " +str(max(Oxygen_values)) +" dex analysed","Observed Spectra of NameStar"),fontsize=10, loc='lower left', scatterpoints=1, markerscale=1,   frameon=False)		
		plt.xlim(O_lim_inf[i]-5,O_lim_sup[i]+5)
		plt.ylim(0.5,1.1)
		plt.ticklabel_format(useOffset=False)	
		plt.xlabel(r'Wavelength ($\AA{}$)',fontsize=12);plt.xticks(fontsize=10)
		plt.ylabel(r'Intensity',fontsize=12);plt.yticks(fontsize=10)
#		plt.text(O_lim_inf[i]-4.8,0.65,r'Best Abundance A(O) = '+str(Oxygen_abundances[i]) +' $\pm$ '+str(best_O_chi[i]) +' dex',fontsize=10,style='italic')
#		plt.show()
		plt.savefig('figure'+str(i) +'.jpg',dpi=300)
		plt.clf
	
	
	print "Chi2.py finished well"
	
	
if __name__=='__main__':
	run()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
