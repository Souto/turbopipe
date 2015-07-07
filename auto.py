#!

'''auto.py: main script to obtain automatic abundances for cool dwarfs
   	    using the turbospectrum code from B. Plez.
   	    
   	    This code consist in obtain the stellar abundance using the chi-square
   	    method to find a best fit between the synthetic and observed spectra.

	'''

__author__ = 'Diogo Souto'
__email__  = 'souto.at.on.br'

 
###################################
#       Importing modules         #
###################################

import numpy as np
import matplotlib.pyplot as plt
import os
import pyconvert
import glob as glob

###################################
#       Auxiliary Methods         #
###################################


###################################
#          Run script 	          #
###################################


def run():	

	# Set the path of the system 
	os.chdir('/home/diogo/turbopipe/turbospec_script/')					

	# Run turbospectrum script - The first element should be Oxygen, because there is several lines of.
	os.system('csh -f auto_dwarfs_O.com')							
	
	# Set the path of the system
	os.chdir('/home/diogo/turbopipe/trash/')
	
	# Searching for synthetic spectra 		
	Oxygen_syn = glob.glob('teste_trash_Oxygen_*.txt')
	
	# The pyconvert.faltbo it's a python version of the faltbo.f90 file from turbospectrum, B. Plez
	for a in Oxygen_syn:
		pyconvert.faltbo(a,"convol_"+a,710.,2);print 'Making a gaussian profile to ', a	
		

	# Open main folder and run chi2 method
	os.chdir('/home/diogo/turbopipe/')	
	os.system('python chi2.py')


	print 'Program finished well'
	return
	
if __name__=='__main__':
	run()
