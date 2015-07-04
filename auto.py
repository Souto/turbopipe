import numpy as np
import matplotlib.pyplot as plt
import os
import pyconvert
import glob as glob

def run():	


	os.chdir('/home/diogo/turbopipe/turbospec_script/')					#set the path of system
	os.system('csh -f auto_dwarfs_O.com')							#run turbospectrum script
	
	os.chdir('/home/diogo/turbopipe/trash/')						#set the path of system
	Oxygen_syn = glob.glob('teste_trash_Oxygen_*.txt')
	for a in Oxygen_syn:
		pyconvert.faltbo(a,"convol"+a,710.,2);print 'Making a gaussian profile to ', a	#runing a gaussian smooth
		
			
	
	print 'Program finished well'
	return
	
if __name__=='__main__':
	run()
