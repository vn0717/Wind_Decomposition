#MPI STARTS HERE
from mpi4py import MPI
#these packages load for each processor
from vort_div import get_wind
import numpy as np



#get MPI information
comm = MPI.COMM_WORLD
#What processor are we on
rank = comm.Get_rank()

#How many processors do we have
size = comm.Get_size()

#0 is the root processor in this case and we will initalize things here
if rank == 0:
	#these packages I want loaded only on the root processor
	import time
	import matplotlib.pyplot as plt
	
	dx = 5000
	dy = 5000
	vort = np.ones((250,250)) * 0.000001
	div = np.ones((250,250)) * 0.000001
	limit = 2000000
	
	#You can only broadcast one variable to processors.  So I put
	#things into a dictonary to make it easier than multiple broadcast
	#calls
	data = {"DX" : dx, "DY" : dy, "vort" : vort, "div" : div, "limit", limit}
	
	#for timing performance
	start_time = time.time()
	
	#for initalizing arrays later when we come back to the 
	#root processor
	size = np.shape(vort)
	

#if we are not on the root processor, the variable we are going to 
#broadcast still needs to exist on each processor.
else:
	data = None

#broadcast data to all processors.  Broadcast is used instead of 
#scatter because we need the whole array for calculations, not
#just part.
data = comm.bcast(data, root=0)
	
#call the subroutine to calcluate the wind decomosition on each processor
out_data = get_wind.vort_div_wind(data["DX"],data["DY"], data["vort"], data["div"], data["limit"])

#take all of the data from the processor and package it up into
#a dictonary to send back to the root processor.  Once again
#you can only use one variable.
send_data = {"U_CHI": get_wind.u_chi, \
"V_CHI": get_wind.v_chi, \
"U_PSI": get_wind.u_psi, \
"V_PSI":get_wind.v_psi, \
"start":out_data[0], \
"end":out_data[1]}

#send to root processor
data = comm.gather(send_data,root=0)

#if this is the root processor gather the data into one array for each variable
#if you don't do this the data is spread across the processors
if rank == 0:
	#initalize final data arrays
	u_psi_final = np.empty(size)
	v_psi_final = np.empty(size)
	u_chi_final = np.empty(size)
	v_chi_final = np.empty(size)
	
	#for the data from each processor, take the variables and 
	#put them into a single array 
	for sub_data in data:
		#one is subtracted since fortran arrays start at 1 rather than 0
		start = sub_data["start"] - 1
		#fortran is inclusive, whereas python isn't, so we don't have to subtract 1
		end = sub_data["end"] 
		u_psi_final[:,start:end] = sub_data["U_PSI"]
		v_psi_final[:,start:end] = sub_data["V_PSI"]
		u_chi_final[:,start:end] = sub_data["U_CHI"]
		v_chi_final[:,start:end] = sub_data["V_CHI"]
	
	#print how long it took
	print(time.time() - start_time)
	#plt.pcolormesh(u_psi_final)
	#plt.show()

#close out the MPI sesson.  This does not stop the need
#to use if statements to select the processor to work with
#If if statements are not used, whatever you are doing will
#be done on all processors.	
MPI.Finalize()
	

