import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
from csv import reader
import subprocess
import time
import pandas as pd
import matplotlib.animation as animation
import scipy.stats as stats
import time

#%%
#DEFINE PLASMA DISCHARGE CONDITIONS

n_e0 = 1e15#electron and ion densities in bulk plasma
n_i0 = 1e15
e_charge = -1.6*1e-19
i_charge = 1.6*1e-19
grain_R = 7*1e-6
m_i = 1.67*1e-27
m_e = 9.11*1e-31
m_D = ((4/3)*np.pi*grain_R**3)*(1.49*1e3)#mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
epsilon_0 = 8.85*1e-12
k_b = 1.38*1e-23
T_e = (2*1.6*1e-19)/k_b
T_i = (0.03*1.6*1e-19)/k_b
lambda_de = ((epsilon_0*k_b*T_e)/(n_e0*(e_charge**2)))**0.5
lambda_di = ((epsilon_0*k_b*T_i)/(n_i0*(e_charge**2)))**0.5
lambda_D = (1/(1/(lambda_de**2) + 1/(lambda_di**2)))**0.5#get lambda_D
container_height = 11*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 25*lambda_D#set radius of contianer ie wall radius
container_length = 10*lambda_D
z_se = 10*lambda_D#distance from bottom of container to the sheath edge
r_se = 25*lambda_D#distance from wall to the sheathe edge
r_se_inv = container_radius - r_se
top_container_graph_mul = 1.5
hist_div_val = 4 #larger less bins, smaller more bins
colour_list = ["red","green","blue","orange","black","brown"]
#%%

status = input("Compile or Run?")
boundry = input("Periodic or Finite?")
layer_plots = input("layers?")

if status == "Compile":
    #COMPILE
    #UNI
    #subprocess.call(["g++", "H:\year 4\computational\MSci-Project-Dusty-Plasma-Crystals\MSci_project.cpp"])
    #LAPTOP
    #os.environ["PROJECT_FILE"]
   if boundry == "Periodic":
      subprocess.call(["g++", "-o", "MSci_project_periodic_HPC", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/MSci_project_periodic_HPC.cpp"])
      """Import stuff to measure how long the code takes to run"""
      start_time = time.time()
      print("start_time =", time.ctime(time.time()))
      subprocess.call("MSci_project_periodic_HPC.exe")
      """prints time taken in minutes"""
      print ("time taken: %s minutes" % ((time.time()-start_time)/60))
   elif boundry == "Finite":
      subprocess.call(["g++", "-o", "MSci_project_HPC", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/MSci_project_HPC.cpp"])
      """Import stuff to measure how long the code takes to run"""
      start_time = time.time()
      print("start_time =", time.ctime(time.time()))
      subprocess.call("MSci_project_HPC.exe")
      """prints time taken in minutes"""
      print ("time taken: %s minutes" % ((time.time()-start_time)/60))
   else:
      print ("run again")


#%%

FILENAME = input("Data file name?")

#%%
data = pd.read_csv(FILENAME)
dust_grain_max = int((len(data.columns))/4)
last_time_val = data["Time_list_0"].iloc[-1]
y_z_se = [z_se/lambda_D]*len(data["Time_list_0"])
temp_ion = [T_i]*len(data["Time_list_0"])

"""
Animation of simulation is produced. Patches are moved on the axes corresponding to the movement of particles
"""

fig = plt.figure(1)
ax = plt.axes(xlim=(container_radius, -container_radius), ylim=(container_radius,-container_radius))
ax.axes.set_aspect('equal')
movie = G.Gas(Animation,ax,frame_rate,balls,frame_limit,container_radius,brownian_mass,Temp,Diatomic)#(self,animate,ax,frame_rate,Balls,frame_limit,container_radius,Brownian_mass,temperature)
anim = animation.FuncAnimation( fig,
                                movie.next_frame,
                                init_func = movie.init_figure, 
                                frames = 1000, 
                                interval = 1,
                                blit = True)