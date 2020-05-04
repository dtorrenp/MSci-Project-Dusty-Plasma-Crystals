import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as stats
from csv import reader
import subprocess
import time
import pandas as pd
import matplotlib.animation as animation
from matplotlib.patches import Circle,Rectangle
from scipy.spatial import Voronoi,voronoi_plot_2d,Delaunay
import os
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,mark_inset
from scipy.fftpack import fft,fftfreq
import matplotlib.mlab as mlab
#%%
#DEFINE PLASMA DISCHARGE CONDITIONS

n_e0 = 1e15#electron and ion densities in bulk plasma
n_i0 = 1e15
e_charge = -1.6*1e-19
grain_R = 7*1e-6
m_i = 1.67*1e-27
m_e = 9.11*1e-31
m_D = ((4/3)*np.pi*grain_R**3)*(1.49*1e3)#mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
epsilon_0 = 8.85*1e-12
k_b = 1.38*1e-23
T_e = (2*1.6*1e-19)/k_b
T_i = (0.03*1.6*1e-19)/k_b
T_n_a = (0.03*1.6*1e-19)/k_b
lambda_de = ((epsilon_0*k_b*T_e)/(n_e0*(e_charge**2)))**0.5
lambda_di = ((epsilon_0*k_b*T_i)/(n_i0*(e_charge**2)))**0.5
lambda_D = (1/(1/(lambda_de**2) + 1/(lambda_di**2)))**0.5#get lambda_D
container_radius = 40*lambda_D#set radius of contianer ie wall radius
container_length = 20*lambda_D
z_se = 29.1*lambda_D#distance from bottom of container to the sheath edge
top_container_graph_mul = 1.5

final_pos_plot = "yes"


#%%
DATA_NAME_1 = "Finite_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_386.499819_frames_500001"#input("Data name_1?[PERIODIC](NOT FILENAME)")
DATA_NAME_2 = "Finite_Dust_grain_max_100_Tn_EV_0.030000_Final_Temperature_400.139865_frames_100002"#input("Data name_2?[FINITE](NOT FILENAME)")

FILENAME_1 = "HPC_Data/"+ DATA_NAME_1 + ".csv"
FILENAME_2 = "HPC_Data/"+ DATA_NAME_2 + ".csv"

NEW_FOLDER = "Figures/" + "Bleendin_Compare_two_dust_grains"# + DATA_NAME_1 + "_and_" + DATA_NAME_2
if str(os.path.exists(NEW_FOLDER)) == "False":
    #os.mkdir("NEW_FOLDER")
    os.makedirs(NEW_FOLDER)

#%%
data_1 = pd.read_csv(FILENAME_1)
data_2 = pd.read_csv(FILENAME_2)

#data_list = [data_1,data_2]

dust_grain_max_1 = int( (len(data_1.columns))/3 - 1)
dust_grain_max_2 = int( (len(data_2.columns))/3 - 1)

#dust_grain_max_list = [dust_grain_max_1,dust_grain_max_2]
frames_len = len(data_1["Time_list"])
y_z_se = [z_se/lambda_D]*len(data_1["Time_list"])
last_time_val = data_1["Time_list"].iloc[-1]
#%%
if final_pos_plot == "yes":

    fig, axs = plt.subplots(1,2)
    axs[0][0].grid()
    axs[0][0].text(0.01, 0.95, "a)", transform=axs[0].transAxes, fontsize=12, va='top')
    axs[0][0].set_xlabel(r"$x/\lambda_D$")
    axs[0][0].set_ylabel(r"$y/\lambda_D$")
    for i in np.arange(dust_grain_max_1):
        last_val_index = np.where(data_1["Time_list"] == last_time_val)
        axs[0][0].plot(data_1["X_" + str(i)].iloc[last_val_index], data_1["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[0][1].set_xlim(-container_radius/lambda_D,container_radius/lambda_D)
    axs[0][1].set_ylim(-container_radius/lambda_D,container_radius/lambda_D)
    
    axs[0][1].grid()
    #axs[1].set_xlabel(r"$z/\lambda_D$")
    axs[0][1].text(0.01, 0.95, "b)", transform=axs[1].transAxes, fontsize=12,  va='top')
    axs[0][1].set_xlabel(r"$x/\lambda_D$")
    axs[0][1].set_ylabel(r"$y/\lambda_D$")
    for i in np.arange(dust_grain_max_2):
        last_val_index = np.where(data_1["Time_list"] == last_time_val)
        axs[0][1].plot(data_2["X_" + str(i)].iloc[last_val_index], data_2["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[0][1].set_xlim(-container_radius/lambda_D,container_radius/lambda_D)
    axs[0][1].set_ylim(-container_radius/lambda_D,container_radius/lambda_D)
    axs[0][1].plot(z_list,charge_list)
    axs[0][1].axvline(29.1,color='black', linestyle='dashed')

    plt.savefig("Figures/brah.png", dpi = 1200)
