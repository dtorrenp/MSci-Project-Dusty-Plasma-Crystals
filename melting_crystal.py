import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
from csv import reader
import subprocess
import time
import pandas as pd
import matplotlib.animation as animation
from matplotlib.patches import Circle,Rectangle
from scipy.spatial import Voronoi,voronoi_plot_2d,Delaunay
import os

plt.rcParams['animation.ffmpeg_path'] = "C:/Users/daniel/Anaconda3/Library/bin/ffmpeg.exe"

#%%
#DEFINE PLASMA DISCHARGE CONDITIONS
dt = 10**(-4)
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
T_n_b = (2.0*1.6*1e-19)/k_b
lambda_de = ((epsilon_0*k_b*T_e)/(n_e0*(e_charge**2)))**0.5
lambda_di = ((epsilon_0*k_b*T_i)/(n_i0*(e_charge**2)))**0.5
lambda_D = (1/(1/(lambda_de**2) + 1/(lambda_di**2)))**0.5#get lambda_D
container_radius = 40*lambda_D#set radius of contianer ie wall radius
container_length = 20*lambda_D
z_se = 29.1*lambda_D#distance from bottom of container to the sheath edge
top_container_graph_mul = 1.5#larger less bins, smaller more bins
colour_list = ["red","green","blue","orange","black","brown"]
Brown_num = 1000
#%%

status = "Run"#input("Compile or Run?")
load = "no"
boundry = "Periodic"#input("Periodic or Finite?")
standard_plots = "no"
layer_plots = "no"#input("layers?")
melting_plots_g = "yes"
melting_plots_snap = "no"
pos_plot_six = "yes"
temp_log_avoid = "yes"

layers = [0.0,50.0]
temp_list = [2*1e4,4*1e4,6*1e4,105291,105728,106424]
temp_lim_index = 105291

x_plot = "no"
y_plot = "no"
z_plot = "yes"
final_pos_plot = "no"
motion_plot = "no" 
temp_plot = "yes"
temp_plot_log = "yes"
speed_plot = "no"
brownian_motion_plot = "no"

voronoi_plot = "yes"
delaunay_plot = "no"
g_plot = "no"
layer_final_pos_plot = "yes"

#%%

DATA_NAME = "Loaded_Periodic_Dust_grain_max_80_Final_Temperature_-nan_frames_210002"#input("Data name?(NOT FILENAME)")#"Finite_Dust_grain_max_5_Final_Temperature_48.356239_frames_501"
FILENAME = "HPC_Data/"+ DATA_NAME + ".csv"
#input("FILE name?")#"HPC_Data/Finite_Dust_grain_max_5_Final_Temperature_48.356239_frames_501.csv"
NEW_FOLDER = "Figures/" + DATA_NAME
if str(os.path.exists(NEW_FOLDER)) == "False":
   os.mkdir(NEW_FOLDER)

#%%
data = pd.read_csv(FILENAME)
print("LOADED...")
dust_grain_max = int( (len(data.columns))/3 - 1)
last_time_val = data["Time_list"].iloc[-1]
frames_len = len(data["Time_list"])
y_z_se = [z_se/lambda_D]*frames_len
temp_ion_a = [T_n_a]*frames_len
temp_ion_b = [T_n_b]*frames_len

#%%
if melting_plots_g ==  "yes":

    plt.figure()
    plt.xlabel(r"$r/ (\Delta r)$")
    plt.ylabel(r"g(r)")

    r_crystal = ((2**0.5)*container_length/2)/lambda_D
    #r_crystal = container_radius/lambda_D
    print("r_crystal = " ,r_crystal)

    for j in np.arange(len(temp_list)):
        layer_data = []
        temp_index = int(round(temp_list[j]))
        print("HPPP = ",temp_index)

        layer_x = []
        layer_y = []

        for i in np.arange(dust_grain_max):
            #if ( (data["Z_" + str(i)].values[temp_index] >= layers[0] ) and (data["Z_" + str(i)].values[temp_index] <= layers[1])):
            layer_x.append(data["X_" + str(i)].iloc[temp_index])
            layer_y.append(data["Y_" + str(i)].iloc[temp_index])        

        #set dr to be 100th
        dr = r_crystal/100

        radial_values = np.arange(dr/2,r_crystal - dr/2 , dr)
        #print("yo: ", radial_values)

        print(len(layer_x))
        dust_density = (np.pi*r_crystal**2)/len(layer_x)

        #need to calculate the distances between all the particles, fucccckkkkkkkkkkk
        inter_grain_distance_list = []
        for t in np.arange(len(layer_x)):
            inter_grain_distance_row = []
            for k in np.arange(len(layer_x)):
                r = ((layer_x[t] - layer_x[k])**2 + (layer_y[t] - layer_y[k])**2)**0.5
                inter_grain_distance_row.append(r)
            inter_grain_distance_list.append(inter_grain_distance_row)
        

        g = []
        for i in np.arange(len(radial_values)):
            N = 0
            for t in np.arange(len(layer_x)):
                for q in np.arange(len(layer_x)):
                    if ((inter_grain_distance_list[t][q] > radial_values[i] - dr/2) & ((inter_grain_distance_list[t][q] <= radial_values[i] + dr/2 ))):
                        N += 1
            g.append(N/(2*np.pi*(radial_values[i])*dr*dust_density))

        plt.plot(radial_values,g, label = "T =" + str(round(data["Temperature_list"].values[temp_index],3)) + "K, t = " + str(round(data["Time_list"].values[temp_index],4)) + "s" )

    plt.grid()
    plt.legend()
    plt.savefig(NEW_FOLDER + "/Melting_Pair_correlation" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

if melting_plots_snap ==  "yes":
        for j in np.arange(len(temp_list)):
            layer_data = []
            temp_index = int(round(temp_list[j]))

            layer_x = []
            layer_y = []

            for i in np.arange(dust_grain_max):
                #if ( (data["Z_" + str(i)].values[temp_index] >= layers[0] ) and (data["Z_" + str(i)].values[temp_index] <= layers[1])):
                layer_x.append(data["X_" + str(i)].iloc[temp_index])
                layer_y.append(data["Y_" + str(i)].iloc[temp_index])

            plt.figure()
            plt.title("Final Positions: T = " + str(round(data["Temperature_list"].values[temp_index],3)))
            plt.plot(layer_x, layer_y,"+" ,color='red')
            plt.xlabel(r"$x/\lambda_D$")
            plt.ylabel(r"$y/\lambda_D$")
            plt.grid()
            plt.xlim(-container_radius/(lambda_D),container_radius/(lambda_D))
            plt.ylim(-container_radius/(lambda_D),container_radius/(lambda_D))   

if pos_plot_six == "yes":

    n_rows = 2
    n_cols = 3
    fig, axs = plt.subplots(n_rows, n_cols)
    
    temp_list_data =  []
    for j in np.arange(len(temp_list)):
        layer_data = []
        temp_index = int(round(temp_list[j]))

        layer_x = []
        layer_y = []

        for i in np.arange(dust_grain_max):
            #if ( (data["Z_" + str(i)].values[temp_index] >= layers[0] ) and (data["Z_" + str(i)].values[temp_index] <= layers[1])):
            layer_x.append(data["X_" + str(i)].iloc[temp_index])
            layer_y.append(data["Y_" + str(i)].iloc[temp_index])
        temp_list_data.append([layer_x,layer_y])


    axs[0, 0].text(0.05, 0.95, "a)", transform=axs[0, 0].transAxes, fontsize=12,  va='top')
    axs[0, 0].plot(temp_list_data[0][0], temp_list_data[0][1],"+" ,color='red')
    axs[0, 0].set_xlabel(r"$x/\lambda_D$")
    axs[0, 0].set_ylabel(r"$y/\lambda_D$")
    axs[0, 0].grid()

    axs[0, 0].text(0.05, 0.95, "a)", transform=axs[0, 0].transAxes, fontsize=12,  va='top')
    axs[0, 0].plot(temp_list_data[0][0], temp_list_data[0][1],"+" ,color='red')
    axs[0, 0].set_xlabel(r"$x/\lambda_D$")
    axs[0, 0].set_ylabel(r"$y/\lambda_D$")
    axs[0, 0].grid()

    #axs[0, 1].set_title("z, Dust grains =  " + str(dust_grain_max_list[1]))
    axs[0, 1].text(0.05, 0.95, "b)", transform=axs[0, 1].transAxes, fontsize=12,  va='top')
    axs[0, 1].set_xlabel(r"$x/\lambda_D$")
    axs[0, 1].set_ylabel(r"$y/\lambda_D$")
    axs[0, 1].plot(temp_list_data[1][0], temp_list_data[1][1],"+" ,color='red')
    axs[0, 1].grid()

    axs[0, 2].text(0.05, 0.95, "c)", transform=axs[0, 2].transAxes, fontsize=12,  va='top')
    axs[0, 2].set_xlabel(r"$x/\lambda_D$")
    axs[0, 2].set_ylabel(r"$y/\lambda_D$")
    axs[0, 2].plot(temp_list_data[2][0], temp_list_data[2][1],"+" ,color='red')
    axs[0, 2].grid()

    axs[1, 0].text(0.05, 0.95, "d)", transform=axs[1, 0].transAxes, fontsize=12,  va='top')
    axs[1, 0].plot(temp_list_data[3][0], temp_list_data[3][1],"+" ,color='red')
    axs[1, 0].set_xlabel(r"$x/\lambda_D$")
    axs[1, 0].set_ylabel(r"$y/\lambda_D$")
    axs[1, 0].grid()

    axs[1, 1].text(0.05, 0.95, "e)", transform=axs[1, 1].transAxes, fontsize=12,  va='top')
    axs[1, 1].plot(temp_list_data[4][0], temp_list_data[4][1],"+" ,color='red')
    axs[1, 1].set_xlabel(r"$x/\lambda_D$")
    axs[1, 1].set_ylabel(r"$y/\lambda_D$")
    axs[1, 1].grid()

    axs[1, 2].text(0.05, 0.95, "f)", transform=axs[1, 2].transAxes, fontsize=12,  va='top')
    axs[1, 2].plot(temp_list_data[5][0], temp_list_data[5][1],"+" ,color='red')
    axs[1, 2].set_xlabel(r"$x/\lambda_D$")
    axs[1, 2].set_ylabel(r"$y/\lambda_D$")
    axs[1, 2].grid()
    
    plt.tight_layout()
    plt.savefig(NEW_FOLDER + "/snap_six_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

if temp_log_avoid == "yes":
    plt.figure()
    plt.plot(data["Time_list"][:temp_lim_index],data["Temperature_list"][:temp_lim_index])
    plt.xlabel(r"Time/s")
    plt.ylabel(r"Temperature/K")
    plt.yscale("log")
    plt.grid()
    plt.savefig(NEW_FOLDER + "/Temperature_lim_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)
plt.show()