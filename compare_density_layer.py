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

#%%
N = 5
z_plot = "no"
layers_calc = "no"
z_plot_six = "yes"

#%%
DATA_NAME_1 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
DATA_NAME_2 = "Periodic_Dust_grain_max_60_Tn_EV_0.030000_Final_Temperature_405.860905_frames_100002"#input("Data name_2?[FINITE](NOT FILENAME)")
DATA_NAME_3 = "Periodic_Dust_grain_max_80_Tn_EV_0.030000_Final_Temperature_402.872104_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
DATA_NAME_4 = "Periodic_Dust_grain_max_100_Tn_EV_0.030000_Final_Temperature_405.202493_frames_100002"#input("Data name_2?[FINITE](NOT FILENAME)")
DATA_NAME_5 = "Periodic_Dust_grain_max_120_Tn_EV_0.030000_Final_Temperature_404.670559_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
DATA_NAME_6 = "Periodic_Dust_grain_max_180_Tn_EV_0.030000_Final_Temperature_399.836619_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")

# DATA_NAME_1 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
# DATA_NAME_2 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_2?[FINITE](NOT FILENAME)")
# DATA_NAME_3 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
# DATA_NAME_4 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_2?[FINITE](NOT FILENAME)")
# DATA_NAME_5 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
# DATA_NAME_6 = "Periodic_Dust_grain_max_40_Tn_EV_0.030000_Final_Temperature_407.929253_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")

FILENAME_1 = "HPC_Data/"+ DATA_NAME_1 + ".csv"
FILENAME_2 = "HPC_Data/"+ DATA_NAME_2 + ".csv"
FILENAME_3 = "HPC_Data/"+ DATA_NAME_3 + ".csv"
FILENAME_4 = "HPC_Data/"+ DATA_NAME_4 + ".csv"
FILENAME_5 = "HPC_Data/"+ DATA_NAME_5 + ".csv"
FILENAME_6 = "HPC_Data/"+ DATA_NAME_6 + ".csv"

NEW_FOLDER = "Figures/" + "Compare_various_dust_grains"# + DATA_NAME_1 + "_and_" + DATA_NAME_2
if str(os.path.exists(NEW_FOLDER)) == "False":
    #os.mkdir("NEW_FOLDER")
    os.makedirs(NEW_FOLDER)

#%%
data_1 = pd.read_csv(FILENAME_1)
data_2 = pd.read_csv(FILENAME_2)
data_3 = pd.read_csv(FILENAME_3)
data_4 = pd.read_csv(FILENAME_4)
data_5 = pd.read_csv(FILENAME_5)
data_6 = pd.read_csv(FILENAME_6)

data_list = [data_1,data_2,data_3,data_4,data_5,data_6]

dust_grain_max_1 = int( (len(data_1.columns))/3 - 1)
dust_grain_max_2 = int( (len(data_2.columns))/3 - 1)
dust_grain_max_3 = int( (len(data_3.columns))/3 - 1)
dust_grain_max_4 = int( (len(data_4.columns))/3 - 1)
dust_grain_max_5 = int( (len(data_5.columns))/3 - 1)
dust_grain_max_6 = int( (len(data_6.columns))/3 - 1)

dust_grain_max_list = [dust_grain_max_1,dust_grain_max_2,dust_grain_max_3,dust_grain_max_4,dust_grain_max_5,dust_grain_max_6]
frames_len = len(data_1["Time_list"])
y_z_se = [z_se/lambda_D]*len(data_1["Time_list"])
last_time_val = data_1["Time_list"].iloc[-1]

#%%
if z_plot == "yes":
    
    for v in np.arange(len(dust_grain_max_list)):
        plt.figure()
        plt.title("z, Dust grains =  " + str(dust_grain_max_list[v]))
        for i in np.arange(dust_grain_max_list[v]):
            last_val_index = np.where(data_list[v]["Time_list"] == last_time_val)
            plt.plot(data_list[v]["Time_list"], data_list[v]["Z_" + str(i)])
            plt.plot(data_list[v]["Time_list"][0], data_list[v]["Z_" + str(i)][0],"+" ,color='blue')
            plt.plot(data_list[v]["Time_list"].iloc[last_val_index], data_list[v]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
        plt.xlabel(r"Time/s")
        plt.ylabel(r"$z/\lambda_D$")
        plt.plot(data_list[v]["Time_list"],y_z_se, "--", color = "black")
        plt.ylim(0,z_se*top_container_graph_mul/lambda_D)
        plt.grid()
        #plt.savefig(NEW_FOLDER + "/Periodic_z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)
#%%
if z_plot_six == "yes":

    n_rows = 2
    n_cols = 3
    fig, axs = plt.subplots(n_rows, n_cols)

    #axs[0, 0].set_title("z, Dust grains =  " + str(dust_grain_max_list[0]))
    axs[0, 0].text(0.05, 0.95, "a)", transform=axs[0, 0].transAxes, fontsize=16, fontweight='bold', va='top')
    for i in np.arange(dust_grain_max_list[0]):
        last_val_index = np.where(data_list[0]["Time_list"] == last_time_val)
        axs[0, 0].plot(data_list[0]["Time_list"], data_list[0]["Z_" + str(i)])
        axs[0, 0].plot(data_list[0]["Time_list"][0], data_list[0]["Z_" + str(i)][0],"+" ,color='blue')
        axs[0, 0].plot(data_list[0]["Time_list"].iloc[last_val_index], data_list[0]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[0, 0].set_xlabel(r"Time/s")
    axs[0, 0].set_ylabel(r"$z/\lambda_D$")
    axs[0, 0].plot(data_list[0]["Time_list"],y_z_se, "--", color = "black")
    axs[0, 0].set_ylim(0,z_se*top_container_graph_mul/lambda_D)
    axs[0, 0].grid()

    #axs[0, 1].set_title("z, Dust grains =  " + str(dust_grain_max_list[1]))
    axs[0, 1].text(0.05, 0.95, "b)", transform=axs[0, 1].transAxes, fontsize=16, fontweight='bold', va='top')
    for i in np.arange(dust_grain_max_list[1]):
        last_val_index = np.where(data_list[1]["Time_list"] == last_time_val)
        axs[0, 1].plot(data_list[1]["Time_list"], data_list[1]["Z_" + str(i)])
        axs[0, 1].plot(data_list[1]["Time_list"][0], data_list[1]["Z_" + str(i)][0],"+" ,color='blue')
        axs[0, 1].plot(data_list[1]["Time_list"].iloc[last_val_index], data_list[1]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[0, 1].set_xlabel(r"Time/s")
    axs[0, 1].set_ylabel(r"$z/\lambda_D$")
    axs[0, 1].plot(data_list[1]["Time_list"],y_z_se, "--", color = "black")
    axs[0, 1].set_ylim(0,z_se*top_container_graph_mul/lambda_D)
    axs[0, 1].grid()

    axs[0, 2].text(0.05, 0.95, "c)", transform=axs[0, 2].transAxes, fontsize=16, fontweight='bold', va='top')
    for i in np.arange(dust_grain_max_list[2]):
        last_val_index = np.where(data_list[2]["Time_list"] == last_time_val)
        axs[0, 2].plot(data_list[2]["Time_list"], data_list[2]["Z_" + str(i)])
        axs[0, 2].plot(data_list[2]["Time_list"][0], data_list[2]["Z_" + str(i)][0],"+" ,color='blue')
        axs[0, 2].plot(data_list[2]["Time_list"].iloc[last_val_index], data_list[2]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[0, 2].set_xlabel(r"Time/s")
    axs[0, 2].set_ylabel(r"$z/\lambda_D$")
    axs[0, 2].plot(data_list[2]["Time_list"],y_z_se, "--", color = "black")
    axs[0, 2].set_ylim(0,z_se*top_container_graph_mul/lambda_D)
    axs[0, 2].grid()

    axs[1, 0].text(0.05, 0.95, "d)", transform=axs[1, 0].transAxes, fontsize=16, fontweight='bold', va='top')
    for i in np.arange(dust_grain_max_list[3]):
        last_val_index = np.where(data_list[3]["Time_list"] == last_time_val)
        axs[1, 0].plot(data_list[3]["Time_list"], data_list[3]["Z_" + str(i)])
        axs[1, 0].plot(data_list[3]["Time_list"][0], data_list[3]["Z_" + str(i)][0],"+" ,color='blue')
        axs[1, 0].plot(data_list[3]["Time_list"].iloc[last_val_index], data_list[3]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[1, 0].set_xlabel(r"Time/s")
    axs[1, 0].set_ylabel(r"$z/\lambda_D$")
    axs[1, 0].plot(data_list[3]["Time_list"],y_z_se, "--", color = "black")
    axs[1, 0].set_ylim(0,z_se*top_container_graph_mul/lambda_D)
    axs[1, 0].grid()

    axs[1, 1].text(0.05, 0.95, "e)", transform=axs[1, 1].transAxes, fontsize=16, fontweight='bold', va='top')
    for i in np.arange(dust_grain_max_list[4]):
        last_val_index = np.where(data_list[4]["Time_list"] == last_time_val)
        axs[1, 1].plot(data_list[4]["Time_list"], data_list[4]["Z_" + str(i)])
        axs[1, 1].plot(data_list[4]["Time_list"][0], data_list[4]["Z_" + str(i)][0],"+" ,color='blue')
        axs[1, 1].plot(data_list[4]["Time_list"].iloc[last_val_index], data_list[4]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[1, 1].set_xlabel(r"Time/s")
    axs[1, 1].set_ylabel(r"$z/\lambda_D$")
    axs[1, 1].plot(data_list[4]["Time_list"],y_z_se, "--", color = "black")
    axs[1, 1].set_ylim(0,z_se*top_container_graph_mul/lambda_D)
    axs[1, 1].grid()

    axs[1, 2].text(0.05, 0.95, "f)", transform=axs[1, 2].transAxes, fontsize=16, fontweight='bold', va='top')
    for i in np.arange(dust_grain_max_list[5]):
        last_val_index = np.where(data_list[5]["Time_list"] == last_time_val)
        axs[1, 2].plot(data_list[5]["Time_list"], data_list[5]["Z_" + str(i)])
        axs[1, 2].plot(data_list[5]["Time_list"][0], data_list[5]["Z_" + str(i)][0],"+" ,color='blue')
        axs[1, 2].plot(data_list[5]["Time_list"].iloc[last_val_index], data_list[5]["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
    axs[1, 2].set_xlabel(r"Time/s")
    axs[1, 2].set_ylabel(r"$z/\lambda_D$")
    axs[1, 2].plot(data_list[5]["Time_list"],y_z_se, "--", color = "black")
    axs[1, 2].set_ylim(0,z_se*top_container_graph_mul/lambda_D)
    axs[1, 2].grid()
    
    plt.tight_layout()

    plt.savefig(NEW_FOLDER + "/comaparing_z_Periodic_" + "_frames_" + str(frames_len) + ".png", dpi = 1200)

#%%
if layers_calc == "yes":
    layers_1 = [0.0,50.0]
    layers_2 = [0.0,23.0,50.0]
    layers_3 = [0.0,23.0,50.0]
    layers_4 = [0.0,23.0,50.0]
    layers_5 = [0.0,50.0]
    layers_list = [layers_1,layers_2,layers_3,layers_4,layers_5]

    layer_num_all = []

    for j in np.arange(len(dust_grain_max_list)):
        layer_num_specific_cyrstal = []
        for v in np.arange(len(layers_list[j]) - 1):
            layer_num = 0
            for i in np.arange(dust_grain_max_list[j]):
                last_val_index = np.where(data_list[j]["Time_list"] == last_time_val)
                if ((data_list[j]["Z_" + str(i)].values[last_val_index][0] >= layers_list[j][v] ) and (data_list[j]["Z_" + str(i)].values[last_val_index][0] <= layers_list[j][v+1])):
                    layer_num += 1
            layer_num_specific_cyrstal.append(layer_num)
        layer_num_all.append(layer_num_specific_cyrstal)

    print(layer_num_all)

#%%
plt.show()