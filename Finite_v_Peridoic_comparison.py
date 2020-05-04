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
#hist_div_val = 10 #larger less bins, smaller more bins
colour_list = ["red","green","blue","orange","black","brown"]
Brown_num = 1000
d_t = 10**(-4)
#%%

boundry_1 = "Periodic"#input("Periodic or Finite?")
boundry_2 = "Finite"#input("Periodic or Finite?")

#temp_plot = "yes"
temp_plot_log = "yes"
speed_plot = "no"

#%%
DATA_NAME_1 = "Periodic_Dust_grain_max_20_Tn_EV_0.030000_Final_Temperature_419.737211_frames_100002"#input("Data name_1?[PERIODIC](NOT FILENAME)")
DATA_NAME_2 = "Finite_Dust_grain_max_20_Tn_EV_0.030000_Final_Temperature_394.866634_frames_100002"#input("Data name_2?[FINITE](NOT FILENAME)")
FILENAME_1 = "HPC_Data/"+ DATA_NAME_1 + ".csv"
FILENAME_2 = "HPC_Data/"+ DATA_NAME_2 + ".csv"
NEW_FOLDER = "Figures/" + "Compare_test_are_you_sure"# + DATA_NAME_1 + "_and_" + DATA_NAME_2
if str(os.path.exists(NEW_FOLDER)) == "False":
    #os.mkdir("NEW_FOLDER")
    os.makedirs(NEW_FOLDER)

#%%
data_1 = pd.read_csv(FILENAME_1)
data_2 = pd.read_csv(FILENAME_2)
dust_grain_max = int( (len(data_1.columns))/3 - 1)
#print(dust_grain_max)
#last_time_val = data["Time_list"].iloc[-1]
#y_z_se = [z_se/lambda_D]*len(data["Time_list"])
temp_ion_a = [T_n_a]*len(data_1["Time_list"])
N = 100000
#print(T_n_a)
#temp_ion_b = [T_n_b]*len(data["Time_list"])
#%%
# if temp_plot == "yes":
#     plt.figure()
#     plt.title("Temperature Comparison")
#     plt.plot(data_1["Time_list"],data_1["Temperature_list"],color = "blue" , label = "Periodic")
#     plt.plot(data_2["Time_list"],data_2["Temperature_list"],color = "red"  , label = "Finite")
#     plt.plot(data_1["Time_list"],temp_ion_a, "--", color = "black",label =  "Neutral")
#     plt.xlabel("Time")
#     plt.ylabel("Temperature(K)")
#     plt.legend()
#     plt.grid()
#     plt.savefig(NEW_FOLDER + "/Comparison_Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data_1["Time_list"])) + ".png", dpi = 1200)

if temp_plot_log == "yes":
    plt.figure()
    plt.title("Log Temperature Comparison")
    plt.plot(data_1["Time_list"],data_1["Temperature_list"],color = "blue" , label = "Periodic")
    plt.plot(data_2["Time_list"],data_2["Temperature_list"],color = "red"  , label = "Finite")
    plt.plot(data_1["Time_list"],temp_ion_a, "--", color = "black",label =  "Neutral")
    plt.yscale("log")
    plt.xlabel("Time")
    plt.ylabel("Temperature(K)")
    plt.grid()
    plt.legend()
    # # plt.savefig(NEW_FOLDER + "/Comparison_Temperature_log_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data_1["Time_list"])) + ".png", dpi = 1200)

    fig, ax = plt.subplots() # create a new figure with a default 111 subplot
    ax.plot(data_1["Time_list"], data_1["Temperature_list"],color = "blue" , label = "Periodic")
    ax.plot(data_2["Time_list"], data_2["Temperature_list"],color = "red"  , label = "Finite")
    ax.plot(data_1["Time_list"],temp_ion_a, "--", color = "black",label =  "Neutral")
    axins = zoomed_inset_axes(ax, 12, loc=9) # zoom-factor: 2.5, location: upper-left
    axins.plot(data_1["Time_list"], data_1["Temperature_list"],color = "blue" , label = "Periodic")
    axins.plot(data_2["Time_list"], data_2["Temperature_list"],color = "red"  , label = "Finite")
    axins.plot(data_1["Time_list"],temp_ion_a, "--", color = "black",label =  "Neutral")
    x1, x2, y1, y2 = 9.5, 10, 300, 500 # specify the limits
    axins.set_yscale('log')
    ax.set_yscale('log')
    axins.set_xlim(x1, x2) # apply the x-limits
    axins.set_ylim(y1, y2) # apply the y-limits
    ax.grid()
    ax.legend()
    ax.set_xlabel("Time")
    ax.set_ylabel("Temperature(K)")
    #plt.yticks(visible=False)
    #plt.xticks(visible=False)
    mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")
    plt.savefig(NEW_FOLDER + "/Comparison_Temperature_log_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data_1["Time_list"])) + ".png", dpi = 1200)

    #whats the wacky behaviour?
    # x = np.linspace(0.0, N*d_t, N)
    temp_data_fourier_1 = data_1["Temperature_list"].values[-N:]
    temp_data_fourier_2 = data_2["Temperature_list"].values[-N:]
    temp_data_fourier_1_start = data_1["Temperature_list"].values[:N]
    temp_data_fourier_2_start = data_2["Temperature_list"].values[:N]
    sig_fft_1 = fft(temp_data_fourier_1)
    sig_fft_2 = fft(temp_data_fourier_2)
    sig_fft_1_start = fft(temp_data_fourier_1_start)
    sig_fft_2_start = fft(temp_data_fourier_2_start)
    # And the power (sig_fft is of complex dtype)
    power_1 = np.abs(sig_fft_1)*d_t
    power_2 = np.abs(sig_fft_2)*d_t
    power_1_start = np.abs(sig_fft_1_start)*d_t
    power_2_start = np.abs(sig_fft_2_start)*d_t
    # The corresponding frequencies
    sample_freq = fftfreq(N, d=d_t)

    # Plot the FFT power
    plt.figure()
    plt.plot(sample_freq, power_1,color = "blue" , label = "Periodic: End")
    plt.plot(sample_freq, power_2,color = "red"  , label = "Finite: End")
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power')
    plt.grid()
    plt.legend()

    plt.figure()
    plt.plot(sample_freq, power_1_start,color = "blue" , label = "Periodic: Start")
    plt.plot(sample_freq, power_2_start,color = "red"  , label = "Finite: Start")
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power')
    plt.grid()
    plt.legend()

if speed_plot == "yes":
    # (0, 4.9808603062591041)
    #plt.hist(data, bins=20, normed=True)


    # plt.figure()
    # plt.title("Speed Comparison")
    # plt.hist(data_1["Speed_list"].dropna(), bins = 20, density=True, color = "blue" , label = "Periodic")
    # plt.hist(data_2["Speed_list"].dropna(), bins = 20, density=True, color = "red"  , label = "Finite")
    # #x = np.linspace(0, 25, 100)
    # #plt.plot(x, maxwell.pdf(x, *params), lw=3,color = "black",label =  "Neutral")
    # #plt.plot(x, maxwell.pdf(x, *params_1), lw=3, color = "blue" , label = "Periodic fit")
    # #plt.plot(x, maxwell.pdf(x, *params_2), lw=3, color = "red"  , label = "Finite fit")
    # plt.ylabel("Frequency")
    # plt.xlabel("Speed")
    # plt.grid()
    # plt.legend()
    # plt.savefig(NEW_FOLDER + "/Comparison_Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data_1["Time_list"])) + ".png", dpi = 1200)

    # plt.figure()
    # plt.title("Speed")
    # plt.hist(data_1["Speed_list"].dropna(), bins = 20, density=True, color = "blue" , label = "Periodic")
    # x = np.linspace(0, 25, 100)
    # plt.plot(x, maxwell.pdf(x, *params), lw=3,color = "black",label =  "Neutral")
    # plt.plot(x, maxwell.pdf(x, *params_1), lw=3, color = "blue" , label = "Periodic fit")
    # plt.ylabel("Frequency")
    # plt.xlabel("Speed")
    # plt.grid()
    # plt.legend()

    maxwell = stats.maxwell
    data = maxwell.rvs(loc=0, scale= (k_b*T_n_a/m_D)**0.5, size=1000)
    params = maxwell.fit(data, floc=0)
    #params_1 = maxwell.fit(data_1["Speed_list"].dropna(), floc=0)
    #params_2 = maxwell.fit(data_2["Speed_list"].dropna(), floc=0)
    #print(data_2["Speed_list"].values().to_list())
    plt.figure()
    plt.title("Speed")
    #n, bins_2, patches = plt.hist(data_2["Speed_list"], bins = 20, density=True, color = "red"  , label = "Finite")
    x = np.linspace(0, 25, 1000)
    plt.plot(x, maxwell.pdf(x, *params), lw=3,color = "black",label =  "Neutral")
    plt.xlabel("Speed")
    plt.ylabel("Frequency")
    plt.grid()
    plt.legend()



# if layer_plots == "yes":
#    layer_data = []
#     for v in np.arange(len(layers) - 1):
#         layer_row = []
#         layer_x = []
#         layer_y = []
#         for i in np.arange(dust_grain_max):
#         last_val_index = np.where(data["Time_list"] == last_time_val)
#         if ( (data["Z_" + str(i)].values[-1][0] >= layers[v] ) and (data["Z_" + str(i)].values[-1][0] <= layers[v+1])):
#             print("YOOOOO: ",data["Z_" + str(i)].values[-1][0],data["Z_" + str(i)].values[-1][0],data["X_" + str(i)].iloc[last_val_index].to_numpy())
#             layer_x.append(data["X_" + str(i)].iloc[last_val_index].to_numpy())
#             layer_y.append(data["Y_" + str(i)].iloc[last_val_index].to_numpy())
#         layer_row.append(layer_x)
#         layer_row.append(layer_y)
#         layer_data.append(layer_row)

#         points = np.c_[layer_row[0], layer_row[1]]

#         if g_plot == "yes":
#         pair_correlation_x_data = layer_row[0]
#         pair_correlation_y_data = layer_row[1]

#         #get max radial posiiton - bit shit but i cnat be botherd
#         r_crystal = abs(max(pair_correlation_x_data))[0]
#         #set dr to be 100th
#         dr = r_crystal/100
#         radial_values = np.arange(0,r_crystal, dr)
#         dust_density = np.pi*r_crystal**2/len(points)

#         #need to calculate the distances between all the particles, fucccckkkkkkkkkkk
#         inter_grain_distance_list = []
#         for t in np.arange(len(points)):
#             inter_grain_distance_row = []
#             for k in np.arange(len(points)):
#                 r = ((points[t][0] - points[k][1])**2 + (points[t][1] - points[k][1])**2)**0.5
#                 inter_grain_distance_row.append(r)
#             inter_grain_distance_list.append(inter_grain_distance_row)
            
#         g = []
#         for i in np.arange(len(radial_values) - 1):
#             N = 0
#             for t in np.arange(len(points)):
#                 for q in np.arange(len(points)):
#                     if ((inter_grain_distance_list[t][q] > radial_values[i]) & ((inter_grain_distance_list[t][q] <= radial_values[i+1]))):
#                     N += 1
#             g.append((N/(2*np.pi*radial_values[i+1]*dr))/dust_density)

#         g = np.asarray(g)
#         radial_values = np.delete(radial_values, 0)

#         plt.figure()
#         plt.title("Pair correlation in layer: " + str(v))
#         plt.xlabel("r/dr")
#         plt.ylabel("g")
#         plt.plot(radial_values,g)
#         plt.grid()
#         plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Pair_correlation" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png", dpi = 1200)

plt.show()