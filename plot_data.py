import numpy as np
import matplotlib.pyplot as plt
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
boundry = "Finite"#input("Periodic or Finite?")
standard_plots = "no"
layer_plots = "no"#input("layers?")
threeD_plot = "no"
anim_mono = "yes"
anim_dual = "no"
anim_dual_overlay = "no"
layers = [0.0,50.0]

x_plot = "no"
y_plot = "no"
z_plot = "no"
final_pos_plot = "no"
motion_plot = "no" 
temp_plot = "no"
temp_plot_log = "no"
speed_plot = "no"
brownian_motion_plot = "yes"

voronoi_plot = "yes"
delaunay_plot = "no"
g_plot = "no"
layer_final_pos_plot = "no"

if status == "Compile":
   if boundry == "Periodic":
      if load == "yes":
         subprocess.call(["g++", "-o", "Load_MSci_project_periodic_HPC", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/Load_MSci_project_periodic_HPC.cpp"])
         """Import stuff to measure how long the code takes to run"""
         start_time = time.time()
         print("start_time =", time.ctime(time.time()))
         subprocess.call("Load_MSci_project_periodic_HPC.exe")
         """prints time taken in minutes"""
         print ("time taken: %s minutes" % ((time.time()-start_time)/60))
      else:
         subprocess.call(["g++", "-o", "MSci_project_periodic_HPC", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/MSci_project_periodic_HPC.cpp"])
         """Import stuff to measure how long the code takes to run"""
         start_time = time.time()
         print("start_time =", time.ctime(time.time()))
         subprocess.call("MSci_project_periodic_HPC.exe")
         """prints time taken in minutes"""
         print ("time taken: %s minutes" % ((time.time()-start_time)/60))
   elif boundry == "Finite":
      if load == "yes":
         subprocess.call(["g++", "-o", "Load_MSci_project_HPC", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/Load_MSci_project_HPC.cpp"])
         """Import stuff to measure how long the code takes to run"""
         start_time = time.time()
         print("start_time =", time.ctime(time.time()))
         subprocess.call("Load_MSci_project_HPC.exe")
         """prints time taken in minutes"""
         print ("time taken: %s minutes" % ((time.time()-start_time)/60))
      else:
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

DATA_NAME = "Finite_Dust_grain_max_100_Tn_EV_0.030000_Final_Temperature_400.139865_frames_100002"#input("Data name?(NOT FILENAME)")#"Finite_Dust_grain_max_5_Final_Temperature_48.356239_frames_501"
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
if standard_plots == "yes":
   if boundry == "Periodic":
      if x_plot == "yes":
         plt.figure()
         plt.title("test - x")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["X_" + str(i)])
            plt.plot(data["Time_list"][0], data["X_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["X_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"Time /s")
         plt.ylabel("x/lambda_D")
         plt.grid()
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.savefig(NEW_FOLDER + "/Periodic_x_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if y_plot == "yes":
         plt.figure()
         plt.title("test - y")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Y_" + str(i)])
            plt.plot(data["Time_list"][0], data["Y_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"Time /s")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.savefig(NEW_FOLDER + "/Periodic_y_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if z_plot == "yes":
         plt.figure()
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Z_" + str(i)])
            plt.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"Time /s")
         plt.ylabel(r"z /$\lambda_D$")
         plt.plot(data["Time_list"],y_z_se, "--", color = "black")
         plt.xlim(0,11)
         plt.ylim(21,26)
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if final_pos_plot == "yes":
         plt.figure()
         plt.title("Final Positions")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("x/lambda_D")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.savefig(NEW_FOLDER + "/Periodic_final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if motion_plot == "yes":   
         plt.figure()
         plt.title("Motion")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"].values == last_time_val)
            plt.plot(data["X_" + str(i)],data["Y_" + str(i)])
            plt.plot(data["X_" + str(i)][0], data["Y_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("x/lambda_D")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))

      if brownian_motion_plot == "yes":   
         plt.figure()
         plt.title("Brownian Motion")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"].values == last_time_val)
            plt.plot(data["X_" + str(i)].values.tolist()[-Brown_num:],data["Y_" + str(i)].values.tolist()[-Brown_num:])
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("x/lambda_D")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.savefig(NEW_FOLDER + "/Periodic_Brownian_motion_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if temp_plot == "yes":
         plt.figure()
         plt.title("Temperature")
         plt.plot(data["Time_list"],data["Temperature_list"])
         plt.plot(data["Time_list"],temp_ion_a, "--", color = "black")
         plt.plot(data["Time_list"],temp_ion_b, "--", color = "black")
         plt.xlabel(r"Time /s")
         plt.ylabel("Temperature(K)")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if temp_plot_log == "yes":
         plt.figure()
         plt.title("Temperature")
         plt.plot(data["Time_list"],data["Temperature_list"])
         plt.plot(data["Time_list"],temp_ion_a, "--", color = "black")
         plt.plot(data["Time_list"],temp_ion_b, "--", color = "black")
         plt.yscale("log")
         plt.xlabel(r"Time /s")
         plt.ylabel("Temperature(K)")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_Temperature_log_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if speed_plot == "yes":
         plt.figure()
         plt.title("Speed")
         plt.hist(data["Speed_list"].dropna(), bins = 20  )
         plt.xlabel("Speed")
         plt.ylabel("Frequency")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

   else:
      theta = np.linspace(0, 2*np.pi, 100)
      x_r_wall = container_radius/lambda_D*np.cos(theta)
      y_r_wall = container_radius/lambda_D*np.sin(theta)

      if x_plot == "yes":
         plt.figure()
         plt.title("test - x")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["X_" + str(i)])
            plt.plot(data["Time_list"][0], data["X_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["X_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"Time /s")
         plt.ylabel("x/lambda_D")
         plt.grid()
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_x_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if y_plot == "yes":
         plt.figure()
         plt.title("test - y")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Y_" + str(i)])
            plt.plot(data["Time_list"][0], data["Y_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"Time /s")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_y_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if z_plot == "yes":
         plt.figure()
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Z_" + str(i)])
            plt.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"Time /s")
         plt.ylabel(r"z /$\lambda_D$")
         plt.plot(data["Time_list"],y_z_se, "--", color = "black")
         plt.ylim(0,z_se*top_container_graph_mul/lambda_D)
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Finite_z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if final_pos_plot == "yes":
         plt.figure()
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"$x/\lambda_D$")
         plt.ylabel(r"$y/\lambda_D$")
         plt.plot(x_r_wall,y_r_wall, color = "black")
         plt.grid()
         plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if motion_plot == "yes":   
         plt.figure()
         plt.title("Motion")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"].values == last_time_val)
            plt.plot(data["X_" + str(i)],data["Y_" + str(i)])
            plt.plot(data["X_" + str(i)][0], data["Y_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("x/lambda_D")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_motion_plot_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if brownian_motion_plot == "yes":   
         plt.figure()
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"].values == last_time_val)
            plt.plot(data["X_" + str(i)].values.tolist()[-Brown_num:],data["Y_" + str(i)].values.tolist()[-Brown_num:])
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel(r"$x/\lambda_D$")
         plt.ylabel(r"$y/\lambda_D$")
         plt.grid()
         plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_Brownian_motion_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)


      if temp_plot == "yes":
         plt.figure()
         plt.title("Temperature")
         plt.plot(data["Time_list"],data["Temperature_list"])
         plt.plot(data["Time_list"],temp_ion_a, "--", color = "black")
         plt.plot(data["Time_list"],temp_ion_b, "--", color = "black")
         plt.yscale("log")
         plt.xlabel(r"Time /s")
         plt.ylabel("Temperature(K)")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Finite_Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

      if speed_plot == "yes":
         plt.figure()
         plt.title("Speed")
         plt.hist(data["Speed_list"].dropna(), bins = 20  )
         plt.xlabel("Speed")
         plt.ylabel("Frequency")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Finite_Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

if layer_plots == "yes":
   layer_data = []

   if boundry == "Periodic":
      for v in np.arange(len(layers) - 1):
         layer_row = []
         layer_x = []
         layer_y = []
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            if ( (data["Z_" + str(i)].values[last_val_index][0] >= layers[v] ) and (data["Z_" + str(i)].values[last_val_index][0] <= layers[v+1])):
               #print(data["X_" + str(i)].iloc[last_val_index])
               layer_x.append(data["X_" + str(i)].iloc[last_val_index].to_numpy())
               layer_y.append(data["Y_" + str(i)].iloc[last_val_index].to_numpy())
         layer_row.append(layer_x)
         layer_row.append(layer_y)
         layer_data.append(layer_row)

         points = np.c_[layer_row[0], layer_row[1]]

         if voronoi_plot == "yes":
            vor = Voronoi(points)
            #fig = voronoi_plot_2d(vor)
            fig_v = plt.figure()
            ax_v = fig_v.add_subplot(111)
            plt.xlabel(r"$x/\lambda_D$")
            plt.ylabel(r"$y/\lambda_D$")
            voronoi_plot_2d(vor, ax = ax_v)
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Voronoi" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

         if delaunay_plot == "yes":
            tri = Delaunay(points)
            plt.figure()
            plt.title("Delaunay_Triangulation layer: " + str(v))
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.triplot(points[:,0], points[:,1], tri.simplices)
            plt.plot(points[:,0], points[:,1], 'o')
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Delaunay_Triangulation" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

         if g_plot == "yes":
            pair_correlation_x_data = layer_row[0]
            pair_correlation_y_data = layer_row[1]

            #get max radial posiiton - bit shit but i cnat be botherd
            r_crystal = abs(max(pair_correlation_x_data))[0]
            #set dr to be 100th
            dr = r_crystal/100
            radial_values = np.arange(0,r_crystal, dr)
            dust_density = np.pi*r_crystal**2/len(points)

            #need to calculate the distances between all the particles, fucccckkkkkkkkkkk
            inter_grain_distance_list = []
            for t in np.arange(len(points)):
               inter_grain_distance_row = []
               for k in np.arange(len(points)):
                  r = ((points[t][0] - points[k][1])**2 + (points[t][1] - points[k][1])**2)**0.5
                  inter_grain_distance_row.append(r)
               inter_grain_distance_list.append(inter_grain_distance_row)
               
            g = []
            for i in np.arange(len(radial_values) - 1):
               N = 0
               for t in np.arange(len(points)):
                  for q in np.arange(len(points)):
                     if ((inter_grain_distance_list[t][q] > radial_values[i]) & ((inter_grain_distance_list[t][q] <= radial_values[i+1]))):
                        N += 1
               g.append((N/(2*np.pi*radial_values[i+1]*dr))/dust_density)

            g = np.asarray(g)
            radial_values = np.delete(radial_values, 0)

            plt.figure()
            plt.title("Pair correlation in layer: " + str(v))
            plt.xlabel("r/dr")
            plt.ylabel("g")
            plt.plot(radial_values,g)
            plt.grid()
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Pair_correlation" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

         if layer_final_pos_plot == "yes":
            plt.figure()
            plt.title("Final Positions layer: " + str(v))
            plt.plot(layer_row[0],layer_row[1],"+" ,color='red')
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.grid()
            plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
            plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

   else:
      theta = np.linspace(0, 2*np.pi, 100)
      x_r_wall = container_radius/lambda_D*np.cos(theta)
      y_r_wall = container_radius/lambda_D*np.sin(theta)

      for v in np.arange(len(layers) - 1):
         layer_row = []
         layer_x = []
         layer_y = []
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            if ( (data["Z_" + str(i)].values[last_val_index][0] >= layers[v] ) and (data["Z_" + str(i)].values[last_val_index][0] <= layers[v+1])):
               #print(data["X_" + str(i)].iloc[last_val_index])
               layer_x.append(data["X_" + str(i)].iloc[last_val_index].to_numpy())
               layer_y.append(data["Y_" + str(i)].iloc[last_val_index].to_numpy())
         layer_row.append(layer_x)
         layer_row.append(layer_y)
         layer_data.append(layer_row)
         #print("hi")
         #print(layer_row)
         points = np.c_[layer_row[0], layer_row[1]]

         if voronoi_plot == "yes":
            vor = Voronoi(points)
            #fig = voronoi_plot_2d(vor)
            fig_v = plt.figure()
            ax_v = fig_v.add_subplot(111)
            plt.xlabel(r"$x/\lambda_D$")
            plt.ylabel(r"$y/\lambda_D$")
            voronoi_plot_2d(vor, ax = ax_v)
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "Voronoi" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

         if delaunay_plot == "yes":
            tri = Delaunay(points)
            plt.figure()
            plt.title("Delaunay_Triangulation layer: " + str(v))
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.triplot(points[:,0], points[:,1], tri.simplices)
            plt.plot(points[:,0], points[:,1], 'o')
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "Delaunay_Triangulation" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)
            
         if g_plot == "yes":
            pair_correlation_x_data = layer_row[0]
            pair_correlation_y_data = layer_row[1]

            #get max radial posiiton - bit shit but i cnat be botherd
            r_crystal = abs(max(pair_correlation_x_data))[0]
            #set dr to be 100th
            dr = r_crystal/100
            radial_values = np.arange(0,r_crystal, dr)
            dust_density = np.pi*r_crystal**2/len(points)

            #need to calculate the distances between all the particles, fucccckkkkkkkkkkk
            inter_grain_distance_list = []
            for t in np.arange(len(points)):
               inter_grain_distance_row = []
               for k in np.arange(len(points)):
                  r = ((points[t][0] - points[k][0])**2 + (points[t][1] - points[k][1])**2)**0.5
                  inter_grain_distance_row.append(r)
               inter_grain_distance_list.append(inter_grain_distance_row)
               
            g = []
            for i in np.arange(len(radial_values) - 1):
               N = 0
               for t in np.arange(len(points)):
                  for q in np.arange(len(points)):
                     if ((inter_grain_distance_list[t][q] > radial_values[i]) & ((inter_grain_distance_list[t][q] <= radial_values[i+1]))):
                        N += 1
               g.append((N/(2*np.pi*radial_values[i+1]*dr))/dust_density)

            g = np.asarray(g)
            radial_values = np.delete(radial_values, 0)

            plt.figure()
            plt.title("Pair correlation in layer: " + str(v))
            plt.xlabel("r/dr")
            plt.ylabel("g")
            plt.plot(radial_values,g)
            plt.grid()
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "Pair_correlation" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

         if layer_final_pos_plot == "yes":
            plt.figure()
            plt.title("Final Positions layer: " + str(v))
            plt.plot(layer_row[0],layer_row[1],"+" ,color='red')
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.grid()
            plt.xlim(-container_radius/(lambda_D),container_radius/(lambda_D))
            plt.ylim(-container_radius/(lambda_D),container_radius/(lambda_D))
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(frames_len) + ".png", dpi = 1200)

if threeD_plot == "yes":
   fig = plt.figure()
   ax = fig.add_subplot(111, projection='3d')
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list"] == last_time_val)
      for k in np.arange(len(layers) - 1):
         if ( (data["Z_" + str(i)].values[last_val_index][0] >= layers[k] ) and (data["Z_" + str(i)].values[last_val_index][0] <= layers[k+1])):
            colour_layer = colour_list[k]
            ax.scatter(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],data["Z_" + str(i)].iloc[last_val_index], color = colour_layer)
   if boundry == "Periodic":
      ax.set_xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
      ax.set_ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   else:
      ax.set_xlim(-container_radius/lambda_D,container_radius/lambda_D)
      ax.set_ylim(-container_radius/lambda_D,container_radius/lambda_D)
   ax.set_zlim(0,top_container_graph_mul*z_se/lambda_D)
   ax.set_xlabel(r"x/$\lambda_D$")
   ax.set_ylabel(r"y/$\lambda_D$")
   ax.set_zlabel(r"z/$\lambda_D$")

if anim_mono == "yes":
   size_mul = 1.2
   f_p_s = 30
   anim_time_lim = 60
   speed_mul = round(frames_len/(f_p_s*anim_time_lim))
   frames_anim = int(round(frames_len/speed_mul))

   fig, (ax1, ax3) = plt.subplots(1, 2, figsize =[10,5])

   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list"] == last_time_val)
      ax3.plot(data["Time_list"], data["Z_" + str(i)])
      ax3.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
      ax3.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
   ax3.plot(data["Time_list"],y_z_se, "--", color = "black")
   ax3.set_xlabel(r"Time")
   ax3.set_ylabel(r"z/$\lambda_D$")
   ax3.set_ylim(0,z_se*top_container_graph_mul/lambda_D)
   ax3.grid()
   time_line = ax3.axvline(x=0.0,linewidth=2, color='r')

   ax1.set_xlabel(r"x/$\lambda_D$")
   ax1.set_ylabel(r"y/$\lambda_D$")
   dust_points = ax1.plot([], [],"o", color = "blue")

   if boundry == "Periodic":
      ax1.set(xlim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul), ylim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul))
      rect1 = Rectangle((-container_length/(2*lambda_D), -container_length/(2*lambda_D)),width =  container_length/lambda_D, height = container_length/lambda_D, facecolor = "none", edgecolor="black", linewidth=1)
      ax1.add_patch(rect1)
      
      text = ax1.text(0 + 0.01, 0 + 0.01,'', transform = ax1.transAxes)

   else:
      ax1.set(xlim=(-container_radius/lambda_D, container_radius/lambda_D), ylim=(-container_radius/lambda_D, container_radius/lambda_D))
      circle1 = Circle((0, 0), container_radius/lambda_D,facecolor = "none", edgecolor="black", linewidth=1)
      ax1.add_patch(circle1)

      text = ax1.text(0 + 0.01, 0 + 0.01,'', transform = ax1.transAxes)

   plt.tight_layout()

   #need vector with all the x points and all the y points
   def get_data_mono(v):
      x = []
      y = []
      for i in np.arange(dust_grain_max):
            x.append(data["X_" + str(i)].iloc[v])
            y.append(data["Y_" + str(i)].iloc[v])

      time = data["Time_list"].iloc[v]
      return [x,y,time]


   def init_mono():
      """initialize animation"""
      dust_points[0].set_data([], [])
      text.set_text("Time = " + str(0) + "s")
      
      return [dust_points[0],text,time_line]

   def animate_mono(i):
      """perform animation step"""
      frame_num = i*speed_mul

      data_anim = get_data_mono(frame_num)
      
      # update pieces of the animation
      dust_points[0].set_data(data_anim[0],data_anim[1])
      text.set_text("Time = " + str(round(data_anim[2],5)) + "s")#updat value of the frame number
      time_line.set_data([data_anim[2], data_anim[2]], [0, 1])

      return [dust_points[0],text,time_line]

   #frames = round(frames/speed_mul) 
   #YOU HAVE TO MANUALLY SET THE LAYERS
   ani = animation.FuncAnimation(fig, animate_mono, interval=10, blit=True, frames = frames_anim, init_func=init_mono)
   FFwriter=animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
   ani.save('Finite_100_dust_grains_100000_frames_with_temperature.mp4', writer=FFwriter)

if anim_dual == "yes":
   size_mul = 1.2
   f_p_s = 30
   anim_time_lim = 30
   speed_mul = round(frames_len/(f_p_s*anim_time_lim))
   frames_anim = int(round(frames_len/speed_mul))

   fig = plt.figure(figsize =[7,7])#figsize = [8,5]
   ax1 = fig.add_subplot(2,2,1)
   ax2 = fig.add_subplot(2,2,2)
   ax3 = fig.add_subplot(2,1,2)

   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list"] == last_time_val)
      ax3.plot(data["Time_list"], data["Z_" + str(i)])
      ax3.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
      ax3.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
   ax3.plot(data["Time_list"],y_z_se, "--", color = "black")
   ax3.set_xlabel(r"Time")
   ax3.set_ylabel(r"z/$\lambda_D$")
   ax3.set_ylim(0,z_se*top_container_graph_mul/lambda_D)
   ax3.grid()
   time_line = ax3.axvline(x=0.0,linewidth=2, color='r')

   ax1.set_xlabel(r"x/$\lambda_D$")
   ax1.set_ylabel(r"y/$\lambda_D$")
   ax1.set_title(r"Lower Crystal Layer")
   ax2.set_xlabel(r"x/$\lambda_D$")
   ax2.set_ylabel(r"y/$\lambda_D$")
   ax2.set_title(r"Upper Crystal Layer")
   dust_points_lower = ax1.plot([], [],"o", color = "green")
   dust_points_upper = ax2.plot([], [],"o", color = "blue")

   if boundry == "Periodic":
      ax1.set(xlim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul), ylim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul))
      rect1 = Rectangle((-container_length/(2*lambda_D), -container_length/(2*lambda_D)),width =  container_length/lambda_D, height = container_length/lambda_D, facecolor = "none", edgecolor="black", linewidth=1)
      ax1.add_patch(rect1)
      
      ax2.set(xlim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul), ylim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul))
      rect2 = Rectangle((-container_length/(2*lambda_D), -container_length/(2*lambda_D)),width =  container_length/lambda_D, height = container_length/lambda_D, facecolor = "none", edgecolor="black", linewidth=1)
      ax2.add_patch(rect2)
      
      text = ax1.text(0 + 0.01, 0 + 0.01,'', transform = ax1.transAxes)

   else:
      ax1.set(xlim=(-container_radius/lambda_D, container_radius/lambda_D), ylim=(-container_radius/lambda_D, container_radius/lambda_D))
      circle1 = Circle((0, 0), container_radius/lambda_D,facecolor = "none", edgecolor="black", linewidth=1)
      ax1.add_patch(circle1)

      ax2.set(xlim=(-container_radius/lambda_D, container_radius/lambda_D), ylim=(-container_radius/lambda_D, container_radius/lambda_D))
      circle2 = Circle((0, 0), container_radius/lambda_D,facecolor = "none", edgecolor="black", linewidth=1)
      ax2.add_patch(circle2)

      text = ax1.text(0 + 0.01, 0 + 0.01,'', transform = ax1.transAxes)

   plt.tight_layout()

   #need vector with all the x points and all the y points
   def get_data_dual(v,layer_bottom,layer_top):
      x_lower = []
      y_lower = []
      x_upper = []
      y_upper = []
      for i in np.arange(dust_grain_max):
         if(data["Z_" + str(i)].iloc[v] < layer_bottom):
            x_lower.append(data["X_" + str(i)].iloc[v])
            y_lower.append(data["Y_" + str(i)].iloc[v])
         else:
            x_upper.append(data["X_" + str(i)].iloc[v])
            y_upper.append(data["Y_" + str(i)].iloc[v])

      time = data["Time_list"].iloc[v]
      return [[x_lower,y_lower],[x_upper,y_upper],time]


   def init_dual():
      """initialize animation"""
      dust_points_lower[0].set_data([], [])
      dust_points_upper[0].set_data([], [])
      text.set_text("Time = " + str(0) + "s")
      
      return [dust_points_lower[0],dust_points_upper[0],text,time_line]

   def animate_dual(i,layer_bottom,layer_top):
      """perform animation step"""
      frame_num = i*speed_mul

      data_anim = get_data_dual(frame_num,layer_bottom,layer_top)
      
      # update pieces of the animation
      dust_points_lower[0].set_data(data_anim[0][0],data_anim[0][1])
      dust_points_upper[0].set_data(data_anim[1][0],data_anim[1][1])
      text.set_text("Time = " + str(round(data_anim[2],5)) + "s")#updat value of the frame number
      time_line.set_data([data_anim[2], data_anim[2]], [0, 1])

      return [dust_points_lower[0],dust_points_upper[0],text,time_line]

   #frames = round(frames/speed_mul) 
   #YOU HAVE TO MANUALLY SET THE LAYERS
   ani = animation.FuncAnimation(fig, animate_dual, interval=10, blit=True, frames = frames_anim, init_func=init_dual,fargs = [layers[1],layers[2]])
   #FFwriter=animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
   #ani.save('fek.mp4', writer=FFwriter)

if anim_dual_overlay == "yes":
   size_mul = 1.3
   f_p_s = 30
   anim_time_lim = 30
   speed_mul = round(frames_len/(f_p_s*anim_time_lim))
   frames_anim = int(round(frames_len/speed_mul))

   fig = plt.figure(figsize =[6,6])#figsize = [8,5]
   ax1 = fig.add_subplot()
   # ax3 = fig.add_subplot(2,1,2)

   # for i in np.arange(dust_grain_max):
   #    last_val_index = np.where(data["Time_list"] == last_time_val)
   #    ax3.plot(data["Time_list"], data["Z_" + str(i)])
   #    ax3.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
   #    ax3.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')

   # ax3.plot(data["Time_list"],y_z_se, "--", color = "black")
   # ax3.set_xlabel(r"Time")
   # ax3.set_ylabel(r"z/$\lambda_D$")
   # ax3.set_ylim(0,z_se*top_container_graph_mul/lambda_D)
   # ax3.grid()
   # time_line = ax3.axvline(x=0.0,linewidth=2, color='r')

   ax1.set_xlabel(r"x/$\lambda_D$")
   ax1.set_ylabel(r"y/$\lambda_D$")
   ax1.set_title(r"Crystal Layers Animation")
   dust_points_lower = ax1.plot([], [],"o", color = "green",label = "Lower Layer")
   dust_points_upper = ax1.plot([], [],"o", color = "blue",label = "Upper Layer")
   ax1.legend(fontsize =  "small")

   if boundry == "Periodic":
      ax1.set(xlim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul), ylim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul))
      rect1 = Rectangle((-container_length/(2*lambda_D), -container_length/(2*lambda_D)),width =  container_length/lambda_D, height = container_length/lambda_D, facecolor = "none", edgecolor="black", linewidth=1)
      ax1.add_patch(rect1)
      
      text = ax1.text(0 + 0.01, 0 + 0.01,'', transform = ax1.transAxes)

   else:
      ax1.set(xlim=(-container_radius/lambda_D, container_radius/lambda_D), ylim=(-container_radius/lambda_D, container_radius/lambda_D))
      circle1 = Circle((0, 0), container_radius/lambda_D,facecolor = "none", edgecolor="black", linewidth=1)
      ax1.add_patch(circle1)

      text = ax1.text(0 + 0.01, 0 + 0.01,'', transform = ax1.transAxes)

   plt.tight_layout()

   #need vector with all the x points and all the y points
   def get_data_dual_overlay(v,layer_bottom,layer_top):
      x_lower = []
      y_lower = []
      x_upper = []
      y_upper = []
      for i in np.arange(dust_grain_max):
         if(data["Z_" + str(i)].iloc[v] < layer_bottom):
            x_lower.append(data["X_" + str(i)].iloc[v])
            y_lower.append(data["Y_" + str(i)].iloc[v])
         else:
            x_upper.append(data["X_" + str(i)].iloc[v])
            y_upper.append(data["Y_" + str(i)].iloc[v])

      time = data["Time_list"].iloc[v]
      return [[x_lower,y_lower],[x_upper,y_upper],time]


   def init_dual_overlay():
      """initialize animation"""
      dust_points_lower[0].set_data([], [])
      dust_points_upper[0].set_data([], [])
      text.set_text("Time = " + str(0) + "s")
      
      return [dust_points_lower[0],dust_points_upper[0],text]

   def animate_dual_overlay(i,layer_bottom,layer_top):
      """perform animation step"""
      frame_num = i*speed_mul

      data_anim = get_data_dual_overlay(frame_num,layer_bottom,layer_top)
      
      # update pieces of the animation
      dust_points_lower[0].set_data(data_anim[0][0],data_anim[0][1])
      dust_points_upper[0].set_data(data_anim[1][0],data_anim[1][1])
      text.set_text("Time = " + str(round(data_anim[2],5)) + "s")#updat value of the frame number
      #time_line.set_data([data_anim[2], data_anim[2]], [0, 1])

      return [dust_points_lower[0],dust_points_upper[0],text]

   #frames = round(frames/speed_mul) 
   #YOU HAVE TO MANUALLY SET THE LAYERS
   ani = animation.FuncAnimation(fig, animate_dual_overlay, interval=10, blit=True, frames = frames_anim, init_func=init_dual_overlay,fargs = [layers[1],layers[2]])
   #FFwriter=animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264'])
   #ani.save('Animation_Periodic_100_100000.mp4', writer=FFwriter)

plt.show()