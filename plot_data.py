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

status = "Run"#input("Compile or Run?")
load = "no"
boundry = "Periodic"#input("Periodic or Finite?")
standard_plots = "no"
layer_plots = "yes"#input("layers?")

x_plot = "yes"
y_plot = "yes"
z_plot = "yes"
final_pos_plot = "yes"
motion_plot = "yes" 
temp_plot = "yes"
speed_plot = "yes"

voronoi_plot = "yes"
delaunay_plot = "yes"
g_plot = "yes"
layer_final_pos_plot = "yes"
threeD_plot = "no"#doesnt work atm
anim = "no"#doesnt work atm

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
DATA_NAME = input("Data name?(NOT FILENAME)")#"Finite_Dust_grain_max_5_Final_Temperature_48.356239_frames_501"
FILENAME = input("FILE name?")#"HPC_Data/Finite_Dust_grain_max_5_Final_Temperature_48.356239_frames_501.csv"
NEW_FOLDER = "Figures/" + DATA_NAME
if str(os.path.exists(NEW_FOLDER)) == "False":
   os.mkdir(NEW_FOLDER)

#%%
data = pd.read_csv(FILENAME)
dust_grain_max = int( (len(data.columns))/3 - 1)
last_time_val = data["Time_list"].iloc[-1]
y_z_se = [z_se/lambda_D]*len(data["Time_list"])
temp_ion = [T_i]*len(data["Time_list"])

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
         plt.xlabel("Time")
         plt.ylabel("x/lambda_D")
         plt.grid()
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.savefig(NEW_FOLDER + "/Periodic_x_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if y_plot == "yes":
         plt.figure()
         plt.title("test - y")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Y_" + str(i)])
            plt.plot(data["Time_list"][0], data["Y_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("Time")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         plt.savefig(NEW_FOLDER + "/Periodic_y_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if z_plot == "yes":
         plt.figure()
         plt.title("test - z")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Z_" + str(i)])
            plt.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("Time")
         plt.ylabel("z/lambda_D")
         plt.plot(data["Time_list"],y_z_se, "--", color = "black")
         plt.ylim(0,container_height/lambda_D)
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

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
         plt.savefig(NEW_FOLDER + "/Periodic_final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

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

      if temp_plot == "yes":
         plt.figure()
         plt.title("Temperature")
         plt.plot(data["Time_list"],data["Temperature_list"])
         plt.plot(data["Time_list"],temp_ion, "--", color = "black")
         plt.yscale("log")
         plt.xlabel("Time")
         plt.ylabel("Temperature(K)")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if speed_plot == "yes":
         plt.figure()
         plt.title("Speed")
         plt.hist(data["Speed_list"].dropna(), bins = int(dust_grain_max/hist_div_val)  )
         plt.ylabel("Speed")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Periodic_Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

   else:
      theta = np.linspace(0, 2*np.pi, 100)
      x_r_se = r_se_inv/lambda_D*np.cos(theta)
      y_r_se = r_se_inv/lambda_D*np.sin(theta)
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
         plt.xlabel("Time")
         plt.ylabel("x/lambda_D")
         plt.grid()
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_x_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if y_plot == "yes":
         plt.figure()
         plt.title("test - y")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Y_" + str(i)])
            plt.plot(data["Time_list"][0], data["Y_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("Time")
         plt.ylabel("y/lambda_D")
         plt.grid()
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_y_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if z_plot == "yes":
         plt.figure()
         plt.title("test - z")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["Time_list"], data["Z_" + str(i)])
            plt.plot(data["Time_list"][0], data["Z_" + str(i)][0],"+" ,color='blue')
            plt.plot(data["Time_list"].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("Time")
         plt.ylabel("z/lambda_D")
         plt.plot(data["Time_list"],y_z_se, "--", color = "black")
         plt.ylim(0,container_height/lambda_D)
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Finite_z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if final_pos_plot == "yes":
         plt.figure()
         plt.title("Final Positions")
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
         plt.xlabel("x/lambda_D")
         plt.ylabel("y/lambda_D")
         plt.plot(x_r_se,y_r_se, "--", color = "black")
         plt.plot(x_r_wall,y_r_wall, color = "black")
         plt.grid()
         plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
         plt.savefig(NEW_FOLDER + "/Finite_final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if temp_plot == "yes":
         plt.figure()
         plt.title("Temperature")
         plt.plot(data["Time_list"],data["Temperature_list"])
         plt.plot(data["Time_list"],temp_ion, "--", color = "black")
         plt.yscale("log")
         plt.xlabel("Time")
         plt.ylabel("Temperature(K)")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Finite_Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if speed_plot == "yes":
         plt.figure()
         plt.title("Speed")
         plt.hist(data["Speed_list"].dropna(), bins = int(dust_grain_max/hist_div_val)  )
         plt.ylabel("Speed")
         plt.grid()
         plt.savefig(NEW_FOLDER + "/Finite_Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

if layer_plots == "yes":
   layers = [0.0,9.9,11]
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
         #print("hi")
         #print(layer_row)
         points = np.c_[layer_row[0], layer_row[1]]

         if voronoi_plot == "yes":
            vor = Voronoi(points)
            #fig = voronoi_plot_2d(vor)
            fig_v = plt.figure()
            ax_v = fig_v.add_subplot(111)
            plt.title("Voronoi layer: " + str(v))
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            voronoi_plot_2d(vor, ax = ax_v)
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Voronoi" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

         if delaunay_plot == "yes":
            tri = Delaunay(points)
            plt.figure()
            plt.title("Delaunay_Triangulation layer: " + str(v))
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.triplot(points[:,0], points[:,1], tri.simplices)
            plt.plot(points[:,0], points[:,1], 'o')
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Delaunay_Triangulation" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

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
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "Pair_correlation" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

         if layer_final_pos_plot == "yes":
            plt.figure()
            plt.title("Final Positions layer: " + str(v))
            plt.plot(layer_row[0],layer_row[1],"+" ,color='red')
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.grid()
            plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
            plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
            plt.savefig(NEW_FOLDER + "/Periodic_layers" + str(v) + "final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")
      
      if threeD_plot == "yes":
         fig = plt.figure()
         ax = fig.add_subplot(111, projection='3d')
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            for k in np.arange(len(layers) - 1):
               if ( (data["Z_" + str(i)].values[last_val_index][0] >= layers[k] ) and (data["Z_" + str(i)].values[last_val_index][0] <= layers[k+1])):
                  colour_layer = colour_list[v]
                  ax.scatter(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],data["Z_" + str(i)].iloc[last_val_index], color = colour_layer)
         ax.set_xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         ax.set_ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
         ax.set_zlim(0,top_container_graph_mul*container_height/lambda_D)
         ax.set_xlabel('X position')
         ax.set_ylabel('Y position')
         ax.set_zlabel('Z position')

   else:
      theta = np.linspace(0, 2*np.pi, 100)
      x_r_se = r_se_inv/lambda_D*np.cos(theta)
      y_r_se = r_se_inv/lambda_D*np.sin(theta)
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
            plt.title("Voronoi layer: " + str(v))
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            voronoi_plot_2d(vor, ax = ax_v)
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "Voronoi" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

         if delaunay_plot == "yes":
            tri = Delaunay(points)
            plt.figure()
            plt.title("Delaunay_Triangulation layer: " + str(v))
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.triplot(points[:,0], points[:,1], tri.simplices)
            plt.plot(points[:,0], points[:,1], 'o')
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "Delaunay_Triangulation" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")
            
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
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "Pair_correlation" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

         if layer_final_pos_plot == "yes":
            plt.figure()
            plt.title("Final Positions layer: " + str(v))
            plt.plot(layer_row[0],layer_row[1],"+" ,color='red')
            plt.xlabel("x/lambda_D")
            plt.ylabel("y/lambda_D")
            plt.grid()
            plt.xlim(-container_radius/(lambda_D),container_radius/(lambda_D))
            plt.ylim(-container_radius/(lambda_D),container_radius/(lambda_D))
            plt.savefig(NEW_FOLDER + "/Finite_layers" + str(v) + "final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list"])) + ".png")

      if threeD_plot == "yes":
         fig = plt.figure()
         ax = fig.add_subplot(111, projection='3d')
         for i in np.arange(dust_grain_max):
            last_val_index = np.where(data["Time_list"] == last_time_val)
            for k in np.arange(len(layers) - 1):
               if ( (data["Z_" + str(i)].values[last_val_index][0] >= layers[k] ) and (data["Z_" + str(i)].values[last_val_index][0] <= layers[k+1])):
                  colour_layer = colour_list[k]
                  ax.scatter(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],data["Z_" + str(i)].iloc[last_val_index], color = colour_layer)
         ax.set_xlim(-container_radius/lambda_D,container_radius/lambda_D)
         ax.set_ylim(-container_radius/lambda_D,container_radius/lambda_D)
         ax.set_zlim(0,top_container_graph_mul*container_height/lambda_D)
         ax.set_xlabel('X position')
         ax.set_ylabel('Y position')
         ax.set_zlabel('Z position')

   if anim == "yes":
      for n in np.arange(len(layers) - 1):
         frames = len(data["Time_list"])
         speed_mul = 100
         size_mul = 1.2

         fig, ax = plt.subplots()
         plt.xlabel("x/lambda_D")
         plt.ylabel("y/lambda_D")
         plt.title("Animation layer: " + str(v))
         dust_points = ax.plot([], [],"o")

         if boundry == "Periodic":
            ax.set(xlim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul), ylim=(-container_length/(2*lambda_D)*size_mul, container_length/(2*lambda_D)*size_mul))
            rect = Rectangle((-container_length/(2*lambda_D), -container_length/(2*lambda_D)),width =  container_length/lambda_D, height = container_length/lambda_D, facecolor = "none", edgecolor="black", linewidth=1)
            ax.add_patch(rect)
            text = plt.text(-container_length/(2*lambda_D)*size_mul + 0.2, container_length/(2*lambda_D)*size_mul - 0.5,"")
         else:
            ax.set(xlim=(-container_radius/lambda_D, container_radius/lambda_D), ylim=(-container_radius/lambda_D, container_radius/lambda_D))
            circle = Circle((0, 0), container_radius/lambda_D,facecolor = "none", edgecolor="black", linewidth=1)
            ax.add_patch(circle)
            text = plt.text(-container_radius/lambda_D + 0.5, container_radius/lambda_D - 2.5,"")

         #need vector with all the x points and all the y points
         def get_data(v,layer_bottom,layer_top):
            x = []
            y = []
            for i in np.arange(dust_grain_max):
               if ( (data["Z_" + str(i)].values[last_val_index][0] >= layer_bottom ) and (data["Z_" + str(i)].values[last_val_index][0] <= layer_top)):
                  x.append(data["X_" + str(i)].iloc[v])
                  y.append(data["Y_" + str(i)].iloc[v])
            time = data["Time_list"].iloc[v]
            return [x,y,time]

         def init():
            """initialize animation"""
            dust_points[0].set_data([], [])
            text.set_text("Time = " + str(0) + "s")
            return [dust_points[0],text]

         def animate(i,layer_bottom,layer_top):
            """perform animation step"""
            data_anim = get_data(speed_mul*i,layer_bottom,layer_top)
            # update pieces of the animation
            dust_points[0].set_data(data_anim[0],data_anim[1])
            text.set_text("Time = " + str(round(data_anim[2],5)) + "s")#updat value of the frame number
            return [dust_points[0],text]
         
         ani = animation.FuncAnimation(fig, animate,interval=10, blit=True, init_func=init, frames = round(frames/speed_mul), fargs = [layers[n],layers[n+1]])
         #ani.save("Animation layer: " + str(v),fps=30) #save command
         plt.close()
plt.show()