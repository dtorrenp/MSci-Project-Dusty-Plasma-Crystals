import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
from csv import reader
import subprocess
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
container_radius = 100*lambda_D#set radius of contianer ie wall radius
container_length = 10*lambda_D
z_se = 10*lambda_D#distance from bottom of container to the sheath edge
r_se = 100*lambda_D#distance from wall to the sheathe edge
r_se_inv = container_radius - r_se

#%%

status = input("Compile or Run?")
boundry = input("Periodic or Finite?")
if status == "Compile":
    #COMPILE
    #UNI
    #subprocess.call(["g++", "H:\year 4\computational\MSci-Project-Dusty-Plasma-Crystals\MSci_project.cpp"])
    #LAPTOP
    #os.environ["PROJECT_FILE"]
   if boundry == "Periodic":
      subprocess.call(["g++", "-o", "MSci_project_periodic", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/MSci_project_periodic.cpp"])
      """Import stuff to measure how long the code takes to run"""
      start_time = time.time()
      print("start_time =", time.ctime(time.time()))
      subprocess.call("MSci_project_periodic.exe")
      """prints time taken in minutes"""
      print ("time taken: %s minutes" % ((time.time()-start_time)/60))
   elif boundry == "Finite":
      subprocess.call(["g++", "-o", "MSci_project", "C:/Users/daniel/Documents/UniWork/4th_Year/MSci-Project-Dusty-Plasma-Crystals/MSci_project.cpp"])
      """Import stuff to measure how long the code takes to run"""
      start_time = time.time()
      print("start_time =", time.ctime(time.time()))
      subprocess.call("MSci_project.exe")
      """prints time taken in minutes"""
      print ("time taken: %s minutes" % ((time.time()-start_time)/60))
   else:
      print ("run again")


#%%

FILENAME = input("Data file name?")

#%%

import pandas as pd
data = pd.read_csv(FILENAME)
dust_grain_max = int((len(data.columns))/4)
last_time_val = data["Time_list_0"].iloc[-1]

#%%
if boundry == "Periodic":
   y_z_se = [z_se/lambda_D]*len(data["Time_list_0"])

   plt.figure(1)
   plt.title("Motion")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)].values == last_time_val)
      plt.plot(data["X_" + str(i)],data["Y_" + str(i)])
      plt.plot(data["X_" + str(i)][0], data["Y_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("x/lambda_D")
   plt.ylabel("y/lambda_D")
   plt.grid()
   plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   plt.savefig("Figures/Periodic_Path_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(2)
   plt.title("test - x")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["Time_list_" +str(i)], data["X_" + str(i)])
      plt.plot(data["Time_list_" +str(i)][0], data["X_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["Time_list_" +str(i)].iloc[last_val_index], data["X_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("Time")
   plt.ylabel("x/lambda_D")
   plt.grid()
   plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   plt.savefig("Figures/Periodic_x_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(3)
   plt.title("test - y")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["Time_list_" +str(i)], data["Y_" + str(i)])
      plt.plot(data["Time_list_" +str(i)][0], data["Y_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["Time_list_" +str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("Time")
   plt.ylabel("y/lambda_D")
   plt.grid()
   plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   plt.savefig("Figures/Periodic_y_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(4)
   plt.title("test - z")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["Time_list_" +str(i)], data["Z_" + str(i)])
      plt.plot(data["Time_list_" +str(i)][0], data["Z_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["Time_list_" +str(i)].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("Time")
   plt.ylabel("z/lambda_D")
   plt.plot(data["Time_list_0"],y_z_se, "--", color = "black")
   plt.ylim(0,container_height/lambda_D)
   plt.grid()
   plt.savefig("Figures/Periodic_z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(5)
   plt.title("Final Positions")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("x/lambda_D")
   plt.ylabel("y/lambda_D")
   plt.grid()
   plt.xlim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   plt.ylim(-container_length/(2*lambda_D),container_length/(2*lambda_D))
   plt.savefig("Figures/Periodic_final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   fig = plt.figure(6)
   plt.title("Final Positions - 3D")
   ax = fig.add_subplot(111, projection='3d')
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      ax.scatter(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],data["Z_" + str(i)].iloc[last_val_index])
   ax.set_xlabel('X position')
   ax.set_ylabel('Y position')
   ax.set_zlabel('Z position') 

   plt.figure(7)
   plt.title("Temperature")
   plt.plot(data["Time_list_0"],data["Temperature_list"])
   plt.xlabel("Time")
   plt.ylabel("Temperature(K)")
   plt.grid()
   plt.savefig("Figures/Periodic_Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(8)
   plt.title("Speed")
   plt.hist(data["Speed_list"].dropna())
   plt.ylabel("Speed")
   plt.grid()
   plt.savefig("Figures/Periodic_Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

else:
   theta = np.linspace(0, 2*np.pi, 100)
   x_r_se = r_se_inv/lambda_D*np.cos(theta)
   y_r_se = r_se_inv/lambda_D*np.sin(theta)
   x_r_wall = container_radius/lambda_D*np.cos(theta)
   y_r_wall = container_radius/lambda_D*np.sin(theta)
   y_z_se = [z_se/lambda_D]*len(data["Time_list_0"])

   plt.figure(1)
   plt.title("Motion")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)].values == last_time_val)
      plt.plot(data["X_" + str(i)],data["Y_" + str(i)])
      plt.plot(data["X_" + str(i)][0], data["Y_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.plot(x_r_se,y_r_se, "--", color = "black")
   plt.plot(x_r_wall,y_r_wall, color = "black")
   plt.xlabel("x/lambda_D")
   plt.ylabel("y/lambda_D")
   plt.grid()
   plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
   plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
   plt.savefig("Figures/Path_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(2)
   plt.title("test - x")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["Time_list_" +str(i)], data["X_" + str(i)])
      plt.plot(data["Time_list_" +str(i)][0], data["X_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["Time_list_" +str(i)].iloc[last_val_index], data["X_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("Time")
   plt.ylabel("x/lambda_D")
   plt.grid()
   plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
   plt.savefig("Figures/x_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(3)
   plt.title("test - y")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["Time_list_" +str(i)], data["Y_" + str(i)])
      plt.plot(data["Time_list_" +str(i)][0], data["Y_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["Time_list_" +str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("Time")
   plt.ylabel("y/lambda_D")
   plt.grid()
   plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
   plt.savefig("Figures/y_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(4)
   plt.title("test - z")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["Time_list_" +str(i)], data["Z_" + str(i)])
      plt.plot(data["Time_list_" +str(i)][0], data["Z_" + str(i)][0],"+" ,color='blue')
      plt.plot(data["Time_list_" +str(i)].iloc[last_val_index], data["Z_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("Time")
   plt.ylabel("z/lambda_D")
   plt.plot(data["Time_list_0"],y_z_se, "--", color = "black")
   plt.ylim(0,container_height/lambda_D)
   plt.grid()
   plt.savefig("Figures/z_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(5)
   plt.title("Final Positions")
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      plt.plot(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],"+" ,color='red')
   plt.xlabel("x/lambda_D")
   plt.ylabel("y/lambda_D")
   plt.plot(x_r_se,y_r_se, "--", color = "black")
   plt.plot(x_r_wall,y_r_wall, color = "black")
   plt.grid()
   plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
   plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)
   plt.savefig("Figures/final_pos_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   fig = plt.figure(6)
   plt.title("Final Positions - 3D")
   ax = fig.add_subplot(111, projection='3d')
   for i in np.arange(dust_grain_max):
      last_val_index = np.where(data["Time_list_" + str(i)] == last_time_val)
      ax.scatter(data["X_" + str(i)].iloc[last_val_index], data["Y_" + str(i)].iloc[last_val_index],data["Z_" + str(i)].iloc[last_val_index])
   ax.set_xlabel('X position')
   ax.set_ylabel('Y position')
   ax.set_zlabel('Z position') 

   plt.figure(7)
   plt.title("Temperature")
   plt.plot(data["Time_list_0"],data["Temperature_list"])
   plt.xlabel("Time")
   plt.ylabel("Temperature(K)")
   plt.grid()
   plt.savefig("Figures/Temperature_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

   plt.figure(8)
   plt.title("Speed")
   plt.hist(data["Speed_list"])
   plt.ylabel("Speed")
   plt.grid()
   plt.savefig("Figures/Speed_dust_grain_max_" + str(dust_grain_max) + "_frames_" + str(len(data["Time_list_0"])) + ".png")

plt.show()