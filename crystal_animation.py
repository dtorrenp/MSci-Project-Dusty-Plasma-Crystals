import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
from csv import reader
import subprocess
import time
import pandas as pd
import matplotlib.animation as animation
import scipy.stats as stats
from matplotlib.patches import Circle,Rectangle
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

#%%

status = "Run"#input("Compile or Run?")
boundry = "Periodic"#input("Periodic or Finite?")
layer_plots = "no"#input("layers?")

FILENAME = "HPC_Data/Dust_grain_max_500_wake_charge_multiplier_0.500000_container_length_0.000404_Final_Termperature_nan_frames_12323.csv"#input("Data file name?")

#%%
data = pd.read_csv(FILENAME)
dust_grain_max = int((len(data.columns))/4)
frames = len(data["Time_list_0"])
speed_mul = 100
size_mul = 1.2

fig, ax = plt.subplots()
plt.xlabel("x/lambda_D")
plt.ylabel("y/lambda_D")
dust_points = ax.plot([], [],"ro", ms=grain_R/lambda_D)
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
def get_data(v):
    x = []
    y = []
    for i in np.arange(dust_grain_max):
        x.append(data["X_" + str(i)].iloc[v])
        y.append(data["Y_" + str(i)].iloc[v])
    time = data["Time_list_" + str(0)].iloc[v]
    return [x,y,time]

def init():
    """initialize animation"""
    dust_points[0].set_data([], [])
    text.set_text("Time = " + str(0) + "s")
    return [dust_points[0],text]

def animate(i):
    """perform animation step"""
    data_anim = get_data(speed_mul*i)
    # update pieces of the animation
    dust_points[0].set_data(data_anim[0],data_anim[1])
    text.set_text("Time = " + str(round(data_anim[2],5)) + "s")#updat value of the frame number
    #print(dust_points)
    #print(text)
    return [dust_points[0],text]

ani = animation.FuncAnimation(fig, animate,interval=10, blit=True, init_func=init, frames = round(frames/speed_mul))

plt.show()
                                                                   


