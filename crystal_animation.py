import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
from csv import reader
import subprocess
import time
import pandas as pd
import matplotlib.animation as animation
import scipy.stats as stats
from matplotlib.patches import Circle
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

FILENAME = "HPC_Data/Dust_grain_max_20_wake_charge_multiplier_0.500000_container_radius_0.001011_Final_Termperature_3933.559427_frames_1001.csv"

#%%
data = pd.read_csv(FILENAME)
dust_grain_max = int((len(data.columns))/4)

fig, ax = plt.subplots()
ax.set(xlim=(-container_radius/lambda_D, container_radius/lambda_D), ylim=(-container_radius/lambda_D, container_radius/lambda_D))

circle = Circle((0, 0), container_radius/lambda_D,facecolor = "none", edgecolor="black", linewidth=3, alpha=0.5)
ax.add_patch(circle)
dust_points = ax.plot([], [],"bo", mec = "black",ms=grain_R/lambda_D)

#need vector with all the x points and all the y points
def get_data(v):
    x = []
    y = []
    for i in np.arange(dust_grain_max):
        x.append(data["X_" + str(i)].iloc[v])
        y.append(data["Y_" + str(i)].iloc[v])
    return [x,y]

def init():
    """initialize animation"""
    dust_points[0].set_data([], [])
    return dust_points

def animate(v):
    """perform animation step"""
    data_anim = get_data(v)
    # update pieces of the animation
    dust_points[0].set_data(data_anim[0],data_anim[1])
    return dust_points

ani = animation.FuncAnimation(fig, animate,interval=10, blit=True, init_func=init)

plt.show()
                                                                   


