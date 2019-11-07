import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
from csv import reader

n_e0 = 1e15#electron and ion densities in bulk plasma
n_i0 = 1e15

g_z = 9.81#gravity
e_charge = -1.6*1e-19
i_charge = 1.6*1e-19
grain_R = 7*1e-6
m_i = 1.67*1e-27
m_e = 9.11*1e-31
m_D = ((4/3)*np.pi*grain_R**3)*(1.49*1e3)#mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
epsilon_0 = 8.85*1e-12
k_b = 1.38*1e-23
mu = (m_i/m_e)#normalization used for convienience

beta = 10**-1#T_i/T_e
T_e = 1e2
T_i = T_e*beta#defined a bit wierdly but basically just in terms of beta and T_e
v_i = (2*k_b*T_i/m_i)**(1/2)
Z = 1#atomic number??
a_0 = 1#intial guess for halley's method

lambda_de = ((epsilon_0*k_b*T_e)/(n_e0*(e_charge**2)))**0.5
lambda_di = ((epsilon_0*k_b*T_i)/(n_i0*(e_charge**2)))**0.5
lambda_D = (1/(1/(lambda_de**2) + 1/(lambda_di**2)))**0.5#get lambda_D

container_height = 11*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 11*lambda_D#set radius of contianer ie wall radius
z_se = 10*lambda_D#distance from bottom of container to the sheath edge
r_se = 10*lambda_D#distance from wall to the sheathe edge
r_se_inv = container_radius - r_se

phi_wall_z = -100 #volts
phi_wall_r = -1 #volts
k_z_restore = -2*phi_wall_z/z_se**2#WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
k_r_restore = -2*phi_wall_r/r_se**2

v_B = (k_b*T_e/m_i)**0.5

alpha = 1e-9#drag coefficient
time_list = [0]
root = 1e-14#defined preciscion of root finding method used to get dust charge
#S = 0.95

dt = 1e-4#time step in rk4, needs to be small enough to be precise but large enough we can actually move the stuff forward in time
#eps = 1e-7
dust_grain_max = 10#dust grain max number
frames = 1e5#number of frames, time taken is not linear as teh longer u run it the more particles it adds hence increases quadratically
temp_min = 1e-5#minimum temperature for it to stop
for_run = True

#%%

with open('test_csv.csv', 'r') as f:
    data = list(reader(f))
print(data)

""" 
#%%

theta = np.linspace(0, 2*np.pi, 100)
x_r_se = r_se_inv/lambda_D*np.cos(theta)
y_r_se = r_se_inv/lambda_D*np.sin(theta)
x_r_wall = container_radius/lambda_D*np.cos(theta)
y_r_wall = container_radius/lambda_D*np.sin(theta)


y_z_se = [z_se/lambda_D]*len(time_list)

#%%

plt.figure(1)
plt.title("test")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].x_history, Final_conditions_dust_grains[i].y_history)
    plt.plot(Final_conditions_dust_grains[i].x_history[0], Final_conditions_dust_grains[i].y_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].x_history[-1], Final_conditions_dust_grains[i].y_history[-1],"+" ,color='red')
plt.plot(x_r_se,y_r_se, "--", color = "black")
plt.plot(x_r_wall,y_r_wall, color = "black")
plt.xlabel("x/lambda_D")
plt.ylabel("y/lambda_D")
plt.grid()
plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

plt.figure(2)
plt.title("test - x")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].time_list,Final_conditions_dust_grains[i].x_history)
    plt.plot(Final_conditions_dust_grains[i].time_list[0],Final_conditions_dust_grains[i].x_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].time_list[-1],Final_conditions_dust_grains[i].x_history[-1],"+" ,color='red')
plt.xlabel("Time")
plt.ylabel("x/lambda_D")
plt.grid()
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

plt.figure(3)
plt.title("test - y")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].time_list,Final_conditions_dust_grains[i].y_history)
    plt.plot(Final_conditions_dust_grains[i].time_list[0],Final_conditions_dust_grains[i].y_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].time_list[-1],Final_conditions_dust_grains[i].y_history[-1],"+" ,color='red')
plt.xlabel("Time")
plt.ylabel("y/lambda_D")
plt.grid()
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

plt.figure(4)
plt.title("test - z")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].time_list,Final_conditions_dust_grains[i].z_history)
    plt.plot(Final_conditions_dust_grains[i].time_list[0],Final_conditions_dust_grains[i].z_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].time_list[-1],Final_conditions_dust_grains[i].z_history[-1],"+" ,color='red')
plt.xlabel("Time")
plt.ylabel("z/lambda_D")
plt.plot(time_list,y_z_se, "--", color = "black")
#horiz_line_data = np.array([40 for i in xrange(len(xs))])
#xs = np.linspace(1,21,200)
#plt.plot(xs, horiz_line_data, 'r--') 
plt.ylim(0,container_height/lambda_D)
plt.grid()

plt.figure(5)
plt.title("Final Positions")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].x_history[-1], Final_conditions_dust_grains[i].y_history[-1],"+" ,color='red')
plt.xlabel("x/lambda_D")
plt.ylabel("y/lambda_D")
plt.plot(x_r_se,y_r_se, "--", color = "black")
plt.plot(x_r_wall,y_r_wall, color = "black")
plt.grid()
plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

fig = plt.figure(6)
ax = fig.add_subplot(111, projection='3d')

for i in np.arange(len(Final_conditions_dust_grains)):
    ax.scatter(Final_conditions_dust_grains[i].x_history[-1], Final_conditions_dust_grains[i].y_history[-1], Final_conditions_dust_grains[i].z_history[-1])
ax.set_xlabel('X position')
ax.set_ylabel('Y position')
ax.set_zlabel('Z position') """