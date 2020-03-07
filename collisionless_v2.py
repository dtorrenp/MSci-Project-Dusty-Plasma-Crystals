import numpy as np
import sys
import matplotlib.pyplot as plt #Graphing module
#import scipy as sp
from scipy.optimize import minimize_scalar
from scipy.special import iv
import matplotlib.animation as animation
np.set_printoptions(threshold=sys.maxsize)

"""
changed ion mass number from 1u to 18u
changed dz from lambda_D/100 to 1/100 as expected (because of normalisation in z_list)
"""

#%%

#SIMULATION PARAMETERS
dust_grain_max_input = 5 #maximum number of dust grains
dt_a = 1.0e-4 #timestep [seconds]
time_limit = 0.5 #maximum runtime [seconds]
frame_req = 5

#VARIABLES
n_e0 = 1.0e15 #electron number density (in bulk plasma)
n_i0 = 1.0e15 #ion number density (in bulk plasma)
n_n0 = 1.0e24 #neutrals number density (in bulk plasma)
Z = 18 #ion atomic number (Ar)
Z_n = 18 #neutrals atomic number (Ar)
grain_R = 7*1e-6 #dust grain radius
dust_grain_density = 2*1e3 #dust density
phi_wall_r = -1000.0 #radial wall potential [volts]
init_speed = 1e-5

#PHYSICAL CONSTANTS
g_z = 9.81 #gravitational field strength
epsilon_0 = 8.85*1e-12 #vacuum permittivity
k_b = 1.38*1e-23 #Boltzmann constant

e_charge = 1.6*1e-19 #electron charge (magnitude)
i_charge = 1.6*1e-19 #ion charge (magnitude)
m_e = 9.11*1e-31 #electron mass
m_i = Z*1.67*1e-27 #ion mass
m_n = m_i #neutrals mass
m_D = ((4.0/3.0)*np.pi*pow(grain_R,3))*dust_grain_density #mass of spherical dust grain
mu = (m_i/m_e) #mass normalisation

T_e = 2.0*(1.6*1e-19)/k_b #electron temperature [K]
T_i = 0.04*(1.6*1e-19)/k_b #ion temperature [K]
beta = T_i/T_e #temperature normalisation

lambda_de = pow(((epsilon_0*k_b*T_e)/(n_e0*(pow(e_charge,2)))),0.5) #electron Debye length
lambda_di = pow(((epsilon_0*k_b*T_i)/(n_i0*(pow(e_charge,2)))),0.5) #ion Debye length
lambda_D = pow((1/(1/(pow(lambda_de,2)) + 1/(pow(lambda_di,2)))),0.5) #total Debye length
wake_potential_below = 1*lambda_D #distance of wake potential beneath dust grain
wake_charge_multiplier = 1.0 #wake charge fraction of dust grain charge
container_height = 100*lambda_D #dust grain drop height
container_height_numerical = 100 #dust grain drop height / Debye lengths
container_radius = 25.0*lambda_D #wall radius

freq = 13.56*1e6 #frequency [Hz] (13.56MHz as quoted from Nitter)
omega = (2*np.pi)*freq #angular frequency [rad/s]
T = 1/freq #period
V_RF = 48.2 #(as quoted from Nitter)
Y_DC_epsilon = 1e-5
Y_epsilon = 0
Q_epsilon = 0
Y_epsilon_bisection = 1e-5
Q_epsilon_bisection = 1e-2
u_epsilon_bisection = 1e-5
sigma = 3e-19
n_s = n_e0*np.exp(-0.5) #electron density at the sheath edge; see Nitter (2)
alpha = 0.134 #((epsilon_0*k_b*T_e)/(n_s*(e_charge**2)) )**0.5*n_n0*sigma
n_d = (dust_grain_max_input)/(np.pi*pow(container_radius,2)*container_height)
#print("alpha = ", alpha)

#%% Producing array coordinates

dT = T/4
dz = 1/100 #in Debye lengths

def produce_z_list():
    """Produces discrete z-coordinates in steps of d_z and units of Debye lengths."""
    z_list = []
    for z in np.arange(0, container_height_numerical, dz):
        z_list.append(z)
    return z_list

def produce_T_list():
    """Produces discrete time coordinates in steps of dT."""
    T_list = []
    for t in np.arange(0,T,dT):
        T_list.append(t)
    return T_list

###
z_list = produce_z_list()
T_list = produce_T_list()

#%%

#STEP 1: Get Y_DC
"""The initial Y_DC sheath potential is used as the basis for both the DC and RF cases. 
This initial potential guess is based on Cameron's collisionless equations."""

V_DC = 0.5*np.log(2*np.pi*m_e/m_i) - np.log(iv(0,V_RF)) #Produces the self-bias V_DC in both the RF case and DC case (where V_RF = 0). V in normalised units of eV/kT.

def f_Y_DC(Y_z_0):
    """RK4 derivative function"""
    return (2*(np.exp(Y_z_0) + (1-2*Y_z_0)**0.5 - 2))**0.5 #Cameron (2.36)

def step_Y_DC(Y_z_0):
    """RK4 step of size dz for Y_DC"""
    k1 = dz * f_Y_DC(Y_z_0)
    k2 = dz * f_Y_DC(Y_z_0 + k1/2)
    k3 = dz * f_Y_DC(Y_z_0 + k2/2)
    k4 = dz * f_Y_DC(Y_z_0 + k3)
    Y_z_1 = Y_z_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
    return Y_z_1

###
Y_DC=[V_DC] #Producing Y_z_init with boundary condition Y_wall = V_DC
Q_DC =[f_Y_DC(V_DC)] #Producing  Q_init i.e. dY_DC/dz
for j in np.arange(len(z_list) - 1):
    Y_DC.append(step_Y_DC(Y_DC[j]))
    Q_DC.append(f_Y_DC(Y_DC[j]))

"""
plt.figure()
plt.plot(z_list, Y_DC)
plt.show()
"""

#%%%

#STEP 2: Solve for next iteration of Y; 1st BC is Y_z(0)=V_DC; 2nd BC of Q at wall must be guessed and iterated with shooting method

def produce_YQ_z_t_wall():
    Y_z_t = []
    Q = []
    for i in np.arange(len(T_list)):
        Y_z_t[i] = [V_DC +  V_RF*np.sin(omega*T_list[i])]
        Q[i] = [f_Y_DC(V_DC +  V_RF*np.sin(omega*T_list[i])]

def f_Y(Q):
    return Q

def g_Q(Y, Y_DC):
    return (np.exp(Y) - (1-2*Y_DC)**(-0.5))

def step_YQ(Y_0, Q_0, Y_DC_0, Y_DC_1):
    Y_DC_half = (Y_DC_0 + Y_DC_1)/2
    k1 = dz * f_Y(Q_0)
    l1 = dz * g_Q(Y_0, Y_DC_0)
    k2 = dz * f_Y(Q_0 + l1/2)
    l2 = dz * g_Q(Y_0 + k1/2, Y_DC_half)
    k3 = dz * f_Y(Q_0 + l2/2)
    l3 = dz * g_Q(Y_0 + k2/2, Y_DC_half)
    k4 = dz * f_Y(Q_0 + l3)
    l4 = dz * g_Q(Y_0 + k3, Y_DC_1)
    Y_1 = Y_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
    Q_1 = Q_0 + (l1 + 2*l2 + 2*l3 + l4)/6.0
    return np.asarray([Y_1, Q_1])

def produce_YQ(Y_init, Q_init, Y_DC):
    Y = [Y_init]
    Q = [Q_init]
    for j in np.arange(len(z_list)-1):
        QY = step_YQ(Y[j], Q[j], Y_DC[j], Y_DC[j+1])
        Y.append(QY[0])
        Q.append(QY[1])
    return [Y, Q]

YQ_z_t_wall = produce_YQ_z_t_wall()

for i in np.arange(len(T_list)):
    Y_init = YQ_z_t_wall[i][0]
    Q_init = YQ_z_t_wall[i][1]
    Y = []
    Q = []
    Y.insert(i, produce_YQ(Y_init, Q_init, Y_DC)[0])
    Q.insert(i, produce_YQ(Y_init, Q_init, Y_DC)[0])

plt.figure()
plt.plot(z_list, Y)
plt.xlim(0,100)
plt.show()

#STEP 2b: We must internally iterate by shooting different Q_init

"""
def bisection(Y_init, Q_init, u):
    YQ = produce_YQ(Y_init, Q_init, u)
    Y = YQ[0]
    Q = YQ[1]
    Q_1 = Q[0]
    while (Q[-1] > Q_epsilon_bisection):
        Q_0 = Q_1
        print (Q_0)
        if Y[-1] >= 0:
            Q_1 = 0.1*Q_0
        else:
            Q_1 = 1.1*Q_0
        YQ = produce_YQ(Y_init, Q_1, u)
        Y = YQ[0]
        Q = YQ[1]
    return (Y, Q)
"""