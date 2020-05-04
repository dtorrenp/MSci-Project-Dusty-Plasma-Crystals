import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
np.set_printoptions(threshold=sys.maxsize)

#%%
#CRITICAL VALUES
dust_grain_max_input = 5 #dust grain max number
dt_a = 1.0e-4
time_limit = 0.5
frame_req = 5

#CONSTANTS TO FUCK ABOUT WITH
n_e0 = 1.0e15 #electron number density (in bulk plasma)
n_i0 = 1.0e15 #ion number density (in bulk plasma)
n_n0 = 1.0e24 #electron number density (in bulk plasma)
Z = 1 #ion atomic number
Z_n = 18 #neutrals atomic number (Ar)
grain_R = 7*1e-6 #dust grain radius
dust_grain_density = 2*1e3 #dust density
phi_wall_r = -1000.0 #radial wall potential [volts]
init_speed = 1e-5

#CONSTANTS DEPENDANT ON ACTUAL PHYSICS
g_z = 9.81#gravity
e_charge = 1.6*1e-19#magnitude of e charge
i_charge = 1.6*1e-19#magnitude of i charge DOES THIS NEED TO CHANGE WHEN USED IN THE ION DRAG?????
m_i = 1.67*1e-27
m_e = 9.11*1e-31
m_n = Z_n*m_i
m_D = ((4.0/3.0)*np.pi*pow(grain_R,3))*dust_grain_density#m_D of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
epsilon_0 = 8.85*1e-12
k_b = 1.38*1e-23
mu = (m_i/m_e)#normalization used for convienience
T_e = 2.0*(1.6*1e-19)/k_b
T_i = 0.03*(1.6*1e-19)/k_b
beta = T_i/T_e
lambda_de = pow(((epsilon_0*k_b*T_e)/(n_e0*(pow(e_charge,2)))),0.5)
lambda_di = pow(((epsilon_0*k_b*T_i)/(n_i0*(pow(e_charge,2)))),0.5)
lambda_D = pow((1/(1/(pow(lambda_de,2)) + 1/(pow(lambda_di,2)))),0.5)#get lambda_D

#SET THESE?
drop_height = 10*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 25.0*lambda_D

dz = lambda_D/2
z_0 = 0
z_se = -10*lambda_D
p_0 = 5*lambda_D
Q = -1.8*1e-15
q_d = -1.8*1e-15

t = 0
z = np.linspace(0,z_se,100)
p = np.linspace(-p_0,p_0,100)

#THIS GIVES THE HIGHT DEPENDANT BUT WE ONLY WANT 1 M PER ENTIRE GRAPH
# phi_wall_z = -100
# k_z_restore = -2.0*phi_wall_z/pow(z_se,2)
# C_s = ((3*k_b*T_e)/m_i)**0.5
# v_i_z = -(C_s**2 + (i_charge*k_z_restore*((z - z_se)**2))/m_i - 2.0*g_z*(z - z_se))**0.5
v = 0
# v_0 = -1.1*C_s
# v_d = v - v_i_z
# M = v_d/C_s

#print(grain_R/lambda_D)
#exit(
M = 1.2

near_plot = "no"
far_plot = "yes"
electric_far_plot  = "yes"
Force_far_plot = "yes"

#%%
if near_plot == "yes":
    #NEAR FIELD PLOT
    phi_near = (1/(1 - M_single**(-2)))  *  ((2*Q)/np.abs(z - z_0 - v*t))  *  np.cos((z - z_0 - v*t)/(lambda_D*(M_single**2 - 1)**0.5))
    plt.figure()
    plt.plot(z/lambda_D,phi_near)
    plt.xlabel(r"z/$\lambda_D$")
    plt.ylabel(r"$\phi$")
    plt.grid()
    plt.savefig("Figures/near_field_wake_Mach_" + str(M_single)+".png")

if far_plot == "yes":
    #FAR FIELD PLOT
    phi_far = np.zeros((len(p),len(z)))
    for i in np.arange(len(phi_far)):
        z_plus = np.abs(z - z_0 - v*t) + np.abs(p[i])*(M**2-1)**0.5
        z_minus = np.abs(z - z_0 - v*t) - np.abs(p[i])*(M**2-1)**0.5
        phi_far[i] = (2*Q/(1 - M**(-2)))*((lambda_D/(2*np.pi*np.abs(p[i])))**0.5)*((1/z_plus)*(np.cos(((z_plus/lambda_D)/(M**2 -1)**0.5) -np.pi/4) - 1/(2**0.5)) +  (1/z_minus)*(np.cos(((z_minus/lambda_D)/(M**2-1)**0.5) + np.pi/4) - 1/(2**0.5)))
    phi_far = phi_far.transpose()

    fig, ax = plt.subplots(1, 1)
    #print(np.shape(p/lambda_D),np.shape(z/lambda_D), np.shape(phi_far))
    cp = ax.contourf(p/lambda_D,z/lambda_D,phi_far)
    cb = fig.colorbar(cp)
    circle = Circle((0, 0), radius = grain_R/lambda_D,color='black')
    ax.add_patch(circle)
    cb.set_label(r"$\rho$")
    ax.set_xlabel(r"$\rho$/$\lambda_D$")
    ax.set_ylabel(r"z/$\lambda_D$")
    plt.savefig("Figures/far_field_wake.png")

if electric_far_plot == "yes":
    #FAR FIELD E PLOT
    E_p_far = np.zeros((len(p),len(z)))
    E_z_far = np.zeros((len(p),len(z)))
    E_mag = np.zeros((len(p),len(z)))
    for i in np.arange(len(E_p_far)):
        z_plus = np.abs(z - z_0 - v*t) + np.abs(p[i])*(M**2-1)**0.5
        z_minus = np.abs(z - z_0 - v*t) - np.abs(p[i])*(M**2-1)**0.5
        A = np.cos((z_plus/lambda_D)/(M**2 - 1)**0.5 - np.pi/4)
        B = np.cos((z_minus/lambda_D)/(M**2 - 1)**0.5 + np.pi/4)
        C = np.sin((z_plus/lambda_D)/(M**2 - 1)**0.5 - np.pi/4)
        D = np.sin((z_minus/lambda_D)/(M**2 - 1)**0.5 + np.pi/4)

        E_p_far[i] = -1*(p[i]/np.abs(p[i]))*(2*Q/(1 - M**(-2)))*((lambda_D/(2*np.pi*np.abs(p[i])))**0.5)*((1/np.abs(p[i]))*(-0.5)*((1/z_plus)*(A - 1/(2**0.5)) + (1/z_minus)*(B - 1/(2**0.5))) + ((-z_plus**(-2))*((M**2 - 1)**0.5)*(A - 1/(2**0.5)) - (z_plus**-1)*(C/lambda_D) + (z_minus**(-2))*((M**2 - 1)**0.5)*(B - 1/(2**0.5)) + (z_minus**-1)*(D/lambda_D)))
        E_z_far[i] = -1*(2*Q/(1 - M**(-2)))*((lambda_D/(2*np.pi*np.abs(p[i])))**0.5)*( -(z_plus**(-2))*(A - 1/(2**0.5)) - (z_plus**(-1))*(C*(1/lambda_D)/((M**2 - 1)**0.5)) - (z_minus**(-2))*(B - 1/(2**0.5)) - (z_minus**(-1))*(D*(1/lambda_D)/((M**2 - 1)**0.5)))
        E_mag[i] = ((E_p_far[i])**2 + (E_z_far[i])**2)**0.5

    #print(np.shape(p/lambda_D),np.shape( z/lambda_D), np.shape(E_p_far), np.shape(E_z_far), np.shape(E_mag))
    E_p_far = E_p_far.transpose()
    E_z_far = E_z_far.transpose()
    E_mag = E_mag.transpose()
    #print(np.shape(p/lambda_D),np.shape( z/lambda_D), np.shape(E_p_far), np.shape(E_z_far), np.shape(E_mag))
    fig, ax = plt.subplots()
    q = ax.quiver(p/lambda_D, z/lambda_D, E_p_far, E_z_far, E_mag)
    circle = Circle((0, 0), radius = grain_R/lambda_D,color='black')
    ax.add_patch(circle)
    ax.set_xlabel(r"$\rho$/$\lambda_D$")
    ax.set_ylabel(r"z/$\lambda_D$")
    plt.savefig("Figures/Electric_field_far_field_wake.png")

if Force_far_plot == "yes":
    #FAR FIELD E PLOT
    F_p_far = np.zeros((len(p),len(z)))
    F_z_far = np.zeros((len(p),len(z)))
    F_mag = np.zeros((len(p),len(z)))
    for i in np.arange(len(p)):
        z_plus = np.abs(z - z_0 - v*t) + np.abs(p[i])*(M**2-1)**0.5
        z_minus = np.abs(z - z_0 - v*t) - np.abs(p[i])*(M**2-1)**0.5
        A = np.cos((z_plus/lambda_D)/(M**2 - 1)**0.5 - np.pi/4)
        B = np.cos((z_minus/lambda_D)/(M**2 - 1)**0.5 + np.pi/4)
        C = np.sin((z_plus/lambda_D)/(M**2 - 1)**0.5 - np.pi/4)
        D = np.sin((z_minus/lambda_D)/(M**2 - 1)**0.5 + np.pi/4)
        F_p_far[i] = -1*(p[i]/np.abs(p[i])) *q_d*(2*Q/(1 - M**(-2))) * ((lambda_D/(2*np.pi*np.abs(p[i])))**0.5)  *((1/np.abs(p[i]))*(-0.5)*((1/z_plus)*(A - 1/(2**0.5)) + (1/z_minus)*(B - 1/(2**0.5))) + ((-z_plus**(-2))*((M**2 - 1)**0.5)*(A - 1/(2**0.5)) - (z_plus**-1)*(C/lambda_D) + (z_minus**(-2))*((M**2 - 1)**0.5)*(B - 1/(2**0.5)) + (z_minus**-1)*(D/lambda_D)))
        F_z_far[i] = -1*q_d*(2*Q/(1 - M**(-2)))*((lambda_D/(2*np.pi*np.abs(p[i])))**0.5)*(-(z_plus**(-2))*(A - 1/(2**0.5)) - (z_plus**(-1))*(C*(1/lambda_D)/((M**2 - 1)**0.5)) - (z_minus**(-2))*(B - 1/(2**0.5)) - (z_minus**(-1))*(D*(1/lambda_D)/((M**2 - 1)**0.5)))
        F_mag[i] = ((F_p_far[i])**2 + (F_z_far[i])**2)**0.5

    F_p_far = F_p_far.transpose()
    F_z_far = F_z_far.transpose()
    F_mag = F_mag.transpose()

    fig, ax = plt.subplots()
    q = ax.quiver(p/lambda_D, z/lambda_D, F_p_far, F_z_far, F_mag)
    circle = Circle((0, 0), radius = grain_R/lambda_D,color='black')
    ax.add_patch(circle)
    ax.set_title(r"Force - Far field")
    ax.set_xlabel(r"$\rho$/$\lambda_D$")
    ax.set_ylabel(r"z/$\lambda_D$")
    plt.savefig("Figures/Force_far_field_wake.png")

plt.show()


