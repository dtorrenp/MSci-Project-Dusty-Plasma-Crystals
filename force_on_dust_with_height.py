import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs

#CONSTANTS TO FUCK ABOUT WITH
n_e0 = 1.0e15#electron and ion densities in bulk plasma
n_i0 = 1.0e15
phi_wall_z = -100.0#volts

#CONSTANTS DEPENDANT ON ACTUAL PHYSICS
g_z = 9.81#gravity
e_charge = 1.6*1e-19#magnitude of e charge
i_charge = 1.6*1e-19#magnitude of i charge DOES THIS NEED TO CHANGE WHEN USED IN THE ION DRAG?????
grain_R = 7*1e-6
m_i = 1.67*1e-27
m_e = 9.11*1e-31
m_D = ((4/3)*np.pi*grain_R**3)*(1.49*1e3)#mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
epsilon_0 = 8.85*1e-12
k_b = 1.38*1e-23
T_e = 2.0*(1.6*1e-19)/k_b
T_i = 0.03*(1.6*1e-19)/k_b
beta = T_i/T_e
lambda_de = (((epsilon_0*k_b*T_e)/(n_e0*((e_charge**2))))**0.5)
lambda_di = (((epsilon_0*k_b*T_i)/(n_i0*((e_charge**2))))**0.5)
lambda_D = ((1/(1/((lambda_de**2)) + 1/((lambda_di**2))))**0.5)#get lambda_D

#related to finding the charge of the dust grains
d_z = lambda_D/1000
low_z = 7.5*lambda_D
high_z = 20.0*lambda_D
a_0_OML = -0.5
a_0_Sheath = -2.5#intial guess for halley's method
root = 1.0e-4#preciscion of root finding method used to get dust charge
points = 10000.0

#container_dust_dist_creation = np.sqrt(container_radius/(2*lambda_D))
z_se = 10.0*lambda_D#distance from bottom of container to the sheath edge
k_z_restore = -2.0*phi_wall_z/(z_se**2)#WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
v_B = ((3*k_b*T_e/m_i)**0.5)


def f_x_OML(x_n):
    return np.sqrt(beta)*(1-x_n/beta) -np.exp(x_n)

def f_x_first_derv_OML(x_n):
    return -np.sqrt(beta)*(1/beta) -np.exp(x_n)

def f_x_second_derv_OML(x_n):
    return -np.exp(x_n)

def calc_x_plus_OML(x_n):
    f_x_0 = f_x_OML(x_n)
    f_x_1 = f_x_first_derv_OML(x_n)
    f_x_2 = f_x_second_derv_OML(x_n)
    x_plus = x_n - ((2*f_x_0*f_x_1)/(2.0*(pow(f_x_1,2))-f_x_0*f_x_2))
    return x_plus

def find_phi_D_norm_OML(a_init, root_prec):
    a_n = a_init
    a_plus = calc_x_plus_OML(a_init)

    while(np.abs(a_plus - a_n) > root_prec):
        a_n = a_plus
        a_plus = calc_x_plus_OML(a_n)
    
    W_0 = (a_n + a_plus)/2
    return W_0


def f_x_Sheath( x_n,  A,  B):
    return A*np.exp(B + x_n) - (2/(2*B - 1))*x_n - 1

def f_x_first_derv_Sheath( x_n,  A,  B):
    return A*np.exp(B + x_n) - (2/(2*B - 1))

def f_x_second_derv_Sheath( x_n,  A,  B):
    return A*np.exp(B + x_n)

def calc_x_plus_Sheath( x_n,  A,  B):
    f_x_0 = f_x_Sheath(x_n,A,B)
    f_x_1 = f_x_first_derv_Sheath(x_n,A,B)
    f_x_2 = f_x_second_derv_Sheath(x_n,A,B)
    x_plus = x_n - ((2*f_x_0*f_x_1)/(2.0*((f_x_1**2))-f_x_0*f_x_2))
    return x_plus

def  find_phi_D_norm_OML_Sheath( a_init, z, root_prec):
    A = np.sqrt((8*k_b*T_e)/((v_B**2)*np.pi*m_e))
    B = ((2*e_charge*phi_wall_z)/(k_b*T_e))*((z/z_se - 1)**2)
    a_n = a_init
    a_plus = calc_x_plus_Sheath(a_init,A,B)

    while(np.abs(a_plus - a_n) > root_prec):
        a_n = a_plus
        a_plus = calc_x_plus_Sheath(a_n,A,B)
    W_0 = (a_n + a_plus)/2
    return W_0

def OML_charge(z):
    if(z > z_se):
        Phi_D_norm = find_phi_D_norm_OML(a_0_OML,root)
    
    else:
        Phi_D_norm = find_phi_D_norm_OML_Sheath(a_0_Sheath,z,root)
    phi_grain =  (k_b*T_e*Phi_D_norm)/(e_charge)
    return 4.0*np.pi*epsilon_0*grain_R*phi_grain

def calc_dust_grain_force( d_z,  low_z,  high_z):
    height_list = []
    force_list = []
    for z in np.linspace(low_z,high_z,points):
        charge = OML_charge(z)
        E = 2*k_z_restore*(z - z_se)
        force = charge*E
        height_list.append(z/lambda_D)
        force_list.append(force)
    return [height_list,force_list]



data = calc_dust_grain_force(d_z,  low_z,  high_z)
sheathe_edge = [z_se/lambda_D]*len(data[1])
grav = [9.81*m_D]*len(data[1])

plt.figure()
plt.grid()
plt.title("Force on dust grain vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("Force")
plt.plot(data[0],data[1])
plt.plot(sheathe_edge,data[1],"-" ,color='grey', linestyle='dashed')
plt.plot(data[0],grav,"-" ,color='black', linestyle='dashed')
plt.savefig("Figures/Force_on_dust_grain_vs_Height.png")
plt.show()
    
