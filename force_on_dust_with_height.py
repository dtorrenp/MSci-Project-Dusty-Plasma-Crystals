import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from scipy.special import erf
from scipy.special import iv

#CONSTANTS TO FUCK ABOUT WITH
n_e0 = 1.0e15#electron and ion densities in bulk plasma
n_i0 = 1.0e15
phi_wall_z = -100.0#volts

#CONSTANTS DEPENDANT ON ACTUAL PHYSICS
Z = 18
g_z = 9.81#gravity
e_charge = 1.6*1e-19#magnitude of e charge
i_charge = 1.6*1e-19#magnitude of i charge DOES THIS NEED TO CHANGE WHEN USED IN THE ION DRAG?????
grain_R = 7*1e-6
m_i = Z*1.67*1e-27
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
a_0_Sheath = -2.5#intial guess for halley's method
root = 1.0e-4#preciscion of root finding method used to get dust charge
#print(low_z,high_z,points)
#print(type(low_z), type(high_z), type(points))
#container_dust_dist_creation = np.sqrt(container_radius/(2*lambda_D))
#k_z_restore = -2.0*phi_wall_z/(z_se**2)#WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
v_B = ((3*k_b*T_e/m_i)**0.5)

#%%
######CHARGE

#freq = 13.56*1e6 #frequency [Hz] (13.56MHz as quoted from Nitter)
#omega = (2*np.pi)*freq #angular frequency [rad/s]
#sigma = 3e-19
#n_s = n_e0*np.exp(-0.5) #electron density at the sheath edge; see Nitter (2)
#alpha = 0.134 #((epsilon_0*k_b*T_e)/(n_s*(e_charge**2)) )**0.5*n_n0*sigma
#T = 1/freq #period
V_RF = 50 #(as quoted from Nitter)
#dT = T/3
dz = 1/100
drop_height_numerical = 50

def produce_z_list():
    """Produces discrete z-coordinates in steps of d_z and units of Debye lengths."""
    z_list = []
    for z in np.arange(0, drop_height_numerical, dz):
        z_list.append(z)
    return z_list

# def produce_T_list():
#     """Produces discrete time coordinates in steps of dT."""
#     T_list = []
#     for t in np.arange(0,T,dT):
#         T_list.append(t)
#     return T_list

###
z_list = produce_z_list()
#T_list = produce_T_list()

#%%

#STEP 0: Get initial Y_z (and also initial Q)
"""The initial Y_z DC sheath potential is used as the basis for both the DC and RF cases. 
This initial potential guess is based on Cameron's collisionless equations."""

V_DC = 0.5*np.log(2*np.pi*m_e/m_i) - np.log(iv(0,V_RF)) #Produces the self-bias V_DC in both the RF case and DC case (where V_RF = 0). V in normalised units of eV/kT.
Y_epsilon = np.abs(V_DC*1e-3)

def f_Y_z_init(Y_z_0):
    """RK4 derivative function"""
    return (2*(np.exp(Y_z_0) + (1-2*Y_z_0)**0.5 - 2))**0.5 #Cameron (2.36)

def step_Y_z_init(Y_z_0):
    """RK4 step of size dz for Y_z_init"""
    k1 = dz * f_Y_z_init(Y_z_0)
    k2 = dz * f_Y_z_init(Y_z_0 + k1/2)
    k3 = dz * f_Y_z_init(Y_z_0 + k2/2)
    k4 = dz * f_Y_z_init(Y_z_0 + k3)
    Y_z_1 = Y_z_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
    return Y_z_1

###
Y_z_init=[V_DC] #Producing Y_z_init with boundary condition Y_wall = V_DC
Q_init = [f_Y_z_init(V_DC)] #Producing  Q_init i.e. dY_z/dz
for j in np.arange(len(z_list) - 1):
    Y_z_init.append(step_Y_z_init(Y_z_init[j]))
    Q_init.append(f_Y_z_init(Y_z_init[j]))

#STEP 0b: Find sheath edge L

###
L = np.where(np.abs(Y_z_init) < Y_epsilon)[0][0]
len_z_rep = len(Y_z_init) - L

for i in np.arange(len_z_rep):
    Y_z_init[L + i] = 0
    Q_init[L + i] = 0

plt.figure()
plt.grid()
plt.title("Electric field on dust grain vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("Electric field")
plt.plot(z_list,Q_init)
plt.savefig("Figures/Electric_field_on_dust_grain_vs_Height.png")

plt.figure()
plt.grid()
plt.title("Electric Potential on dust grain vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("Electric Potential")
plt.plot(z_list,Y_z_init)
plt.savefig("Figures/Electric_potential_on_dust_grain_vs_Height.png")
#%%

#ION VELOCITY

def v_i_z_calc(Y,z):
    if z < z_list[L]:
        return -(v_B**2 - 2*i_charge*Y/m_i)**0.5
    else:
        return 0

v_i_list = []
for i in np.arange(len(z_list)):
    v_i_list.append(v_i_z_calc(Y_z_init[i],z_list[i]))

plt.figure()
plt.grid()
plt.title("Ion velocity vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("Ion Velcocity/v_B")
plt.plot(z_list,np.asarray(v_i_list)/v_B)
plt.savefig("Figures/Ion_vel_norm_vs_Height.png")

#%%

######CHARGE

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

def  find_phi_D_norm_OML_Sheath(a_init,root_prec,Phi_DC):
    A = np.sqrt((8*k_b*T_e)/((v_B**2)*np.pi*m_e))
    a_n = a_init
    a_plus = calc_x_plus_Sheath(a_init,A,Phi_DC)

    while(np.abs(a_plus - a_n) > root_prec):
        a_n = a_plus
        a_plus = calc_x_plus_Sheath(a_n,A,Phi_DC)
    W_0 = (a_n + a_plus)/2
    return W_0

Phi_D_norm_list = []
def OML_charge(Phi_DC):
    
    Phi_D_norm = find_phi_D_norm_OML_Sheath(a_0_Sheath,root,Phi_DC)
    Phi_D_norm_list.append(Phi_D_norm)
    phi_grain =  (k_b*T_e*Phi_D_norm)/(e_charge)
    #print(Phi_DC,Phi_D_norm,phi_grain,4.0*np.pi*epsilon_0*grain_R*phi_grain)
    return 4.0*np.pi*epsilon_0*grain_R*phi_grain

charge_list = []
for i in np.arange(len(z_list)):
    charge_list.append( OML_charge(Y_z_init[i])   )

plt.figure()
plt.grid()
plt.title("charge on dust grain vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("charge")
plt.plot(z_list,charge_list, label = "collisionless")
plt.savefig("Figures/Charge_on_dust_grain_vs_Height.png")
plt.legend()

plt.figure()
plt.grid()
plt.title("phiD vs phi DC")
plt.xlabel("phi DC")
plt.ylabel("phi D")
plt.plot(Y_z_init,Phi_D_norm_list)
plt.savefig("Figures/Charge_on_dust_grain_vs_phi_DC.png")
plt.legend()


#%%

#ION DRAG

def ion_drag(z,charge_z,v_i):#IT SHOULD DOUBEL BACK SURELY COS OF THE CHARGE FLIP?
    #charge_z = 0# THE COLLECTION TERM DOMIANTES DRAMATICALLY
    z_ion = (charge_z*i_charge/(4*np.pi*epsilon_0))*(1/(grain_R*k_b*T_e))
    t_ion = beta
    n = n_e0 
    m = m_i
    v_ti = -(k_b*T_i/m_i)**0.5
    u = v_i/v_ti
    R = (charge_z*i_charge)/((4*np.pi*epsilon_0)*k_b*T_i*(1+u**2))
    beta_bar = R/lambda_D#I AM NOT USIG THE IMPORVED LAMBDA MAY BE SOEMTHIGN TO TRY OUT
    LAMBDA = np.abs((beta_bar + 1)/(beta_bar + (grain_R/lambda_D)))
    #print(z,z_ion,charge_z,u,R,LAMBDA)
    F_i = -(2*np.pi)**0.5*grain_R**2*n*m*v_ti**2*( (np.pi/2)**0.5*erf(u/2**0.5)*(1+u**2+(1-u**-2)*(1+2*z_ion*t_ion) + 4*z_ion**2*t_ion**2*u**-2*np.log(LAMBDA)) + u**-1*(1 + 2*z_ion*t_ion+u**2 -  4*z_ion**2*t_ion**2*np.log(LAMBDA))*np.exp(-u**2/2)  )
    #F_i = -(2*np.pi)**0.5*grain_R**2*n*m*v_ti**2*(  u**-1*(1 + 2*z_ion*t_ion+u**2 -  4*z_ion**2*t_ion**2*np.log(LAMBDA))*np.exp(-u**2/2)   )
    #print(u**(-0.5), z_ion*t_ion,np.exp(-
    # u**2/2) )
    #print(F_i)
    return F_i

ion_drag_list = []
for i in np.arange(len(z_list)):
    if i >= L:
        ion_drag_list.append(0)
    else:
        ion_drag_list.append(ion_drag(z_list[i],charge_list[i],v_i_list[i]))

#print(ion_drag_list)
plt.figure()
plt.grid()
plt.title("Ion drag vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("Ion drag")
plt.plot(z_list,ion_drag_list)
plt.savefig("Figures/ion_drag_on_dust_grain_vs_Height.png")

#%%
#####FORCE
force_list = []
for i in np.arange(len(z_list)):
    force_list.append(charge_list[i]*Q_init[i] - g_z*m_D)

grav = [9.81*m_D]*len(z_list)

plt.figure()
plt.grid()
plt.title("Force on dust grain vs Height")
plt.xlabel("z/lambda_D")
plt.ylabel("Force")
plt.plot(z_list,force_list)
plt.plot(z_list,grav,"-" ,color='black', linestyle='dashed')
plt.savefig("Figures/Force_on_dust_grain_vs_Height.png")

plt.show()
