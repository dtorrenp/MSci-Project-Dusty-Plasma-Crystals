import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs
from scipy.special import erf
from scipy.special import iv
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
#CONSTANTS TO FUCK ABOUT WITH
n_e0 = 1.0e15#electron and ion densities in bulk plasma
n_i0 = 1.0e15

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
v_B = ((3*k_b*T_e/m_i)**0.5)
V_RF = 100 #(as quoted from Nitter)
dz = 1/100
drop_height_numerical = 40
#%%
######CHARGE
def produce_z_list():
    """Produces discrete z-coordinates in steps of d_z and units of Debye lengths."""
    z_list = []
    for z in np.arange(0, drop_height_numerical, dz):
        z_list.append(z)
    return z_list

###
z_list = produce_z_list()

#%%

#STEP 0: Get initial Y_z (and also initial Q)
"""The initial Y_z DC sheath potential is used as the basis for both the DC and RF cases. 
This initial potential guess is based on Cameron's collisionless equations."""

def inital_Y_z(V_RF):

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

    return [Y_z_init,Q_init]

#%%
######CHARGE

def calc_charge(Y_z_init):
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
        charge_list.append(OML_charge(Y_z_init[i]))

    return charge_list

#%%
def calc_sheath_force(Q_z_list,charge_list):

    sheath_force = []

    for i in np.arange(len(z_list)):
        E_field = -((k_b*T_e)/(e_charge*lambda_D))*Q_z_list[i]
        sheath_force.append(charge_list[i]*E_field)
    
    return sheath_force



#%%

V_RF_list = [50,55,60,65,70]

Y_z_list = []
Q_z_list = []
charge_list = []
sheath_force_list = []

for i in np.arange(len(V_RF_list)):
    Y_z_init_and_Q_init = inital_Y_z(V_RF_list[i])
    Y_z_list.append(Y_z_init_and_Q_init[0])
    Q_z_list.append(Y_z_init_and_Q_init[1])
    charge_list.append(calc_charge(Y_z_init_and_Q_init[0]))
    sheath_force_list.append(  calc_sheath_force(Q_z_list[i],charge_list[i])      )

plt.figure()
plt.grid()
plt.title(r"Q on dust grain vs Height")
plt.xlabel(r"z/$\lambda_D$")
plt.ylabel(r"Q")
for i in np.arange(len(V_RF_list)):
    plt.plot(z_list,Y_z_list[i], label = r"$V_{RF}$ = " + str(V_RF_list[i]))
plt.legend()

plt.figure()
plt.grid()
plt.title(r"Normalised Electric Potential on dust grain vs Height")
plt.xlabel(r"z/$\lambda_D$")
plt.ylabel(r"Normalised Electric Potential")
for i in np.arange(len(V_RF_list)):
    plt.plot(z_list,Y_z_list[i], label = r"$V_{RF}$ = " + str(V_RF_list[i]))
plt.legend()

plt.figure()
plt.grid()
plt.title(r"Charge on dust grain vs Height")
plt.xlabel(r"z/$\lambda_D$")
plt.ylabel(r"Charge")
for i in np.arange(len(V_RF_list)):
    plt.plot(z_list,charge_list[i], label = r"$V_{RF}$ = " + str(V_RF_list[i]))
plt.legend()

grav = [9.81*m_D]*len(z_list)
#print(-9.81*m_D)
#print(sheath_force_list[3])

plt.figure()
plt.grid()
plt.title(r"Sheath Force on dust grain vs Height")
plt.xlabel(r"z/$\lambda_D$")
plt.ylabel(r"Sheath Force")
plt.plot(z_list,grav,"-" ,color='black', linestyle='dashed', label = "Gravity")
for i in np.arange(len(V_RF_list)):
    plt.plot(z_list,sheath_force_list[i], label = r"$V_{RF}$ = " + str(V_RF_list[i]))
plt.legend()

fig, ax1 = plt.subplots() # create a new figure with a default 111 subplot
ax1.set_title(r"Sheath Force on dust grain vs Height")
ax1.set_xlabel(r"z/$\lambda_D$")
ax1.set_ylabel(r"Sheath Force")
ax1.plot(z_list,grav,"-" ,color='black', linestyle='dashed', label = "Gravity")
for i in np.arange(len(V_RF_list)):
    ax1.plot(z_list,sheath_force_list[i], label = r"$V_{RF}$ = " + str(V_RF_list[i]))

axins = zoomed_inset_axes(ax1, 20, loc=5) # zoom-factor: 2.5, location: #
axins.plot(z_list,grav,"-" ,color='black', linestyle='dashed', label = "Gravity")
for i in np.arange(len(V_RF_list)):
    axins.plot(z_list,sheath_force_list[i], label = r"$V_{RF}$ = " + str(V_RF_list[i]))
x1, x2, y1, y2 = 18.5, 20, -2*10**-10, 2*10**-10 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
mark_inset(ax1, axins, loc1=1, loc2=4, fc="none", ec="0.5")
ax1.grid()
axins.grid()
#axins.set_aspect(1.0)
ax1.legend(ncol = 3)
plt.tight_layout()

plt.show()
