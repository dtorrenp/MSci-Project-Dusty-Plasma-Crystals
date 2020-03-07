import numpy as np
import sys
import matplotlib.pyplot as plt#import module used to produce graphs
#import scipy as sp
from scipy.optimize import minimize_scalar
from scipy.special import iv
import matplotlib.animation as animation
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
#wake_potential_below = 1*lambda_D
#wake_charge_multiplier = 1.0
#drop_height = 10*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 25.0*lambda_D

T = 1/(13.6*1e6)
dT = T/3
dz = lambda_D/10
omega = (2*np.pi)/T#TAKEN FROM NITTER PAPER 13.6M Hz
V_RF = -50

low_z = 0
high_z = 100*lambda_D

Y_DC_epsilon = 1e-5
Y_epsilon = 0
Q_epsilon = 0
Y_epsilon_bisection = 1e-5
Q_epsilon_bisection = 1e-5
#u_epsilon_bisection = 1e-5
#sigma = 3e-19
#n_s = n_i0*np.exp(0.5)
#alpha = 0.134 #((epsilon_0*k_b*T_e)/(n_s*(e_charge**2)) )**0.5*n_n0*sigma
#n_d = (dust_grain_max_input)/(np.pi*pow(container_radius,2)*drop_height)
#print("alpha = ", alpha)

#%%

def produce_z_vals():
    z_list = []
    for z in np.arange(low_z,high_z,dz):
        z_list.append(z)
    return z_list

def produce_T_vals():
    T_list = []
    for t in np.arange(0,T,dT):
        T_list.append(t)
    return T_list
#%%
#produce z and t
z_list = produce_z_vals()
T_list = produce_T_vals()
#%%
#produce the dc wall voltage
def produce_V_DC(V_RF):
    V_DC = 0.5*np.log(2*np.pi*m_e/m_i) - np.log(iv(0,V_RF))
    return V_DC

V_DC = produce_V_DC(V_RF)
#%%
#produce Y_z_DC and inital Q guess, no time dependance of electrons

def produce_Y_z_init(Y_Wall_DC):
    Y_z_init = [Y_Wall_DC]
    for v in np.arange(len(z_list) - 1):
        Y_z_init.append(step_Y_z_DC(Y_z_init[v]))
    return np.asarray(Y_z_init)

def step_Y_z_DC(Y):
    k1 = f_der_Y_DC(Y)*dz
    k2 = f_der_Y_DC(Y + k1/2)*dz
    k3 = f_der_Y_DC(Y + k2/2)*dz
    k4 = f_der_Y_DC(Y + k3)*dz
    Y_1 = Y + (k1 + 2*k2 + 2*k3 + k4)/6.0
    return Y_1

def f_der_Y_DC(Y):
    return (2*(np.exp(Y) + (1-2*Y)**0.5 - 2))**0.5

print("V_DC = ",V_DC)
Y_DC_z_no_e = produce_Y_z_init(V_DC)
Q_init_t = np.asarray([f_der_Y_DC(Y_DC_z_no_e[0])]*len(T_list))

plt.figure()
plt.grid()
plt.title("Electric Potential: RF_collisionless")
plt.xlabel("z")
plt.ylabel("Y")
plt.plot(np.asarray(z_list)/lambda_D,Y_DC_z_no_e, label = "original")

plt.show()
exit()
#%%input into poissons equation for a second time now with time dependant electrons

def f_der_Y_z(Y_z_n,Q_z_n,Y_DC_n):
    a = Q_z_n
    b = np.exp(Y_z_n) - (1-2*Y_DC_n)**0.5
    return np.asarray([a,b])

def step_Y_z(Y_z_0,Q_z_0,Y_DC_z):
    Y_DC_z_half = (Y_DC_z[1] + Y_DC_z[0])/dz
    k1 = f_der_Y_z(Y_z_0,Q_z_0, Y_DC_z[0])*dz
    k2 = f_der_Y_z(Y_z_0 + k1[0]/2,Q_z_0 + k1[1]/2, Y_DC_z_half)*dz
    k3 = f_der_Y_z(Y_z_0 + k2[0]/2,Q_z_0 + k2[1]/2, Y_DC_z_half)*dz
    k4 = f_der_Y_z(Y_z_0 + k3[0],Q_z_0 + k3[1], Y_DC_z[1])*dz
    Y_z_1 = Y_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
    Q_z_1 = Q_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
    return np.asarray([Y_z_1,Q_z_1])

def produce_Y_z(Y_z_init,Q_z_init,Y_DC_vals):
    Y_z = [Y_z_init]
    Q_z = [Q_z_init]
    for v in np.arange(len(z_list) - 1):
        Y_DC_select = [Y_DC_vals[v],Y_DC_vals[v+1]]
        a = step_Y_z(Y_z[v],Q_z[v],Y_DC_select)
        Y_z.append(a[0])
        Q_z.append(a[1])
    return np.asarray([Y_z,Q_z])

def produce_Y_z_t(Y_DC_z,Q_init_t):
    Y_z_t=[]
    Q_z_t=[]
    for i in np.arange(len(T_list)):
        print("HEEYEY")
        Y_init_val = Y_DC_z[0] +  V_RF*np.sin(omega*T_list[i]) 
        Q_init_val = Q_init_t[i]
        res = produce_Y_z(Y_init_val,Q_init_val,Y_DC_z)
        Y_last = res[0][-1]
        print(Y_last)
        #BRACKET SOLUTION ABOUT
        if(Y_last < Y_epsilon):
            while(Y_last < Y_epsilon):
                Q_init_val += Q_init_val*0.1
                res = produce_Y_z(Y_init_val,Q_init_val,Y_DC_z)
                Y_last = res[0][-1]
                print(Y_last)
            a = Q_init_val - Q_init_val*0.1
            b = Q_init_val
        else:
            while(Y_last > Y_epsilon):
                Q_init_val -= Q_init_val*0.1
                res = produce_Y_z(Y_init_val,Q_init_val,Y_DC_z)
                Y_last = res[0][-1]
                print(Y_last)
            a = Q_init_val
            b = Q_init_val + Q_init_val*0.1
        #a produces Y_a which will be negative
        #BISECTION
        print("bisection")
        Y_last_a = Y_last
        Y_last_epsilon = np.abs(Y_last_a)
        c = a
        Y_last_c = 0
        while(Y_last_epsilon > Y_epsilon_bisection):
            c = (a+b)/2
            res = produce_Y_z(Y_init_val,c,Y_DC_z)
            Y_last_c = res[0][-1]
            if (Y_last_c*Y_last_a > 0):
                a = c
                Y_last_a = Y_last_c
            else:
                b = c
            print(Y_last_c)
            Y_last_epsilon = np.abs(Y_last_c)
        Q_init_t[i] = c
        final_Y_z = produce_Y_z(Y_init_val,c,Y_DC_z)
        Y_z_t.append(final_Y_z[0])
        Q_z_t.append(final_Y_z[1])
    return np.asarray([Y_z_t,Q_z_t])

def produce_Y_DC_converge(Y_DC_z_init,Q_init_t):
    #run the Y_z_t once
    a = produce_Y_z_t(Y_DC_z_init,Q_init_t)
    Y_z_t_converge = a[0]
    Q_z_t_converge = a[1]
    Y_DC_converge = produce_Y_DC(Y_z_t_converge)
    Q_DC_converge = produce_Y_DC(Q_z_t_converge)
    return [Y_DC_converge,Q_DC_converge,Y_z_t_converge,Q_z_t_converge,]
#%%
#average over time
def produce_Y_DC(Y):
    Y_DC = []
    for v in np.arange(len(z_list)):
        a = []
        for i in np.arange(len(T_list)):
            a.append(Y[i][v])
        Y_DC.append(np.mean(a))
    return np.asarray(Y_DC)
#%%
#call the functions
results = produce_Y_DC_converge(Y_DC_z_no_e,Q_init_t)
Y_DC_z = results[0]
Q_DC_z = results[1]
Y_z_t = results[2]
Q_z_t = results[3]

plt.figure()
plt.grid()
plt.title("Electric Potential: RF_collisionless")
plt.xlabel("z")
plt.ylabel("Y")
for i in np.arange(len(T_list)):
    plt.plot(np.asarray(z_list)/lambda_D,Y_z_t[i] , label = "t = " + str(T_list[i]))
plt.plot(np.asarray(z_list)/lambda_D,Y_DC_z, label = "time averaged" )
plt.plot(np.asarray(z_list)/lambda_D,Y_DC_z_no_e, label = "original")
plt.legend()
plt.savefig("Figures/Electric_potential_vs_z_RF_collisionless.png")

plt.figure()
plt.grid()
plt.title("Electric Field: RF_collisionless")
plt.xlabel("z")
plt.ylabel("E")
for i in np.arange(len(T_list)):
    plt.plot(np.asarray(z_list)/lambda_D,Q_z_t[i] , label = "t = " + str(T_list[i]))
plt.plot(np.asarray(z_list)/lambda_D,Q_DC_z, label = "time averaged" )
plt.legend()
plt.savefig("Figures/Electric_Field_vs_z_RF_collisionless.png")

plt.show()


#%%
# def collisionless():
#     def produce_V_DC(V_RF):
#         V_DC = 0.5*np.log(2*np.pi*m_e/m_i) - np.log(iv(0,V_RF))
#         return V_DC

#     def produce_Y_z_init(Y_Wall_DC):
#         Y_z_init = []
#         for v in np.arange(len(z_list)):
#             Y_z_init.append(Y_Wall_DC *(z_list[v]/z_list[-1] - 1)**2)
#         return np.asarray(Y_z_init)

#     def produce_Y_z_t(Y_DC,Q_init):
#         Y_z_t=[]
#         for i in np.arange(len(T_list)):
#             Y_init_val = Y_DC[0] +  V_RF*np.sin(omega*T_list[i]) 
#             Q_init_val = Q_init[i]
#             res = produce_Y_z(Y_init_val,Q_init_val,Y_DC)
#             Y_last = res[0][-1]
#             while(Y_last < Y_epsilon):
#                 Q_init_val += Q_init_val*0.5
#                 res = produce_Y_z(Y_init_val,Q_init_val,Y_DC)
#                 Y_last = res[0][-1]
#             #BISECTION
#             a = Q_init_val
#             b = Q_init_val/2
#             Y_last_a = Y_last
#             Y_last_epsilon = np.abs(Y_last)
#             Y_last_c = 0
#             while(Y_last_epsilon > Y_epsilon_bisection):
#                 c = (a+b)/2
#                 res = produce_Y_z(Y_init_val,c,Y_DC)
#                 Y_last_c = res[0][-1]
#                 if (Y_last_c*Y_last_a > 0):
#                     a = c
#                     Y_last_a = Y_last_c
#                 else:
#                     b = c
#                 Y_last_epsilon = np.abs(Y_last_c)
#             Q_init[i] = b
#             Y_z_t.append(produce_Y_z(Y_init_val,b,Y_DC)[0])
#         return np.asarray([Y_z_t,Q_init])
    
#     def produce_Y_z(Y_z_init,Q_z_init,Y_DC_vals):
#         Y_z = [Y_z_init]
#         Q_z = [Q_z_init]
#         for v in np.arange(len(z_list) - 1):
#             Y_DC_select = [Y_DC_vals[v],Y_DC_vals[v+1]]
#             a = step_Y_z(Y_z[v],Q_z[v],Y_DC_select)
#             Y_z.append(a[0])
#             Q_z.append(a[1])
#         return np.asarray([Y_z,Q_z])

#     def step_Y_z(Y_z_0,Q_z_0,Y_DC_z):
#         Y_DC_z_half = (Y_DC_z[1] + Y_DC_z[0])/dz
#         k1 = f_der_Y_z(Y_z_0,Q_z_0, Y_DC_z[0])*dz
#         k2 = f_der_Y_z(Y_z_0 + k1[0]/2,Q_z_0 + k1[1]/2, Y_DC_z_half)*dz
#         k3 = f_der_Y_z(Y_z_0 + k2[0]/2,Q_z_0 + k2[1]/2, Y_DC_z_half)*dz
#         k4 = f_der_Y_z(Y_z_0 + k3[0],Q_z_0 + k3[1], Y_DC_z[1])*dz
#         Y_z_1 = Y_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
#         Q_z_1 = Q_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
#         return np.asarray([Y_z_1,Q_z_1])

#     def f_der_Y_z(Y_z_n,Q_z_n,Y_DC_n):
#         a = Q_z_n
#         b = np.exp(Y_z_n) - (1-2*Y_DC_n)**0.5
#         return np.asarray([a,b])
    
#     def produce_Y_DC_converge(Y_DC_z_init,Q_init):
#         print("1")
#         a = produce_Y_z_t(Y_DC_z_init,Q_init)
#         print("2")
#         Y_z_t_converge = a[0]
#         Q_init_converge = a[1]
#         Y_DC_converge = produce_Y_DC(Y_z_t_converge)
#         print("3")
#         Y_difference = produce_Y_difference(0,Y_DC_converge)
#         print("4")
#         reps = 0
#         while(Y_difference > Y_DC_epsilon):
#             reps+=1
#             print("5")
#             Y_DC_temp = np.copy(Y_DC_converge)
#             a = produce_Y_z_t(Y_DC_converge,Q_init_converge)
#             print("6")
#             Y_z_t_converge = a[0]
#             Q_init_converge = a[1]
#             Y_DC_converge = produce_Y_DC(Y_z_t_converge)
#             print("7")
#             Y_difference = produce_Y_difference(Y_DC_temp,Y_DC_converge)
#             print("8")

#         return np.asarray([Y_z_t_converge,Y_DC_converge,Q_init_converge])

#     def produce_Y_DC(Y):
#         Y_DC = []
#         for v in np.arange(len(z_list)):
#             a = []
#             for i in np.arange(len(T_list)):
#                 a.append(Y[i][v])
#             Y_DC.append(np.mean(a))
#         return np.asarray(Y_DC)

#     def produce_Y_difference(Y_bar_temp, Y_bar): 
#         a = Y_bar_temp -  Y_bar
#         return np.linalg.norm(a)/np.linalg.norm(Y_bar)
#     #produce z and t
#     z_list = produce_z_vals()
#     T_list = produce_T_vals()
#     #produce the dc wall voltage
#     V_DC = produce_V_DC(V_RF)
#     #produce Y_z_DC, no time dependance of electrons
#     Y_z_init = produce_Y_z_init(V_DC)

#     init_grad = ((Y_z_init[1]-Y_z_init[0])/dz)

#     Q_init = np.asarray([1e-3*init_grad]*len(T_list))
#     results = produce_Y_DC_converge(Y_z_init,Q_init)
#     Y_z_t = results[0]
#     Y_DC_z = results[1]
#     Q_final = results[2]

#     def produce_electric_field(Y,Q_init_val):
#         Q = [Q_init_val]
#         for v in np.arange(1,len(z_list) - 1):
#             Q.append((Y[v+1] - Y[v-1])/(2*dz))
#         Q.append(Q[-1])
#         return Q

#     E = []
#     for i in np.arange(len(T_list)):
#         E.append(produce_electric_field(Y_z_t[i],Q_final[i]))
#     E_average = produce_Y_DC(E)

#     plt.figure()
#     plt.grid()
#     plt.title("Electric Potential: RF_collisionless")
#     plt.xlabel("z")
#     plt.ylabel("Y")
#     for i in np.arange(len(T_list)):
#         plt.plot(np.asarray(z_list),Y_z_t[i] , label = "t = " + str(T_list[i]))
#     plt.plot(np.asarray(z_list),Y_DC_z, label = "time averaged" )
#     plt.plot(np.asarray(z_list),Y_z_init, label = "original")
#     plt.legend()
#     plt.savefig("Figures/Electric_potential_vs_z_RF_collisionless.png")

#     plt.figure()
#     plt.grid()
#     plt.title("Electric Field: RF_collisionless")
#     plt.xlabel("z")
#     plt.ylabel("E")
#     for i in np.arange(len(T_list)):
#         plt.plot(np.asarray(z_list),E[i] , label = "t = " + str(T_list[i]))
#     plt.plot(np.asarray(z_list),E_average, label = "time averaged" )
#     plt.legend()
#     plt.savefig("Figures/Electric_Field_vs_z_RF_collisionless.png")
#     return [np.asarray(z_list),Y_DC_z]