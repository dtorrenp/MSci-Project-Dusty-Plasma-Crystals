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
wake_potential_below = 1*lambda_D
wake_charge_multiplier = 1.0
drop_height = 10*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 25.0*lambda_D

T = 1/(13.6*1e6)
dT = T/3
dz = lambda_D/100
omega = (2*np.pi)/T#TAKEN FROM NITTER PAPER 13.6M Hz
V_RF = 50
low_z = 0
high_z = drop_height
Y_DC_epsilon = 1e-5
Y_epsilon = 0
Q_epsilon = 0
Y_epsilon_bisection = 1e-5
Q_epsilon_bisection = 1e-5
u_epsilon_bisection = 1e-5
sigma = 3e-19
n_s = n_i0*np.exp(0.5)
alpha = 0.134 #((epsilon_0*k_b*T_e)/(n_s*(e_charge**2)) )**0.5*n_n0*sigma
n_d = (dust_grain_max_input)/(np.pi*pow(container_radius,2)*drop_height)
print("alpha = ", alpha)

#%%

def produce_z_vals():
    z_list = []
    for z in np.arange(low_z,high_z,dz):
        z_list.append(z/lambda_D)
    return z_list

def produce_T_vals():
    T_list = []
    for t in np.arange(0,T,dT):
        T_list.append(t)
    return T_list

#%%
def collisionless():
    def produce_V_DC(V_RF):
        V_DC = 0.5*np.log(2*np.pi*m_e/m_i) - np.log(iv(0,V_RF))
        return V_DC

    def produce_Y_z_init(Y_Wall_DC):
        Y_z_init=[Y_Wall_DC]
        for v in np.arange(len(z_list) - 1):
            Y_z_init.append(step_Y_z_init(Y_z_init[v]))
        return Y_z_init

    def step_Y_z_init(Y_z_0):
        k1 = f_der_Y_z_init(Y_z_0)*dz
        k2 = f_der_Y_z_init(Y_z_0 + k1/2)*dz
        k3 = f_der_Y_z_init(Y_z_0 + k2/2)*dz
        k4 = f_der_Y_z_init(Y_z_0 + k3)*dz
        Y_z_1 = Y_z_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
        return Y_z_1

    def f_der_Y_z_init(Y_z_n):
        return (2*(np.exp(Y_z_n) + (1-2*Y_z_n)**0.5 -2))**0.5

    def produce_Y_z_t(Y_DC,Q_init):
        Y_z_t=[]
        for i in np.arange(len(T_list)):
            Y_init_val = Y_DC[0] +  V_RF*np.sin(omega*T_list[i]) 
            Q_init_val = Q_init[i]
            res = produce_Y_z(Y_init_val,Q_init_val,Y_DC)
            Y_last = res[0][-1]
            while(Y_last < Y_epsilon):
                Q_init_val += Q_init_val*0.5
                res = produce_Y_z(Y_init_val,Q_init_val,Y_DC)
                Y_last = res[0][-1]
            #BISECTION
            a = Q_init_val
            b = Q_init_val/2
            Y_last_a = Y_last
            Y_last_epsilon = np.abs(Y_last)
            Y_last_c = 0
            while(Y_last_epsilon > Y_epsilon_bisection):
                c = (a+b)/2
                res = produce_Y_z(Y_init_val,c,Y_DC)
                Y_last_c = res[0][-1]
                if (Y_last_c*Y_last_a > 0):
                    a = c
                    Y_last_a = Y_last_c
                else:
                    b = c
                Y_last_epsilon = np.abs(Y_last_c)
            Q_init[i] = b
            Y_z_t.append(produce_Y_z(Y_init_val,b,Y_DC)[0])
        return np.asarray([Y_z_t,Q_init])
    
    def produce_Y_z(Y_z_init,Q_z_init,Y_DC_vals):
        Y_z = [Y_z_init]
        Q_z = [Q_z_init]
        for v in np.arange(len(z_list) - 1):
            Y_DC_select = [Y_DC_vals[v],Y_DC_vals[v+1]]
            a = step_Y_z(Y_z[v],Q_z[v],Y_DC_select)
            Y_z.append(a[0])
            Q_z.append(a[1])
        return np.asarray([Y_z,Q_z])

    def step_Y_z(Y_z_0,Q_z_0,Y_DC_z):
        Y_DC_z_half = (Y_DC_z[1] + Y_DC_z[0])/dz
        k1 = f_der_Y_z(Y_z_0,Q_z_0, Y_DC_z[0])*dz
        k2 = f_der_Y_z(Y_z_0 + k1[0]/2,Q_z_0 + k1[1]/2, Y_DC_z_half)*dz
        k3 = f_der_Y_z(Y_z_0 + k2[0]/2,Q_z_0 + k2[1]/2, Y_DC_z_half)*dz
        k4 = f_der_Y_z(Y_z_0 + k3[0],Q_z_0 + k3[1], Y_DC_z[1])*dz
        Y_z_1 = Y_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
        Q_z_1 = Q_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
        return np.asarray([Y_z_1,Q_z_1])

    def f_der_Y_z(Y_z_n,Q_z_n,Y_DC_n):
        a = Q_z_n
        b = np.exp(Y_z_n) - (1-2*Y_DC_n)**0.5
        return np.asarray([a,b])
    
    def produce_Y_DC_converge(Y_DC_z_init,Q_init):

        a = produce_Y_z_t(Y_DC_z_init,Q_init)
        Y_z_t_converge = a[0]
        Q_init_converge = a[1]
        Y_DC_converge = produce_Y_DC(Y_z_t_converge)
        Y_difference = produce_Y_difference(0,Y_DC_converge)
        reps = 0
        while(Y_difference > Y_DC_epsilon):
            reps+=1
            Y_DC_temp = np.copy(Y_DC_converge)
            a = produce_Y_z_t(Y_DC_converge,Q_init_converge)
            Y_z_t_converge = a[0]
            Q_init_converge = a[1]
            Y_DC_converge = produce_Y_DC(Y_z_t_converge)
            Y_difference = produce_Y_difference(Y_DC_temp,Y_DC_converge)

        return np.asarray([Y_z_t_converge,Y_DC_converge,Q_init_converge])

    def produce_Y_DC(Y):
        Y_DC = []
        for v in np.arange(len(z_list)):
            a = []
            for i in np.arange(len(T_list)):
                a.append(Y[i][v])
            Y_DC.append(np.mean(a))
        return np.asarray(Y_DC)

    def produce_Y_difference(Y_bar_temp, Y_bar): 
        a = Y_bar_temp -  Y_bar
        return np.linalg.norm(a)/np.linalg.norm(Y_bar)
    
    z_list = produce_z_vals()
    T_list = produce_T_vals()
    V_DC = produce_V_DC(V_RF)
    Y_z_init = produce_Y_z_init(V_DC)
    Q_init = np.asarray([0.1]*len(T_list))
    results = produce_Y_DC_converge(Y_z_init,Q_init)
    Y_z_t = results[0]
    Y_DC_z = results[1]
    Q_final = results[2]

    def produce_electric_field(Y,Q_init_val):
        Q = [Q_init_val]
        for v in np.arange(1,len(z_list) - 1):
            Q.append((Y[v+1] - Y[v-1])/(2*dz))
        Q.append(Q[-1])
        return Q

    E = []
    for i in np.arange(len(T_list)):
        E.append(produce_electric_field(Y_z_t[i],Q_final[i]))
    E_average = produce_Y_DC(E)

    plt.figure()
    plt.grid()
    plt.title("Electric Potential: RF_collisionless")
    plt.xlabel("z")
    plt.ylabel("Y")
    for i in np.arange(len(T_list)):
        plt.plot(np.asarray(z_list),Y_z_t[i] , label = "t = " + str(T_list[i]))
    plt.plot(np.asarray(z_list),Y_DC_z, label = "time averaged" )
    plt.plot(np.asarray(z_list),Y_z_init, label = "original")
    plt.legend()
    plt.savefig("Figures/Electric_potential_vs_z_RF_collisionless.png")

    plt.figure()
    plt.grid()
    plt.title("Electric Field: RF_collisionless")
    plt.xlabel("z")
    plt.ylabel("E")
    for i in np.arange(len(T_list)):
        plt.plot(np.asarray(z_list),E[i] , label = "t = " + str(T_list[i]))
    plt.plot(np.asarray(z_list),E_average, label = "time averaged" )
    plt.legend()
    plt.savefig("Figures/Electric_Field_vs_z_RF_collisionless.png")
    return [np.asarray(z_list),Y_DC_z]

def collisional():
    #STEP 0 GET Y initial
    def produce_V_DC(V_RF):
        V_DC = 0.5*np.log(2*np.pi*m_e/m_i) - np.log(iv(0,V_RF))
        return V_DC

    def produce_Y_z_init(Y_Wall_DC):
        Y_z_init=[Y_Wall_DC]
        for v in np.arange(len(z_list) - 1):
            Y_z_init.append(step_Y_z_init(Y_z_init[v]))
        return Y_z_init

    def step_Y_z_init(Y_z_0):
        k1 = f_der_Y_z_init(Y_z_0)*dz
        k2 = f_der_Y_z_init(Y_z_0 + k1/2)*dz
        k3 = f_der_Y_z_init(Y_z_0 + k2/2)*dz
        k4 = f_der_Y_z_init(Y_z_0 + k3)*dz
        Y_z_1 = Y_z_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
        return Y_z_1

    def f_der_Y_z_init(Y_z_n):
        return (2*(np.exp(Y_z_n) + (1-2*Y_z_n)**0.5 -2))**0.5

    #STEP 1 GET U
    def f_der_u(u_0,Q):
        #print(-(1/u_0)*Q, alpha*u_0)
        return -(1/u_0)*Q + alpha*u_0

    def step_u(u_n,Q_n,h):
        Q_n_half = (Q_n[0] + Q_n[1])/2
        k1 = f_der_u(u_n,Q_n[0])*h
        k2 = f_der_u(u_n + k1/2,Q_n_half)*h
        k3 = f_der_u(u_n + k2/2,Q_n_half)*h
        k4 = f_der_u(u_n + k3,Q_n[1])*h
        u_1 = u_n + (k1 + 2*k2 + 2*k3 + k4)/6.0
        return u_1

    def produce_u_z(Q_row,u_0_init):
        u = [u_0_init]
        for v in np.arange(len(z_list) - 1):
            Q_select = [Q_row[v], Q_row[v+1]]
            u.append(step_u(u[v],Q_select,dz))
        return u

    def produce_u(Q_row,u_0):
        #print("u_0 = ",u_0)
        u = produce_u_z(Q_row,u_0)
        #print("u[-1] + 1 = ",u[-1] + 1)
        while(u[-1] + 1 < 0):
            #print(u[-1])
            u_0 -= u_0*0.5
            u = produce_u_z(Q_row,u_0)
            #print("u[-1] + 1 = ",u[-1] + 1)
        #BISECTION
        a = u_0
        b = 2*u_0
        u_root_a = u[-1] + 1
        u_root_c = 0
        u_root_epsilon = np.abs(u_root_a)
        while(u_root_epsilon > u_epsilon_bisection):
            c = (a+b)/2
            u = produce_u_z(Q_row,c)
            u_root_c = u[-1] + 1
            if (u_root_c*u_root_a > 0):
                a = c
                u_root_a = u_root_c
            else:
                b = c
            u_root_epsilon = np.abs(u_root_c)
            #print("u root = ", u_root_epsilon)
        return u

    #STEP 2 GET Y
    def produce_Y_z_t(Y_DC_wall,Q_init,u):
        Y_z_t=[]
        Q_z_t=[]
        for i in np.arange(len(T_list)):
            Y_init_val = Y_DC_wall +  V_RF*np.sin(omega*T_list[i]) 
            Q_init_val = Q_init[i]
            res = produce_Y_z(Y_init_val,Q_init_val,u)
            Q_last = res[1][-1]
            while((Q_last - alpha) < Q_epsilon):
                Q_init_val += Q_init_val*0.5
                res = produce_Y_z(Y_init_val,Q_init_val,u)
                Q_last = res[1][-1]
            a = Q_init_val
            b = Q_init_val/2
            Q_root_a = Q_last - alpha
            Q_root_epsilon = np.abs(Q_last - alpha)
            Q_root_c = 0
            while(Q_root_epsilon > Q_epsilon_bisection):
                c = (a+b)/2
                res = produce_Y_z(Y_init_val,c,u)
                Q_root_c = res[1][-1] - alpha
                if (Q_root_c*Q_root_a > 0):
                    a = c
                    Q_root_a = Q_root_c
                else:
                    b = c
                Q_root_epsilon = np.abs(Q_root_c)
                #print("Q_root_epsilon = ",Q_root_epsilon)
            Q_init[i] = b
            results_final = produce_Y_z(Y_init_val,b,u)
            Y_z_t.append(results_final[0])
            Q_z_t.append(results_final[1])
        return np.asarray([Y_z_t,Q_z_t])
    
    def produce_Y_z(Y_z_init,Q_z_init,u_vals):
        Y_z = [Y_z_init]
        Q_z = [Q_z_init]
        for v in np.arange(len(z_list) - 1):
            u_select = [u_vals[v],u_vals[v+1]]
            a = step_Y_z(Y_z[v],Q_z[v],u_select)
            Y_z.append(a[0])
            Q_z.append(a[1])
        return np.asarray([Y_z,Q_z])

    def step_Y_z(Y_z_0,Q_z_0,u_z):
        u_z_half = (u_z[1] + u_z[0])/dz
        k1 = f_der_Y_z(Y_z_0,Q_z_0, u_z[0])*dz
        k2 = f_der_Y_z(Y_z_0 + k1[0]/2,Q_z_0 + k1[1]/2, u_z_half)*dz
        k3 = f_der_Y_z(Y_z_0 + k2[0]/2,Q_z_0 + k2[1]/2, u_z_half)*dz
        k4 = f_der_Y_z(Y_z_0 + k3[0],Q_z_0 + k3[1], u_z[1])*dz
        Y_z_1 = Y_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
        Q_z_1 = Q_z_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
        return np.asarray([Y_z_1,Q_z_1])

    def f_der_Y_z(Y_z_n,Q_z_n,u_n):
        a = Q_z_n
        b = np.exp(Y_z_n) + 1/u_n - (-n_d*m_i/u_n)/(e_charge*n_s)
        #b = np.exp(Y_z_n) + 1/u_n
        return np.asarray([a,b])

    #MSICALEANEOUS
    def produce_Q_DC(Y_z):
        Q = []
        for v in np.arange(len(z_list)-1):
            Q.append((Y_z[v+1]-Y_z[v])/dz)
        Q.append(alpha)
        return np.asarray(Q)

    def produce_Y_DC(Y):
        Y_DC = []
        for v in np.arange(len(z_list)):
            a = []
            for i in np.arange(len(T_list)):
                a.append(Y[i][v])
            Y_DC.append(np.mean(a))
        return np.asarray(Y_DC)
    
    def produce_Y_difference(Y_bar_temp, Y_bar): 
        a = Y_bar_temp -  Y_bar
        return np.linalg.norm(a)/np.linalg.norm(Y_bar)

    #REPEAT THAT SHIT
    def produce_Y_DC_converge(Y_DC_z_init,Q_init,u_init):
        Q_row = produce_Q_DC(Y_DC_z_init)
        u = produce_u(Q_row,u_init)
        u_init = u[0]
        a = produce_Y_z_t(Y_DC_z_init[0],Q_init,u)
        Y_z_t_converge = a[0]
        Q_z_t_converge = a[1]
        Q_init = Q_z_t_converge[:,0]
        Y_DC_converge = produce_Y_DC(Y_z_t_converge)
        Q_DC_converge = produce_Y_DC(Q_z_t_converge)
        Y_difference = produce_Y_difference(0,Y_DC_converge)
        reps = 0
        while(Y_difference > Y_DC_epsilon):
            reps+=1
            Y_DC_temp = np.copy(Y_DC_converge)
            u = produce_u(Q_DC_converge,u_init)
            u_init = u[0]
            a = produce_Y_z_t(Y_DC_converge[0],Q_init,u)
            Y_z_t_converge = a[0]
            Q_z_t_converge = a[1]
            Q_init = Q_z_t_converge[:,0]
            Y_DC_converge = produce_Y_DC(Y_z_t_converge)
            Q_DC_converge = produce_Y_DC(Q_z_t_converge)
            Y_difference = produce_Y_difference(Y_DC_temp,Y_DC_converge)
            #print("Y_difference = ", Y_difference)
        return np.asarray([Y_z_t_converge,Y_DC_converge,Q_DC_converge])

    z_list = produce_z_vals()
    T_list = produce_T_vals()
    V_DC = produce_V_DC(V_RF)
    Y_z_init = produce_Y_z_init(V_DC)
    Q_init = np.asarray([0.01]*len(T_list))
    u_init = -1
    results = produce_Y_DC_converge(Y_z_init,Q_init,u_init)
    Y_z_t = results[0]
    Y_DC_z = results[1]
    Q_final = results[2]

    def produce_electric_field(Y,Q_init_val):
        Q = [Q_init_val]
        for v in np.arange(1,len(z_list) - 1):
            Q.append((Y[v+1] - Y[v-1])/(2*dz))
        Q.append(Q[-1])
        return Q

    E = []
    for i in np.arange(len(T_list)):
        E.append(produce_electric_field(Y_z_t[i],Q_final[i]))
    E_average = produce_Y_DC(E)

    plt.figure()
    plt.grid()
    plt.title("Electric Potential: RF_collisional")
    plt.xlabel("z")
    plt.ylabel("Y")
    for i in np.arange(len(T_list)):
        plt.plot(np.asarray(z_list),Y_z_t[i] , label = "t = " + str(T_list[i]))
    plt.plot(np.asarray(z_list),Y_DC_z, label = "time averaged" )
    plt.plot(np.asarray(z_list),Y_z_init, label = "original")
    plt.legend()
    plt.savefig("Figures/Electric_potential_vs_z_RF_collisional.png")

    plt.figure()
    plt.grid()
    plt.title("Electric Field: RF_collisional")
    plt.xlabel("z")
    plt.ylabel("E")
    for i in np.arange(len(T_list)):
        plt.plot(np.asarray(z_list),E[i] , label = "t = " + str(T_list[i]))
    plt.plot(np.asarray(z_list),E_average, label = "time averaged" )
    plt.legend()
    plt.savefig("Figures/Electric_Field_vs_z_RF_collisional.png")
    return [np.asarray(z_list),Y_DC_z]
a = collisionless()
b = collisional()

plt.figure()
plt.grid()
plt.title("Electric Potential")
plt.xlabel("z")
plt.ylabel("Y")
plt.plot(a[0],a[1], label = "collisionless" )
plt.plot(b[0],b[1], label = "collisional" )
plt.legend()
plt.savefig("Figures/Electric_potential.png")

plt.show()