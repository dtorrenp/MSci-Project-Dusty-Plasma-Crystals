import numpy as np
import sys
import matplotlib.pyplot as plt#import module used to produce graphs
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
drop_height = 9.8*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 25.0*lambda_D#set radius of contianer ie wall radius
coulomb_limit = 5

#related to finding the charge of the dust grains

T = 1/(13.6*1e6)
dT = T/10
dz = lambda_D/10
omega = (2*np.pi)/T#TAKEN FROM NITTER PAPER 13.6M Hz
Y_DC = -50.0#NITTER
phi_wall_z = Y_DC #TRY THIS FOR NOW SEE HOW IT WORKS
Y_ref = 48.0#NITTER
n_s = n_i0*np.exp(0.5)
p_d = (dust_grain_max_input*m_D)/(np.pi*pow(container_radius,2)*drop_height)
root = dz/10#preciscion of root finding method used to get dust charge
Y_epsilon = 1e-3#VARY THIS

low_z = 0
high_z = 12*lambda_D
sigma = 3e-19
alpha = ((epsilon_0*k_b*T_e)/(n_i0*np.exp(0.5)*(e_charge**2)))**0.5*n_n0*sigma
#alpha_them = ((epsilon_0*k_b*T_e)/(6.1e15*np.exp(0.5)*(e_charge**2)))**0.5*3.3e21*sigma
#print(alpha,alpha_them)
#%%
sheath = "RF" #RF
collision = "collisionless" #collisionless

#%%

z_list = []
T_list = []
time_potential_list = []
potential_list = []


#%%

#CALC Y DISTRIBUTION USING NITTER 
def produce_z_vals(dz,low_z,high_z):
    for z in np.arange(low_z,high_z,dz):
        z_list.append(z)

def produce_T_vals(d_T,T):
    for t in np.arange(0,T,dT):
        T_list.append(t)

def produce_Y_0_DC():
    Y_0 = []
    for v in np.arange(len(z_list)):
        Y_0.append(Y_DC)
    return np.asarray(Y_0) 

# Y is an array of discretised Y(z) arrays at various t
def produce_Y_0_RF():
    Y_0 = []
    for i in np.arange(len(T_list)):
        Y_0_row = []
        for v in np.arange(len(z_list)-1):
            Y_0_row.append(Y_DC + Y_ref*np.sin(omega*T_list[i]))
        Y_0_row.append(0) #boundary condition
        Y_0.append(Y_0_row)
    return np.asarray(Y_0) 

# Y_bar is time-averaged Y(z,t) over time period
def produce_Y_bar(Y):
    Y_bar = []
    for v in np.arange(len(z_list)):
        a = Y[:,v]
        Y_bar.append(np.mean(a))
    return np.asarray(Y_bar)

def produce_Q_DC(Y_z,h):
    Q = []
    for v in np.arange(len(z_list)-1):
        Q.append((Y_z[v+1]-Y_z[v])/h)
    if collision == "collisionless":
        Q.append(0)
    else:
        Q.append(alpha)
    return np.asarray(Q)

def produce_Q_RF(Y_t_z,h):
    Q = []
    for i in np.arange(len(T_list)):
        Q_row = []
        for v in np.arange(len(z_list)-1):
            Q_row.append((Y_t_z[i][v+1]-Y_t_z[i][v])/h)
        if collision == "collisionless":
            Q_row.append(0)
        else:
            Q_row.append(alpha)
        Q.append(Q_row)
    print ("HEllo", len(Q), len(Q[0]))
    return np.asarray(Q)

def produce_u_0_collisionless(Y_z):
    u_0 = []
    #print ("hey")
    #rint (Y_z)
    for v in np.arange(len(z_list)-1):
        u_0.append(-(1-2*Y_z[v])**0.5)
    u_0.append(-1)
    #print (u_0)
    #print (len(u_0))
    return np.asarray(u_0)

def produce_u_0_collisional(Q):
    #print(Q)
    u_0 = []
    for v in np.arange(len(z_list)-1):
        #print("hiiii")
        #print(np.exp(-2*alpha*v))
        #print( (1/alpha)*(Q[v]))
        u_0.append((np.exp(-2*alpha*v) + (1/alpha)*(Q[v]))**0.5)
    u_0.append(-1)
    return np.asarray(u_0)

def f_der(Y_row, Q, u,  h):
    f = []
    dY_dz = []
    dQ_dz = []
    #print (u)
    for v in np.arange(len(z_list)):#NOT SURE ABOUT LIMITS???
        dY_dz.append(Q[v])
        #print (np.exp(Y_row[v]), 1/u[v], p_d/(e_charge*n_s))
        dQ_dz.append(np.exp(Y_row[v]) + 1/u[v] - p_d/(e_charge*n_s))
    f.append(dY_dz)
    f.append(dQ_dz)
    return np.asarray(f)

def step(Y_row, Q_row, u,  h):
    k1 = f_der(Y_row,Q_row,u,h)*h
    #print (k1, Y_row, Q_row, u)
    #print (len(Q_row), len(k1[0]))
    k2 = f_der(Y_row + k1[0]*0.5, Q_row + k1[1]*0.5,u + h/2,h)*h
    k3 = f_der( Y_row + k2[0]*0.5, Q_row+ k2[1]*0.5,u + h/2,h)*h
    k4 = f_der( Y_row + k3[0], Q_row+k3[1],u + h,h)*h
    a = np.asarray([Y_row,Q_row]) + (k1 + 2*k2 + 2*k3 + k4)/6.0
    return a

def produce_Y_Q_new_DC(Y_z,Q_z,u,h):
    return np.asarray(step(Y_z,Q_z,u,h))

def produce_Y_Q_new_RF(Y_z_t,Q_z_t,u,h):
    Y_new = []
    Q_new = []
    #print (len(Y_z_t), len(Q_z_t), len(Y_z_t[0]), len(Q_z_t[0]))
    for i in np.arange(len(T_list)):
        b = step(Y_z_t[i],Q_z_t[i],u,h)
        Y_new.append(b[0])
        Q_new.append(b[1])
    return np.asarray([Y_new,Q_new])

def produce_Y_difference(Y_bar_temp, Y_bar):  
    return np.linalg.norm(Y_bar_temp -  Y_bar)

def find_Y(h):
    if collision == "collisionless":
        if sheath == "DC":
            Y = produce_Y_0_DC()
            Q = produce_Q_DC(Y,h)
        else:
            Y = produce_Y_0_RF()
            #print(Y)
            Q = produce_Q_RF(Y,h)
            Y_bar = produce_Y_bar(Y)

        Y_difference = 2*Y_epsilon
        while(Y_difference > Y_epsilon):
            #print("hey")
            if sheath == "DC":
                Y_bar_temp = Y
                u = produce_u_0_collisionless(Y)
                #print (len(Q))
                a = produce_Y_Q_new_DC(Y,Q,u,h)
                Y = a[0]
                Q = a[1]
                Y_difference = produce_Y_difference(Y_bar_temp,Y)
            else:
                Y_bar_temp = Y_bar
                u = produce_u_0_collisionless(Y_bar)
                a = produce_Y_Q_new_RF(Y,Q,u,h)
                Y = a[0]
                Q = a[1]
                Y_bar = produce_Y_bar(Y)
                Y_difference = produce_Y_difference(Y_bar_temp,Y_bar)
    else:
        if sheath == "DC":
            Y = produce_Y_0_DC()
            Q = produce_Q_DC(Y,h)
        else:
            Y = produce_Y_0_RF()
            Q = produce_Q_RF(Y,h)
            Y_bar = produce_Y_bar(Y)
        Y_difference = 2*Y_epsilon

        while(Y_difference > Y_epsilon):
            print("hey")
            if sheath == "DC":
                Y_bar_temp = Y
                u = produce_u_0_collisional(Q)
                a = produce_Y_Q_new_DC(Y,Q,u,h)
                Y = a[0]
                Q = a[1]
                Y_difference = produce_Y_difference(Y_bar_temp,Y)
            else:
                Y_bar_temp = Y_bar
                Q_bar = produce_Y_bar(Q)
                u = produce_u_0_collisional(Q_bar)
                a = produce_Y_Q_new_RF(Y,Q,u,h)
                Y = a[0]
                Q = a[1]
                Y_bar = produce_Y_bar(Y)
                Y_difference = produce_Y_difference(Y_bar_temp,Y_bar)
    return Y

#%%
produce_z_vals(dz,low_z,high_z)
produce_T_vals(dT,T)
data = find_Y(dz)
print("done")
#print(data)

plt.figure()
plt.grid()
plt.title("Electric Potential:" + str(collision) + " " + str(sheath) )
plt.xlabel("z")
plt.ylabel("Y")
if sheath == "DC":
    plt.plot(np.asarray(z_list)/lambda_D,data)
else:
    for i in np.arange(len(data)):
        plt.plot(np.asarray(z_list)/lambda_D,data[i] , label = "t = " + str(T_list[i]))
plt.savefig("Figures/Electric_potential_vs_z_" + str(collision) + "_" + str(sheath) + ".png")
plt.show()
        
    
