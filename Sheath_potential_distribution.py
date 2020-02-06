import numpy as np
import matplotlib.pyplot as plt#import module used to produce graphs

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
dust_grain_density = 1.49*1e3 #dust density
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
dz = lambda_D/100
omega = (2*np.pi)/T#TAKEN FROM NITTER PAPER 13.6M Hz
Y_DC = -50.0#NITTER
phi_wall_z = Y_DC #TRY THIS FOR NOW SEE HOW IT WORKS
Y_ref = 50.0#NITTER
n_s = n_i0*np.exp(0.5)
p_d = (dust_grain_max_input*m_D)/(np.pi*pow(container_radius,2)*drop_height)
root = dz/10#preciscion of root finding method used to get dust charge
Y_epsilon = 1e-12#VARY THIS

low_z = 7.5*lambda_D
high_z = 20.0*lambda_D

z_list = []
T_list = []
time_potential_list = []

#%%
    
#CALC Y DISTRIBUTION USING NITTER 
def produce_z_vals(dz,low_z,high_z):
    for z in np.arange(low_z,high_z,dz):
        z_list.append(z)

def produce_T_vals(d_T,T):
    for t in np.arange(0,T,dT):
        T_list.append(t)

# Y is an array of discretised Y(z) arrays at various t
def produce_Y_0():
    #rows are for a given pos z and inside is for a range of t values
    Y_0 = []
    for i in np.arange(len(T_list)):
        Y_0_row = []
        for v in np.arange(len(z_list)):
            Y_0_row.append(Y_DC + Y_ref*np.sin(omega*T_list[i]) )
        Y_0.append(Y_0_row)
    #print(len(Y_0), len(Y_0[0]))
    return np.asarray(Y_0) 

def produce_u_0():
    u_0 = []
    for v in np.arange(len(z_list)):
        u_0.append(-1)
    return np.asarray(u_0)

# Y_bar is time-averaged Y(z,t) over time period
def produce_Y_bar(Y):
    Y_bar = []
    #print("ybar")
    #print(len(Y), len(Y[0]))
    for v in np.arange(len(z_list)):
        #print(v)
        a = Y[:,v]
        #print(len(a))
        Y_bar.append(np.mean(a))
    return np.asarray(Y_bar)

def f_der_u_new(Y_bar, u_0,  h):
    f = []
    f.append((-1/u_0[0])*(Y_bar[1] - Y_bar[0])/(2*h))
    #print(len(Y_bar), len(u_0))
    for v in np.arange(1, len(z_list)-1):
        #print(v)
        f.append((-1/u_0[v])*(Y_bar[v+1] - Y_bar[v-1])/(2*h))
    f.append(-1)
    return np.asarray(f)

def step_u_new(Y_bar, u_0,  h):
    k1 = f_der_u_new(Y_bar,u_0,h)*h
    k2 = f_der_u_new(Y_bar + h/2,u_0 + k1*0.5,h)*h
    k3 = f_der_u_new(Y_bar + h/2,u_0 + k2*0.5,h)*h
    k4 = f_der_u_new(Y_bar + h,u_0 + k3,h)*h
    return u_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0

def produce_Q(Y_row,  h):
    Q = []
    Q.append((Y_row[2] -2*Y_row[1] + Y_row[0])/pow(h,2))
    for v in np.arange(1,len(z_list)-1):
        Q.append((Y_row[v+1] -2*Y_row[v] + Y_row[v-1])/pow(h,2))
    Q.append((Y_row[-1] -2*Y_row[-2] + Y_row[-3])/pow(h,2))
    return np.asarray(Q)

def f_der_Y_new(Y_row, Q, u,  h):
    f = []
    dY_dz = []
    dQ_dz = []
    for v in np.arange(len(z_list)):#NOT SURE ABOUT LIMITS???
        dY_dz.append(Q[v])
        dQ_dz.append(np.exp(Y_row[v]) + 1/u[v] - p_d/(e_charge*n_s))
    
    f.append(dY_dz)
    f.append(dQ_dz)
    #print("hii")
    #print(len(f), len(f[0]))
    return np.asarray(f)

def step_Y_new(Y_row, u,  h):
    Q = produce_Q(Y_row,h)#CURRETNLY RECALCULATING FOR NO FUCKING REASON
    k1 = f_der_Y_new( Y_row,Q,u,h)*h
    k2 = f_der_Y_new( Y_row + k1[0]*0.5, Q+ k1[1]*0.5,u + h/2,h)*h#SHOULD THESE BE h or h/2???????
    k3 = f_der_Y_new( Y_row + k2[0]*0.5, Q+ k2[1]*0.5,u + h/2,h)*h
    k4 = f_der_Y_new( Y_row + k3[0], Q+k3[1],u + h,h)*h
    a = Y_row + (k1 + 2*k2 + 2*k3 + k4)/6.0
    return a[0]

def produce_Y_new(Y, u,  h):
    Y_new = []
    #print("hee")
    #print(len(Y))
    #print(len(Y[0]))
    for i in np.arange(len(T_list)):
        Y_new.append(step_Y_new( Y[i] ,u,h) )
    #print("asdasdsa")
    #print(len(Y_new))
    #print(len(Y_new[0]))
    return np.asarray(Y_new)

def produce_Y_difference(Y_bar_temp, Y_bar):  
    return np.linalg.norm(Y_bar_temp -  Y_bar)

def find_Y(h):
    Y = produce_Y_0()
    u = produce_u_0()
    Y_bar = produce_Y_bar(Y)
    Y_difference = 2*Y_epsilon
    while(Y_difference > Y_epsilon):
        Y_bar_temp = Y_bar
        u = step_u_new(Y_bar,u,h)
        Y = produce_Y_new(Y,u,h)
        Y_bar = produce_Y_bar(Y)
        Y_difference = produce_Y_difference(Y_bar_temp,Y_bar)
        print(Y_difference)
    return Y

produce_z_vals(dz,low_z,high_z)
produce_T_vals(dT,T)
data = find_Y(dz)

plt.figure()
plt.grid()
plt.title("Electric Potential")
plt.xlabel("z")
plt.ylabel("Y")
for i in np.arange(len(data)):
    plt.plot(z_list,data[i] , label = "t = " + str(T_list[i]) )
plt.savefig("Figures/Electric_potential_vs_z.png")
plt.show()
        
    
