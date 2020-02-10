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
drop_height = 9.8*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 25.0*lambda_D#set radius of contianer ie wall radius
coulomb_limit = 5

z_se = drop_height
V_DC = -100
k_z_restore = -2.0*V_DC/((z_se/lambda_D)**2)
Y_DC = (e_charge*V_DC)/(k_b*T_e)
#Y_DC = -50.0
#related to finding the charge of the dust grains
C_1 = 0.5*np.log(2*np.pi*m_e/m_i) - Y_DC
C_2 = np.exp(C_1)
def root_find_bessel(x):
    return iv(0,x) - C_2
result = minimize_scalar(root_find_bessel,bounds=(0, 200), method='bounded')
Y_ref = result.x
print("Y_DC = ", Y_DC)
print("Y_ref = ", Y_ref)

T = 1/(13.6*1e6)
dT = T/3
dz = lambda_D/2
omega = (2*np.pi)/T#TAKEN FROM NITTER PAPER 13.6M Hz
phi_wall_z = Y_DC #TRY THIS FOR NOW SEE HOW IT WORKS
n_s = n_i0*np.exp(0.5)
n_d = (dust_grain_max_input)/(np.pi*pow(container_radius,2)*drop_height)
root = dz/10#preciscion of root finding method used to get dust charge
Y_epsilon = 1e-3#THIS HAS ALMOST NO EFFECT WHY??
low_z = 0
high_z = drop_height
sigma = 3e-19
alpha = ((epsilon_0*k_b*T_e)/(n_i0*np.exp(0.5)*(e_charge**2)))**0.5*n_n0*sigma
print("alpha = ", alpha)
plots = "yes"
anim = "no"

#%%
#run_list = [["DC","collisionless"],["DC","collisional"]]
run_list = [["DC","collisionless"],["RF","collisionless"],["DC","collisional"],["RF","collisional"]]

def produce_z_vals(dz,low_z,high_z):
    z_list = []
    for z in np.arange(low_z,high_z,dz):
        z_list.append(z/lambda_D)
    return z_list

def produce_T_vals(d_T,T):
    T_list = []
    for t in np.arange(0,T,dT):
        T_list.append(t)
    #print("T_list = ",T_list)
    return T_list

#%%
def run_the_thing(sheath,collision):

    #CALC Y DISTRIBUTION USING NITTER 
    def produce_Y_0_DC():
        Y_0 = []
        for v in np.arange(len(z_list)):
            #d = -(e_charge/(k_b*T_e))*(k_z_restore/2)*(z_list[v] - z_se/lambda_D)**2
            d = Y_DC
            Y_0.append(d)
        return np.asarray(Y_0) 

    # Y is an array of discretised Y(z) arrays at various t
    def produce_Y_0_RF():
        Y_0 = []
        for i in np.arange(len(T_list)):
            Y_0_row = []
            for v in np.arange(len(z_list)-1):
                #d = -(e_charge/(k_b*T_e))*(k_z_restore/2)*(z_list[v] - z_se/lambda_D)**2
                d = Y_DC
                Y_0_row.append(d + Y_ref*np.sin(omega*T_list[i]))
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
                Q_row.append( (Y_t_z[i][v+1]-Y_t_z[i][v])/h  )
            if collision == "collisionless":
                Q_row.append(0)
            else:
                Q_row.append(alpha)
            Q.append(Q_row)
        return np.asarray(Q)

    def produce_u_0_collisionless(Y_z):
        u_0 = []
        for v in np.arange(len(z_list)-1):
            u_0.append(-(1-2*Y_z[v])**0.5)
        u_0.append(-1)
        return np.asarray(u_0)

    def f_der_u(u_0,Q):
        return -(1/u_0)*Q + alpha*u_0

    def step_u(u_0,Q,h):
        k1 = f_der_u(u_0,Q)*h
        k2 = f_der_u(u_0 + k1/2,Q)*h
        k3 = f_der_u(u_0 + k2/2,Q)*h
        k4 = f_der_u(u_0 + k3,Q)*h
        u_1 = u_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0
        return u_1

    def produce_u_0_collisional(Q_row,h):
        u = [-1]
        for v in np.arange(len(z_list) - 1):
            u.append(step_u(u[v],Q_row[v],-h))
        return np.flip(u)#its done backwards so need to flip it

    def f_der(Y_0, Q_0, u_z):
        p_d = -(n_d*m_D)/u_z
        print(p_d)
        a = np.asarray([Q_0, np.exp(Y_0) + 1/u_z - p_d/(e_charge*n_s)])
        print(np.exp(Y_0), 1/u_z, -p_d/(e_charge*n_s))
        return a

    def step(Y_0, Q_0, u_z, h):
        print(collision,sheath)
        k1 = f_der(Y_0,Q_0,u_z)*h
        k2 = f_der(Y_0 + k1[0]*0.5, Q_0 + k1[1]*0.5,u_z)*h
        k3 = f_der(Y_0 + k2[0]*0.5, Q_0 + k2[1]*0.5,u_z)*h
        k4 = f_der(Y_0 + k3[0], Q_0 + k3[1],u_z)*h
        a = Y_0 + (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6.0
        b = Q_0 + (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6.0
        return np.asarray([a,b])

    def produce_Y_Q_new_DC(Y_z,Q_z,u_z,h):
        Y_z = [0]
        if collision == "collisionless":
            Q_z = [0]
        else:
            Q_z = [alpha]
        for v in np.arange(len(z_list)-1):
            a = step(Y_z[v], Q_z[v], u_z[v], -h)
            Y_z.append(a[0])
            Q_z.append(a[1])
        return np.asarray([np.flip(Y_z),np.flip(Q_z)])

    def produce_Y_Q_new_RF(Y_z_t,Q_z_t,u,h):
        Y_new = []
        Q_new = []
        #print(Y_z_t[0][5],Y_z_t[2][5])
        for i in np.arange(len(T_list)):
            b = produce_Y_Q_new_DC(Y_z_t[i],Q_z_t[i],u,h)
            Y_new.append(b[0])
            Q_new.append(b[1])
        return np.asarray([Y_new,Q_new])

    def produce_Y_difference(Y_bar_temp, Y_bar): 
        a = Y_bar_temp -  Y_bar
        return np.linalg.norm(a)/np.linalg.norm(Y_bar)

    def find_Y(h):
        print(collision,sheath)
        if collision == "collisionless":
            if sheath == "DC":
                Y = produce_Y_0_DC()
                Q = produce_Q_DC(Y,h)
            else:
                Y = produce_Y_0_RF()
                Q = produce_Q_RF(Y,h)
                Y_bar = produce_Y_bar(Y)

            Y_difference = 2*Y_epsilon
            reps = 0
            while(Y_difference > Y_epsilon):
                reps+=1
                if sheath == "DC":
                    Y_bar_temp = np.copy(Y)
                    u = produce_u_0_collisionless(Y)
                    a = produce_Y_Q_new_DC(Y,Q,u,h)
                    Y = a[0]
                    Q = a[1]
                    Y_difference = produce_Y_difference(Y_bar_temp,Y)
                    print("Y_diff = ",Y_difference)
                else:
                    Y_bar_temp = np.copy(Y_bar)
                    u = produce_u_0_collisionless(Y_bar)
                    a = produce_Y_Q_new_RF(Y,Q,u,h)
                    Y = a[0]
                    Q = a[1]
                    Y_bar = produce_Y_bar(Y)
                    Y_difference = produce_Y_difference(Y_bar_temp,Y_bar)
                    print("Y_diff = ",Y_difference)
            print("reps = ", reps)
        else:
            if sheath == "DC":
                Y = produce_Y_0_DC()
                Q = produce_Q_DC(Y,h)
            else:
                Y = produce_Y_0_RF()
                Q = produce_Q_RF(Y,h)
                Y_bar = produce_Y_bar(Y)

            Y_difference = 2*Y_epsilon
            reps = 0
            while(Y_difference > Y_epsilon):
                reps+=1
                if sheath == "DC":
                    Y_bar_temp = np.copy(Y)
                    u = produce_u_0_collisional(Q,h)
                    a = produce_Y_Q_new_DC(Y,Q,u,h)
                    Y = a[0]
                    Q = a[1]
                    Y_difference = produce_Y_difference(Y_bar_temp,Y)
                    print("Y_diff = ",Y_difference)
                else:
                    Y_bar_temp = np.copy(Y_bar)
                    Q_bar = produce_Y_bar(Q)
                    u = produce_u_0_collisional(Q_bar,h)
                    a = produce_Y_Q_new_RF(Y,Q,u,h)
                    Y = a[0]
                    Q = a[1]
                    Y_bar = produce_Y_bar(Y)
                    Y_difference = produce_Y_difference(Y_bar_temp,Y_bar)
                    print("Y_diff = ",Y_difference)
            print("reps = ", reps)
        return Y
    
    z_list = produce_z_vals(dz,low_z,high_z)
    T_list = produce_T_vals(dT,T)
    data = find_Y(dz)

    return data

#%%
data_plot = []
for l in run_list:
    data_plot.append(run_the_thing(l[0],l[1]))
z_vals = produce_z_vals(dz,low_z,high_z)
T_vals = produce_T_vals(dT,T)

#print(len(data_plot),len(data_plot[0]), len(data_plot[1]))
print("done")

# print("hey")
# print(data_plot[0])
# print("hey")
# print(data_plot[1])
# print("hey")
# print(data_plot[2])
# print("hey")
# print(data_plot[3])

def produce_Y_bar_out(Y):
    Y_bar = []
    for v in np.arange(len(z_vals)):
        a = Y[:,v]
        Y_bar.append(np.mean(a))
    return np.asarray(Y_bar)

def produce_Y_0_DC_out():
    Y_0 = []
    for v in np.arange(len(z_vals)):
        #d = -(k_z_restore/2)*(z_vals[v] - z_se/lambda_D)**2
        d = Y_DC
        Y_0.append(d)
    return np.asarray(Y_0) 

original_data = produce_Y_0_DC_out()
time_av_RF_collisionless = produce_Y_bar_out(data_plot[1])
time_av_RF_collisional = produce_Y_bar_out(data_plot[3])

if plots == "yes":
    plt.figure()
    plt.grid()
    plt.title("Electric Potential:DC_collisionless")
    plt.xlabel("z")
    plt.ylabel("Y")
    plt.plot(np.asarray(z_vals),data_plot[0])
    plt.savefig("Figures/Electric_potential_vs_z_DC_collisionless.png")

    plt.figure()
    plt.grid()
    plt.title("Electric Potential:RF_collisionless")
    plt.xlabel("z")
    plt.ylabel("Y")
    for i in np.arange(len(T_vals)):
        plt.plot(np.asarray(z_vals),data_plot[1][i] , label = "t = " + str(T_vals[i]))
    plt.plot(np.asarray(z_vals),time_av_RF_collisionless, label = "time averaged" )
    plt.legend()
    plt.savefig("Figures/Electric_potential_vs_z_RF_collisionless.png")

    plt.figure()
    plt.grid()
    plt.title("Electric Potential:DC_collisional")
    plt.xlabel("z")
    plt.ylabel("Y")
    #plt.ylim(-100,100)
    plt.plot(np.asarray(z_vals),data_plot[2])
    plt.savefig("Figures/Electric_potential_vs_z_DC_collisional.png")

    plt.figure()
    plt.grid()
    plt.title("Electric Potential:RF_collisional")
    plt.xlabel("z")
    plt.ylabel("Y")
    for i in np.arange(len(T_vals)):
        plt.plot(np.asarray(z_vals),data_plot[3][i] , label = "t = " + str(T_vals[i]))
    plt.plot(np.asarray(z_vals),time_av_RF_collisional, label = "time averaged" )
    plt.legend()
    plt.savefig("Figures/Electric_potential_vs_z_RF_collisional.png")

    plt.figure()
    plt.grid()
    plt.title("Electric Potentials")
    plt.xlabel("z")
    plt.ylabel("Y")
    plt.plot(np.asarray(z_vals),data_plot[0], label = "time averaged DC collisionless" )
    plt.plot(np.asarray(z_vals),time_av_RF_collisionless, label = "time averaged RF collisionless")
    plt.plot(np.asarray(z_vals),data_plot[2], label = "time averaged DC collisional"  )
    plt.plot(np.asarray(z_vals),time_av_RF_collisional, label = "time averaged RF collisional"  )
    #plt.plot(np.asarray(z_vals)/lambda_D,original_data, label = "original input DC" )
    plt.legend()
    #plt.ylim(-100,100)
    plt.savefig("Figures/Electric_potential_vs_z_all.png")

if anim == "yes":
    speed_mul = 100

    fig, ax = plt.subplots()
    ax.set_xlabel("z/lambda_D")
    ax.set_ylabel("Y")
    ax.set_title("Animation")
    ax.set_xlim(0,np.asarray(z_vals[-1])/lambda_D)
    ax.set_ylim(-100,100)
    line_1, = ax.plot([], [], lw=2)
    line_2, = ax.plot([], [], lw=2)
    plt.grid()
    text = plt.text(z_vals[0]/lambda_D + 0.2, 100 - 5,"")

    def init():
        line_1.set_data([], [])
        line_2.set_data([], [])
        text.set_text("Time = " + str(0) + "s")
        return [line_1,line_2,text]

    def animate(i):
        line_1.set_data(np.asarray(z_vals)/lambda_D, data_plot[1][i] )  # update the data.
        line_2.set_data(np.asarray(z_vals)/lambda_D, data_plot[3][i] )  # update the data.
        text.set_text("Time = " + str(T_vals[i]) + "s")#updat value of the frame number
        return [line_1,line_2,text]

    ani = animation.FuncAnimation(fig, animate,init_func=init,frames = len(T_vals), interval=100, blit=True)

plt.show()
        
    
