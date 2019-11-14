#%%

"""Import the jazz"""
import numpy as np
import time
import random
import matplotlib.pyplot as plt#import module used to produce graphs
from mpl_toolkits.mplot3d import Axes3D
import itertools as it
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay
#plt.ioff()
#%%
"""Import stuff to measure how long the code takes to run"""
start_time = time.time()
print("start_time =", time.ctime(time.time()))
#%%
#random.seed(30)

#RK45 Parameters
#a_2 = 1/5
#a_3 = 3/10
#a_4 = 3/5
#a_5 = 1
#a_6 = 7/8
#
#c_1 = 37/378
#c_2 = 0
#c_3 = 250/621
#c_4 = 125/594
#c_5 = 0
#c_6 = 512/1771
#
#c_1_s = 2825/27648
#c_2_s = 0
#c_3_s = 18575/48384
#c_4_s = 13525/55296
#c_5_s = 277/14336
#c_6_s = 1/4
#
#b_21 = 1/5
#b_31 = 3/40
#b_32 = 9/40
#b_41 = 3/10
#b_42 = -9/10
#b_43 = 6/5
#b_51 = -11/54
#b_52 = 5/2
#b_53 = -70/27
#b_54 = 35/27
#b_61 = 1631/55296
#b_62 = 175/512
#b_63 = 575/13824
#b_64 = 44275/110592
#b_65 = 253/4096

#%%
#DEFINE PLASMA DISCHARGE CONDITIONS

n_e0 = 1e15#electron and ion densities in bulk plasma
n_i0 = 1e15

g_z = 9.81#gravity
e_charge = -1.6*1e-19
i_charge = 1.6*1e-19
grain_R = 7*1e-6
m_i = 1.67*1e-27
m_e = 9.11*1e-31
m_D = ((4/3)*np.pi*grain_R**3)*(1.49*1e3)#mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
epsilon_0 = 8.85*1e-12
k_b = 1.38*1e-23
mu = (m_i/m_e)#normalization used for convienience

beta = 1.5*1e-2#T_i/T_e
T_e = (2*1.6*1e-19)/k_b
print("T_e = ",T_e)
T_i = T_e*beta#defined a bit wierdly but basically just in terms of beta and T_e
print("T_i = ",T_i)
v_i = (2*k_b*T_i/m_i)**(1/2)
Z = 1#atomic number??
a_0 = 1#intial guess for halley's method

lambda_de = ((epsilon_0*k_b*T_e)/(n_e0*(e_charge**2)))**0.5
lambda_di = ((epsilon_0*k_b*T_i)/(n_i0*(e_charge**2)))**0.5
lambda_D = (1/(1/(lambda_de**2) + 1/(lambda_di**2)))**0.5#get lambda_D
print("lambda_D =",lambda_D )

container_height = 11*lambda_D#drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
container_radius = 100*lambda_D#set radius of contianer ie wall radius
z_se = 10*lambda_D#distance from bottom of container to the sheath edge
r_se = 100*lambda_D#distance from wall to the sheathe edge
r_se_inv = container_radius - r_se

wake_potential_below = 0.5*lambda_D
wake_charge_multiplier = 1e-3

phi_wall_z = -100 #volts
phi_wall_r = -1 #volts
k_z_restore = -2*phi_wall_z/z_se**2#WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
k_r_restore = -2*phi_wall_r/r_se**2

v_B = (k_b*T_e/m_i)**0.5
print("v_B =",v_B)

alpha = 1e-9#drag coefficient
time_list = [0]
root = 1e-14#defined preciscion of root finding method used to get dust charge
#S = 0.95

dt = 1e-4#time step in rk4, needs to be small enough to be precise but large enough we can actually move the stuff forward in time
#eps = 1e-7
dust_grain_max = 5 #dust grain max number
print("grain number = ", dust_grain_max)
frames = 1e5 #number of frames, time taken is not linear as teh longer u run it the more particles it adds hence increases quadratically
print("frames = ",frames)
temp_min = 1e-5#minimum temperature for it to stop
for_run = True 
print("for_run =",for_run)

#%%

class Dust_Container:
    """Define container dust contianer, this includes all properties independant of dust particles,method to find the surface potetetial of dust grains and list of dust grains"""
    def __init__(self,container_radius, container_height, m_D, grain_R):
            self.container_radius = container_radius
            self.container_height = container_height
            self.m_D = m_D
            self.grain_R = grain_R
            self.init_W_vec = np.array([0,0,self.container_height,0,0,0])#its dumb but need inital vector for first dust grain W_vec = [x,y,z,vx,vy,vz]
            self.phi_grain = self.OML_find_surface_potential()#defined surface potential of the dust grains, this is a properties of the plasma not the dust grain(however the dust grain charge is a property of the dust grain as it requires radius of dust grain)
            self.Dust_grain_list = np.array([Dust_grain(self.m_D , self.grain_R , self.produce_q_D() , self.prod_W_vec(), time_list[-1])])#create the first dust grain
            self.v_squared_sum = (self.Dust_grain_list[0].calc_speed())**2
            self.comb_list = list(it.combinations(self.Dust_grain_list,2))
            self.temperature = 0
    
    def OML_find_surface_potential(self):
        """Find the dust grain sufrace potential using OML model"""
        def find_W_H(a_0,z,root_prec):
            """root finding method using "halley's method" like newtons method but third order """
            def f_x(x_n,z):
                return x_n*np.exp(x_n) - z
            
            def f_x_first_derv(x_n):
                return np.exp(x_n)*(1 + x_n)
            
            def f_x_second_derv(x_n):
                return np.exp(x_n)*(2 + x_n)
            
            def calc_x_plus(x_n,z):
                f_x_0 = f_x(x_n,z)
                f_x_1 = f_x_first_derv(x_n)
                f_x_2 = f_x_second_derv(x_n)
                x_plus = x_n - ((2*f_x_0*f_x_1)/(2*(f_x_1**2)-f_x_0*f_x_2))
                return x_plus
        
            a_n = a_0
            a_plus = calc_x_plus(a_0,z)
        
            while(np.abs(a_plus - a_n) > root_prec):
                a_n = a_plus
                a_plus = calc_x_plus(a_n,z)
            
            W_0 = (a_n + a_plus)/2
            return W_0

        """"basically we have an equation that has the normalized surface potetnail in it and we rootfind for it, using the lambert W function form though tbh i agree with coppins this isnt actually necesary"""
        k = (((mu*beta)**0.5)/Z)*np.exp(beta/Z)
        W_0_H = find_W_H(a_0,k,root)
        n_float = W_0_H - (beta/Z)
        """"n float is the normalized surface potential"""
        phi_grain = (k_b*T_e*n_float)/(e_charge)# unormalize it
        
        return phi_grain
    
    def produce_q_D(self):
        """ get dust charge using the surface potential, note the use of the grain radius"""
        return 4*np.pi*epsilon_0*grain_R*self.phi_grain

    def prod_W_vec(self):
        """when creating new dust grains gives them random x,y positions so the dust grains dont just stack on 0,0"""
        x_0 = -container_radius/6 + random.random()*container_radius/3#create random number centred about 0 with width container_radius/3
        y_0 = -container_radius/6 + random.random()*container_radius/3
        z_0 = self.container_height
        vx_0 = 0#random.random()*2e-10 - 1e-10
        vy_0 = 0#random.random()*2e-10 - 1e-10
        vz_0 = 0#random.random()*2e-8 - 1e-8
        #print("dust_init_pos=",x_0/lambda_D,y_0/lambda_D,self.container_height/lambda_D)
        return np.array([x_0,y_0,z_0,vx_0,vy_0,vz_0])#W_vec
    
    def calc_temperature(self):
        self.temperature = self.m_D*(self.v_squared_sum/len(self.Dust_grain_list))/(3*k_b)
    
    def inter_particle_ca(self):
        for i in self.comb_list:

            r_01 =  i[0].W_vec[0:3] - i[1].W_vec[0:3]
            r_01_mag = np.linalg.norm(r_01)
            r_01_pos =  i[0].wake_pos - i[1].W_vec[0:3]
            r_01_pos_mag = np.linalg.norm(r_01_pos)
            r_10_pos = i[1].wake_pos - i[0].W_vec[0:3]
            r_10_pos_mag = np.linalg.norm(r_10_pos)

            force_c_pos_01 = np.array([0,0,0])
            force_c_pos_10 = np.array([0,0,0])
            #force_c = ((i[0].charge*i[1].charge)/(4*np.pi*epsilon_0))*(r_01/(np.linalg.norm(r_01))**3)
            #force_c_pos_01 = ((i[0].charge*i[1].wake_charge)/(4*np.pi*epsilon_0))*(r_01_pos/(np.linalg.norm(r_01_pos))**3)
            #force_c_pos_10 = ((i[1].charge*i[0].wake_charge)/(4*np.pi*epsilon_0))*(r_10_pos/(np.linalg.norm(r_10_pos))**3)

            force_c = -((i[0].charge*i[1].charge)/(4*np.pi*epsilon_0))* np.exp((i[1].grain_R/lambda_D) - (r_01_mag/lambda_D)) * (1/(r_01_mag**3) + 1/(lambda_D*(r_01_mag**2)))* r_01
            
            if(i[1].W_vec[2] < z_se):
                
                    #M = ((v_B**2 + (i_charge*k_z_restore*(W_vec_f[2] - z_se)**2)/m_i -2*g_z*(W_vec_f[2] - z_se))**0.5)/v_B
                    #L_s = lambda_D*((M**2 - 1)**0.5)
                    #force_c_pos_01 = i[0].charge*i[1].wake_charge * (2/(1-(M**-2))) * ((np.cos(z/L_s)/z**2) + (np.sin(z/L_s)/z))
                    
                    force_c_pos_01 = -((i[0].charge*i[1].wake_charge)/(4*np.pi*epsilon_0))* np.exp((i[1].grain_R/lambda_D) - (r_01_pos_mag/lambda_D)) * (1/(r_01_pos_mag**3) + 1/(lambda_D*(r_01_pos_mag**2)))* r_01_pos
                    
            if(i[0].W_vec[2] < z_se):
                    force_c_pos_10 = -((i[1].charge*i[0].wake_charge)/(4*np.pi*epsilon_0))* np.exp((i[0].grain_R/lambda_D) - (r_10_pos_mag/lambda_D)) * (1/(r_10_pos_mag**3) + 1/(lambda_D*(r_10_pos_mag**2)))* r_10_pos
            
            i[0].a_c = i[0].a_c + force_c/i[0].mass + force_c_pos_01/i[0].mass
            i[1].a_c = i[1].a_c - force_c/i[1].mass + force_c_pos_10/i[1].mass
            
            
            
#v_i_z = (v_B**2 + (i_charge*k_z_restore*(W_vec_f[2] - z_se)**2)/m_i -2*g_z*(W_vec_f[2] - z_se))**0.5
    
    def next_frame(self):
        """The big daddy of the code, this functions loops over advancing the simulation by a time step dt"""
        
        #dust_list_copy = self.get_dust_list()
        self.inter_particle_ca()
            
        for i in np.arange(len(self.Dust_grain_list)):
            """advance via the rk4 and add the predicted postiosn to the momentary list """
            #print("hi ac =",self.Dust_grain_list[i].a_c)
            self.Dust_grain_list[i].time_list.append(time_list[-1]+ dt)
            self.Dust_grain_list[i].step()
            self.Dust_grain_list[i].a_c = np.array([0,0,0])
            self.Dust_grain_list[i].x_history.append(self.Dust_grain_list[i].W_vec[0]/lambda_D)#add the new postions to the history file
            self.Dust_grain_list[i].y_history.append(self.Dust_grain_list[i].W_vec[1]/lambda_D)
            self.Dust_grain_list[i].z_history.append(self.Dust_grain_list[i].W_vec[2]/lambda_D)
            self.v_squared_sum += (self.Dust_grain_list[i].calc_speed())**2

            if (self.Dust_grain_list[-1] == self.Dust_grain_list[i]) and (self.Dust_grain_list[i].W_vec[2] < z_se) and (len(self.Dust_grain_list) < dust_grain_max):
                """ if the last dust grain added has reached the lower sheathe electrode add another dust grain unless we have reached the dust grain number cap"""
                self.Dust_grain_list = np.append(self.Dust_grain_list,Dust_grain(self.m_D , self.grain_R , self.produce_q_D() , self.prod_W_vec(),time_list[-1]))
                self.comb_list = list(it.combinations(self.Dust_grain_list,2))
        self.calc_temperature()
        self.v_squared_sum = 0

                
#            if (self.Dust_grain_list[i].W_vec[2] <= self.Dust_grain_list[i].grain_R) or ((self.container_radius - (self.Dust_grain_list[i].W_vec[0]**2 + self.Dust_grain_list[i].W_vec[1]**2)**0.5) <= self.Dust_grain_list[i].grain_R):
#                """if statements checks whether grain is within a radius of the container walls, its never got there so not sure if the code acutally works"""
#                self.Dust_grain_list[i].to_be_deleted = True
#        
#        """add in to delete particels that hit the wall"""
#        copy_for_del_bollocks = self.get_dust_list()
#        for i in np.arange(len(self.Dust_grain_list)):
#            if self.Dust_grain_list[i].to_be_deleted == True:
#                copy_for_del_bollocks = np.delete(copy_for_del_bollocks,i)
#        if len(copy_for_del_bollocks) != len(self.Dust_grain_list):
#            self.Dust_grain_list = copy_for_del_bollocks 
#            self.comb_list = list(it.combinations(self.Dust_grain_list,2))

#%%
class Dust_grain(Dust_Container):
    """dust grain class handles all properties relating to the dust grain itself"""
    def __init__(self, m , grain_R , q , W_vec, time_init):
         self.mass = m
         self.grain_R = grain_R
         self.charge = q
         self.W_vec = W_vec
         #self.k_z_restore = -(1e-5)/self.charge#affects how strong the radial and vertical forces are due to sheathe
         #self.k_r_restore = -(1e-5)/self.charge
         self.a_c = np.array([0,0,0])
         
         self.time_list = [time_init]
         self.to_be_deleted = False
         self.x_history = [self.W_vec[0]/lambda_D]#records all x,yz positosn of the dust grain over history
         self.y_history = [self.W_vec[1]/lambda_D]
         self.z_history = [self.W_vec[2]/lambda_D]
         self.wake_pos = np.array([self.W_vec[0],self.W_vec[1],self.W_vec[2] - wake_potential_below])
         self.wake_charge = np.abs(self.charge)*wake_charge_multiplier
    
    def calc_speed(self):
        return np.linalg.norm(self.W_vec[3:6])

    """lots of rk4 stuff below, start at the bottom then work your way up, if you know what i mean... """
    
    def f_der(self, W_vec_f):
        x_diff_new = W_vec_f[3]
        y_diff_new = W_vec_f[4]
        z_diff_new = W_vec_f[5]
        
        radial = ((W_vec_f[0])**2 + (W_vec_f[1])**2)**0.5#radial distance from center
        
        """x,y components"""
        if radial < r_se_inv:
            vx_diff_new = - (alpha*W_vec_f[3])/self.mass + self.a_c[0]#drag and coloumb force
            vy_diff_new = - (alpha*W_vec_f[4])/self.mass + self.a_c[1]
        else:
            v_i_r = (v_B**2 + (i_charge*k_r_restore*(np.abs(radial - r_se_inv))**2)/m_i)**0.5
            a_i_r = (np.pi*self.grain_R**2*m_i*n_i0*v_i_r**2)/self.mass
            rad_force_mag = (self.charge/self.mass)*k_r_restore*np.abs(radial - r_se_inv)
            vx_diff_new = rad_force_mag*(W_vec_f[0]/radial) - (alpha*W_vec_f[3])/self.mass + self.a_c[0] + (W_vec_f[0]/radial)*a_i_r #drag, sheathe and coloumb force and ion drag force
            vy_diff_new = rad_force_mag*(W_vec_f[1]/radial) - (alpha*W_vec_f[4])/self.mass + self.a_c[1] + (W_vec_f[1]/radial)*a_i_r
            #print("x accel = ",rad_force_mag*(W_vec_f[0]/radial),- (alpha*W_vec_f[3])/self.mass, self.a_c[0] ,+ (W_vec_f[0]/radial)*a_i_r)
        if W_vec_f[2] > z_se:
            vz_diff_new = - g_z - (alpha*W_vec_f[5])/self.mass + self.a_c[2]#drag, gravity, coloumb force and ion drag force
        else:
            v_i_z = (v_B**2 + (i_charge*k_z_restore*(W_vec_f[2] - z_se)**2)/m_i -2*g_z*(W_vec_f[2] - z_se))**0.5
            #print("v_i_z =",v_i_z)
            #print(v_B**2, (i_charge*k_z_restore*(W_vec_f[2] - z_se)**2)/m_i,  -2*g_z*(W_vec_f[2] - z_se))
            vz_diff_new = (self.charge/self.mass)*k_z_restore*(W_vec_f[2] - z_se) - g_z - (alpha*W_vec_f[5])/self.mass  + self.a_c[2] + (np.pi*self.grain_R**2*m_i*n_i0*v_i_z**2)/self.mass#drag, sheathe, gravity, coloumb force and ion drag force
            #print("z accel = ", - g_z, - (alpha*W_vec_f[5])/self.mass, self.a_c[2], (np.pi*self.grain_R**2*m_i*n_i0*v_i_z**2)/self.mass)
        f = np.array([x_diff_new, y_diff_new, z_diff_new, vx_diff_new, vy_diff_new, vz_diff_new])
        return f
    
    """ rk4 is love, rk4 is life"""
    
    def k_1(self):
        res = self.f_der(self.W_vec)
        k_1 = dt*res
        return k_1
    
    def k_2(self,k_1):
        W_vec_k2 = self.W_vec + k_1/2
        res = self.f_der(W_vec_k2)
        k_2 = dt*res
        return k_2
    
    def k_3(self,k_2):
        W_vec_k3 = self.W_vec + k_2/2
        res = self.f_der(W_vec_k3)
        k_3 = dt*res
        return k_3
    
    def k_4(self,k_3):
        W_vec_k4 = self.W_vec + k_3
        res = self.f_der(W_vec_k4)
        k_4 = dt*res
    
        return k_4
    
    def step(self):
        
        k1 = self.k_1()
        k2 = self.k_2(k1)
        k3 = self.k_3(k2)
        k4 = self.k_4(k3)
        #print("hi")
        self.W_vec += (k1 + 2*k2 + 2*k3 + k4)/6
        self.wake_pos = np.array([self.W_vec[0],self.W_vec[1],self.W_vec[2] - wake_potential_below])  
#    def k_1_45(self,dt):
#        return dt*self.f_der(self.W_vec)
#    
#    def k_2_45(self,dt,k_1):
#        W_vec_k2 = self.W_vec + b_01*k_1
#        return dt*self.f_der(W_vec_k2)
#    
#    def k_3_45(self,dt,k_1,k_2):
#        W_vec_k3 = self.W_vec + b_31*k_1 + b_32*k_2
#        return dt*self.f_der(W_vec_k3)
#    
#    def k_4_45(self,dt,k_1,k_2,k_3):
#        W_vec_k4 = self.W_vec + b_41*k_1 + b_42*k_2 + b_43*k_3
#        return dt*self.f_der(W_vec_k4)
#    
#    def k_5_45(self,dt,k_1,k_2,k_3,k_4):
#        W_vec_k5 = self.W_vec + b_51*k_1 + b_52*k_2 + b_53*k_3 + b_54*k_4
#        return dt*self.f_der(W_vec_k5)
#    
#    def k_6_45(self,dt,k_1,k_2,k_3,k_4,k_5):
#        W_vec_k6 = self.W_vec + b_61*k_1 + b_62*k_2 + b_63*k_3 + b_64*k_4
#        return dt*self.f_der(W_vec_k6)
#    
#    def RK_fomulae(self,dt):
#        k1 = self.k_1_45(dt)
#        k2 = self.k_2_45(dt,k1)
#        k3 = self.k_3_45(dt,k1,k2)
#        k4 = self.k_4_45(dt,k1,k2,k3)
#        k5 = self.k_5_45(dt,k1,k2,k3,k4)
#        k6 = self.k_6_45(dt,k1,k2,k3,k4,k5)
#        W_vec_5 = self.W_vec + c_1*k1 + c_2*k2 + c_3*k3 + c_4*k4 + c_5*k5 + c_6*k6
#        W_vec_4 = self.W_vec + c_1_s*k1 + c_2_s*k2 + c_3_s*k3 + c_4_s*k4 + c_5_s*k5 + c_6_s*k6
#        return [W_vec_5,W_vec_4]
#
#    def step_rk45(self,dt):
#        delta_0 = np.absolute(eps*self.W_vec)
#        rk = self.RK_fomulae(dt)
#        W_vec_5 = rk[0]
#        W_vec_4 = rk[1]
#        delta_1 = np.absolute(W_vec_5 - W_vec_4)
#        R = np.amax(np.absolute(np.divide(delta_1,delta_0)))
#        if(R > 1):
#            dt_loop = S*dt*R**(-0.25)
#            while(R > 1):
#                rk = self.RK_fomulae(dt_loop)
#                W_vec_5 = rk[0]
#                W_vec_4 = rk[1]
#                delta_1 =  np.absolute(W_vec_5 - W_vec_4)
#                R = np.amax(np.absolute(np.divide(delta_1,delta_0)))
#                dt_loop = S*dt_loop*R**(-0.25)
#            dt_new = dt_loop
#        else:
#            dt_new = S*dt*(R**(-0.2))
#        return [W_vec_5,dt_new]
    
    
#%%
"""create the dusty plasma container"""
Dusty_plasma_crystal = Dust_Container(container_radius, container_height, m_D, grain_R)
if for_run:
    for i in np.arange(frames):
        """ do the loop advancing by dt each time"""
        time_list.append(time_list[-1]+ dt)
        Dusty_plasma_crystal.next_frame()
else:
    #Dust_average_speed = Dusty_plasma_crystal.calc_average_speed()
    Dusty_plasma_crystal.calc_temperature()
    while((len(Dusty_plasma_crystal.Dust_grain_list) <= dust_grain_max) and (len(time_list) < frames)):
        if ( (Dusty_plasma_crystal.temperature < temp_min) and (len(Dusty_plasma_crystal.Dust_grain_list) == dust_grain_max) ):
            break
        time_list.append(time_list[-1]+ dt)
        Dusty_plasma_crystal.next_frame()

#print("heeyyy")
#print(Dusty_plasma_crystal.temperature)

"""get data out"""
Final_conditions_dust_grains = Dusty_plasma_crystal.Dust_grain_list
vel_list = []
for i in np.arange(len(Final_conditions_dust_grains)):
    vel_list.append([Final_conditions_dust_grains[i].W_vec[3:6]])
    
#%%

theta = np.linspace(0, 2*np.pi, 100)
x_r_se = r_se_inv/lambda_D*np.cos(theta)
y_r_se = r_se_inv/lambda_D*np.sin(theta)
x_r_wall = container_radius/lambda_D*np.cos(theta)
y_r_wall = container_radius/lambda_D*np.sin(theta)


y_z_se = [z_se/lambda_D]*len(time_list)


#vor_pos_list = []
#for i in np.arange(len(Final_conditions_dust_grains)):
#    vor_pos_list.append([Final_conditions_dust_grains[i].x_history[-1],Final_conditions_dust_grains[i].y_history[-1]])
#vor = Voronoi(vor_pos_list)


#tri_pos_list = []
#for i in np.arange(len(Final_conditions_dust_grains)):
#    tri_pos_list.append([Final_conditions_dust_grains[i].x_history[-1],Final_conditions_dust_grains[i].y_history[-1],Final_conditions_dust_grains[i].z_history[-1]])
#tri = Delaunay(tri_pos_list)

#%%

plt.figure(1)
plt.title("test")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].x_history, Final_conditions_dust_grains[i].y_history)
    plt.plot(Final_conditions_dust_grains[i].x_history[0], Final_conditions_dust_grains[i].y_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].x_history[-1], Final_conditions_dust_grains[i].y_history[-1],"+" ,color='red')
plt.plot(x_r_se,y_r_se, "--", color = "black")
plt.plot(x_r_wall,y_r_wall, color = "black")
plt.xlabel("x/lambda_D")
plt.ylabel("y/lambda_D")
plt.grid()
plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

plt.figure(2)
plt.title("test - x")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].time_list,Final_conditions_dust_grains[i].x_history)
    plt.plot(Final_conditions_dust_grains[i].time_list[0],Final_conditions_dust_grains[i].x_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].time_list[-1],Final_conditions_dust_grains[i].x_history[-1],"+" ,color='red')
plt.xlabel("Time")
plt.ylabel("x/lambda_D")
plt.grid()
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

plt.figure(3)
plt.title("test - y")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].time_list,Final_conditions_dust_grains[i].y_history)
    plt.plot(Final_conditions_dust_grains[i].time_list[0],Final_conditions_dust_grains[i].y_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].time_list[-1],Final_conditions_dust_grains[i].y_history[-1],"+" ,color='red')
plt.xlabel("Time")
plt.ylabel("y/lambda_D")
plt.grid()
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

plt.figure(4)
plt.title("test - z")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].time_list,Final_conditions_dust_grains[i].z_history)
    plt.plot(Final_conditions_dust_grains[i].time_list[0],Final_conditions_dust_grains[i].z_history[0],"+" ,color='blue')
    plt.plot(Final_conditions_dust_grains[i].time_list[-1],Final_conditions_dust_grains[i].z_history[-1],"+" ,color='red')
plt.xlabel("Time")
plt.ylabel("z/lambda_D")
plt.plot(time_list,y_z_se, "--", color = "black")
#horiz_line_data = np.array([40 for i in xrange(len(xs))])
#xs = np.linspace(1,21,200)
#plt.plot(xs, horiz_line_data, 'r--') 
plt.ylim(0,container_height/lambda_D)
plt.grid()

plt.figure(5)
plt.title("Final Positions")
for i in np.arange(len(Final_conditions_dust_grains)):
    plt.plot(Final_conditions_dust_grains[i].x_history[-1], Final_conditions_dust_grains[i].y_history[-1],"+" ,color='red')
plt.xlabel("x/lambda_D")
plt.ylabel("y/lambda_D")
plt.plot(x_r_se,y_r_se, "--", color = "black")
plt.plot(x_r_wall,y_r_wall, color = "black")
plt.grid()
plt.xlim(-container_radius/lambda_D,container_radius/lambda_D)
plt.ylim(-container_radius/lambda_D,container_radius/lambda_D)

fig = plt.figure(6)
ax = fig.add_subplot(111, projection='3d')

for i in np.arange(len(Final_conditions_dust_grains)):
    ax.scatter(Final_conditions_dust_grains[i].x_history[-1], Final_conditions_dust_grains[i].y_history[-1], Final_conditions_dust_grains[i].z_history[-1])
ax.set_xlabel('X position')
ax.set_ylabel('Y position')
ax.set_zlabel('Z position') 

#"""Voronoi plot out, only 2d"""
#voronoi_plot_2d(vor)
#
#"""Delaunay tesselation"""
##plt.triplot(tri_pos_list[:,0], tri_pos_list[:,1], tri.simplices.copy())
##plt.plot(tri_pos_list[:,0], tri_pos_list[:,1], 'o')

print("steps = ",len(time_list))
print("Dust Grains =", len(Final_conditions_dust_grains))
print("time elapsed = ", time_list[-1],"s")
#print(Final_conditions_dust_grains[0].x_history[0]-Final_conditions_dust_grains[0].x_history[-1],Final_conditions_dust_grains[0].y_history[0]-Final_conditions_dust_grains[0].y_history[-1],Final_conditions_dust_grains[0].z_history[0]-Final_conditions_dust_grains[0].z_history[-1])

#%%
"""prints time taken in minutes"""
print ("time taken: %s minutes" % ((time.time()-start_time)/60))

#%%
plt.show()
