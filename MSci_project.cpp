#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream> 

#include <algorithm>
#include <vector>

using namespace std;

const float n_e0 = 1e15;//electron and ion densities in bulk plasma
const float n_i0 = 1e15;

const float g_z = 9.81;//gravity
const float e_charge = -1.6*1e-19;
const float i_charge = 1.6*1e-19;
const float grain_R = 7*1e-6;
const float m_i = 1.67*1e-27;
const float m_e = 9.11*1e-31;
const float m_D = ((4/3)*np.pi*pow(grain_R,3))*(1.49*1e3);//mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
const float epsilon_0 = 8.85*1e-12;
const float k_b = 1.38*1e-23;
const float mu = (m_i/m_e);//normalization used for convienience

const float beta = pow(10,-1);//T_i/T_e
const float T_e = 1e2;
const float T_i = T_e*beta;//floatined a bit wierdly but basically just in terms of beta and T_e
const float v_i = pow((2*k_b*T_i/m_i),0.5);
const float Z = 1;//atomic number??
const float a_0 = 1;//intial guess for halley's method

const float lambda_de = pow(((epsilon_0*k_b*T_e)/(n_e0*(pow(e_charge,2)))),0.5);
const float lambda_di = pow(((epsilon_0*k_b*T_i)/(n_i0*(pow(e_charge,2)))),0.5);
const float lambda_D = pow((1/(1/(pow(lambda_de,2)) + 1/(pow(lambda_di,2)))),0.5);//get lambda_D

const float container_height = 11*lambda_D;//drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
const float container_radius = 11*lambda_D;//set radius of contianer ie wall radius
const float z_se = 10*lambda_D;//distance from bottom of container to the sheath edge
const float r_se = 10*lambda_D;//distance from wall to the sheathe edge
const float r_se_inv = container_radius - r_se;

const float phi_wall_z = -100;//volts
const float phi_wall_r = -1;//volts
const float k_z_restore = -2*phi_wall_z/pow(z_se,2);//WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
const float k_r_restore = -2*phi_wall_r/pow(r_se,2);

const float v_B = pow((k_b*T_e/m_i),0.5);

const float alpha = 1e-9;//drag coefficient
float time_list[1]; //HOW LONG??
const float root = 1e-14;//floatined preciscion of root finding method used to get dust charge

const float dt = 1e-4;//time step in rk4, needs to be small enough to be precise but large enough we can actually move the stuff forward in time
const float dust_grain_max = 50;//dust grain max number
const float frames = 1e4;//number of frames, time taken is not linear as teh longer u run it the more particles it adds hence increases quadratically
const float speed_min = 1e-15;
const bool for_run = true;

class Dust_Container{
    public:
    Dust_Container(float container_radius, float container_height, float m_D, float grain_R){
        //DO I EVEN NEED THE NEXT COUPLE LINES?
            container_radius = container_radius;
            container_height = container_height;
            m_D = m_D;
            grain_R = grain_R;
            float phi_grain = self.OML_find_surface_potential();
            float Dust_grain_list [dust_grain_max];
            Dust_grain_list[0] = Dust_grain(m_D , grain_R , produce_q_D() , prod_W_vec(), time_list[-1])
            //float comb_list = list(it.combinations(self.Dust_grain_list,2))
    }

    vector<float> time_list;

    float OML_find_surface_potential(){
        //Find the dust grain sufrace potential using OML model
        float find_W_H(a_0,z,root_prec){
            //root finding method using "halley's method" like newtons method but third order
            float f_x(x_n,z){
                return float x_n*np.exp(x_n) - z
            }
            
            float f_x_first_derv(x_n){
                return float np.exp(x_n)*(1 + x_n)
            }
            float f_x_second_derv(x_n){
                return float np.exp(x_n)*(2 + x_n)
            }
            float calc_x_plus(x_n,z){
                f_x_0 = f_x(x_n,z)
                f_x_1 = f_x_first_derv(x_n)
                f_x_2 = f_x_second_derv(x_n)
                x_plus = x_n - ((2*f_x_0*f_x_1)/(2*(pow(f_x_1,2))-f_x_0*f_x_2))
                return x_plus
            }
            float a_n = a_0
            float a_plus = calc_x_plus(a_0,z)
        
            while(np.abs(a_plus - a_n) > root_prec){
                a_n = a_plus
                a_plus = calc_x_plus(a_n,z)
            }
            float W_0 = (a_n + a_plus)/2
            return W_0
        }

        //basically we have an equation that has the normalized surface potetnail in it and we rootfind for it, using the lambert W function form though tbh i agree with coppins this isnt actually necesary
        float k = ((pow((mu*beta),0.5))/Z)*np.exp(beta/Z)
        float W_0_H = find_W_H(a_0,k,root)
        float n_float = W_0_H - (beta/Z)
        ///n float is the normalized surface potential
        float phi_grain = (k_b*T_e*n_float)/(e_charge)// unormalize it
        
        return phi_grain
    }

    float produce_q_D(){
        // get dust charge using the surface potential, note the use of the grain radius//
        return float 4*np.pi*epsilon_0*grain_R*.phi_grain
    }
    float prod_W_vec(){
        //when creating new dust grains gives them random x,y positions so the dust grains dont just stack on 0,0//
        float x_0 = -container_radius/6 + random()*container_radius/3//create random number centred about 0 with width container_radius/3
        float y_0 = -container_radius/6 + random()*container_radius/3
        float z_0 = container_height
        float vx_0 = random()*2e-10 - 1e-10
        float vy_0 = random()*2e-10 - 1e-10
        float vz_0 = random()*2e-10 - 1e-10

        return np.array([x_0,y_0,z_0,vx_0,vy_0,vz_0])#W_vec
    }

    float calc_average_speed(){
        float speed = 0
        for( int i = 0; i < sizeof(in Dust_grain_list); i+=1){
            float a = i.W_vec[36]
            float b = np.linalg.norm(a)
            float speed += b
        } 
        return speed/sizeof(Dust_grain_list)
    }  

    /*
    float inter_particle_ca(){
        for i in .comb_list
            r_21 = i[0].W_vec[03] - i[1].W_vec[03]
            force_c = ((i[0].charge*i[1].charge)/(4*np.pi*epsilon_0))*(r_21/(np.linalg.norm(r_21))**3)
            i[0].a_c = i[0].a_c + force_c/i[0].mass
            i[1].a_c = i[1].a_c - force_c/i[1].mass
    }
    */

    void next_frame(){
    //The big daddy of the code, this functions loops over advancing the simulation by a time step dt
        //inter_particle_ca()
            
        for (i = 0; i < sizeof(in Dust_grain_list); i+=1){
            //advance via the rk4 and add the predicted postiosn to the momentary list
            Dust_grain_list[i].step();
            Dust_grain_list[i].time_list.push_back(time_list[-1]+ dt);
            Dust_grain_list[i].a_c[3] = {0,0,0}//NOT SURE ABOUT THIS
            Dust_grain_list[i].x_history.push_back(Dust_grain_list[i].W_vec[0]/lambda_D )
            Dust_grain_list[i].y_history.push_back( Dust_grain_list[i].W_vec[1]/lambda_D)
            Dust_grain_list[i].z_history.push_back(Dust_grain_list[i].W_vec[2]/lambda_D)

            if ( (Dust_grain_list[-1] == Dust_grain_list[i]) && (Dust_grain_list[i].W_vec[2] < z_se) && (sizeof(Dust_grain_list) < dust_grain_max) ){
                """ if the last dust grain added has reached the lower sheathe electrode add another dust grain unless we have reached the dust grain number cap"""
                Dust_grain_list = np.append(Dust_grain_list,Dust_grain(m_D , grain_R , produce_q_D() , prod_W_vec(),time_list[-1]))//FIXXXXXX
                comb_list = list(it.combinations(Dust_grain_list,2))//FIXXXXX
            }
        }
    }

}

class Dust_grain(){
    //dust grain class handles all properties relating to the dust grain it
    public:
    Dust_grain(float m , float grain_R , float q , float W_vec, float time_init){
        mass = m
        grain_R = grain_R
        charge = q
        W_vec = W_vec
        float  a_c[3] = {0,0,0}
        vector<float> time_list = {time_init};
        bool to_be_deleted = False
        vector<float> x_history = {W_vec[0]/lambda_D}
        vector<float> y_history = {W_vec[1]/lambda_D}
        vector<float> z_history = {W_vec[2]/lambda_D}
    }
    
    float f_der(float W_vec_f){
        float x_diff_new = W_vec_f[3]
        float y_diff_new = W_vec_f[4]
        float z_diff_new = W_vec_f[5]
        float vx_diff_new;
        float vy_diff_new;
        float vz_diff_new;
        float radial = pow((pow(W_vec_f[0],2) + pwo(W_vec_f[1],2)),0.5);
        
        //x,y components
        if (radial < r_se_inv){
            vx_diff_new = - (alpha*W_vec_f[3])/mass + a_c[0]//drag and coloumb force
            vy_diff_new = - (alpha*W_vec_f[4])/mass + a_c[1]//FIX THIISSSSSSSS
        }
        else{
            float v_i_r = (pow(v_B,2) + pow((i_charge*k_r_restore*pow((np.abs(radial - r_se_inv)),2)/m_i),0.5)
            float a_i_r = (np.pi*pow(grain_R,2)*m_i*n_i0*pow(v_i_r,2))/mass
            float rad_force_mag = (charge/mass)*k_r_restore*np.abs(radial - r_se_inv)
            vx_diff_new = rad_force_mag*(W_vec_f[0]/radial) - (alpha*W_vec_f[3])/mass + a_c[0] + (W_vec_f[0]/radial)*a_i_r//drag, sheathe and coloumb force and ion drag force
            vy_diff_new = rad_force_mag*(W_vec_f[1]/radial) - (alpha*W_vec_f[4])/mass + a_c[1] + (W_vec_f[1]/radial)*a_i_r
        }

        if (W_vec_f[2] > z_se){
            vz_diff_new = - g_z - (alpha*W_vec_f[5])/mass + a_c[2]//drag, gravity, coloumb force and ion drag force
        }
        else{
            float v_i_z = pow((pow(v_B,2) + (i_charge*k_z_restore*pow((W_vec_f[2] - z_se),2))/m_i -2*g_z*(W_vec_f[2] - z_se)),0.5);
            vz_diff_new = (charge/mass)*k_z_restore*(W_vec_f[2] - z_se) - g_z - (alpha*W_vec_f[5])/mass  + a_c[2] + (np.pi*pow(grain_R,2)*m_i*n_i0*pow(v_i_z,2))/mass;//drag, sheathe, gravity, coloumb force and ion drag force;
        }
        float f[6] = {x_diff_new, y_diff_new, z_diff_new, vx_diff_new, vy_diff_new, vz_diff_new};

        return f;
    }
    
    //rk4 is love, rk4 is life
    
    float k_1(){
        float res = f_der(W_vec);
        float k_1 = dt*res;
        return k_1;
    }
    
    float k_2(float k_1){
        float W_vec_k2 = W_vec + k_1/2;
        float res = f_der(W_vec_k2);
        float k_2 = dt*res;
        return k_2;
    }
    
    float k_3(float k_2){
        float W_vec_k3 = W_vec + k_2/2;
        float res = f_der(W_vec_k3);
        float k_3 = dt*res;
        return k_3;
    }
    
    float k_4(float k_3){
        float W_vec_k4 = W_vec + k_3;
        float res = f_der(W_vec_k4);
        float k_4 = dt*res;
        return k_4;
    }
    
    void step(){
        float k1 = k_1();
        float k2 = k_2(k1);
        float k3 = k_3(k2);
        float k4 = k_4(k3);
        W_vec = W_vec + (k1 + 2*k2 + 2*k3 + k4)/6;
    }
}

int main(){
    //create the dusty plasma container
    Dusty_plasma_crystal = Dust_Container(container_radius, container_height, m_D, grain_R)
    if (for_run){
        for(int i = 0; i < frames; i++){
            """ do the loop advancing by dt each time"""
            time_list.push_back(time_list[-1]+ dt)
            Dusty_plasma_crystal.next_frame()
        }
    }
    else{
        Dust_average_speed = Dusty_plasma_crystal.calc_average_speed()
        while( sizeof(Dusty_plasma_crystal.Dust_grain_list) <= dust_grain_max){
            Dust_average_speed = Dusty_plasma_crystal.calc_average_speed()
            if (Dust_average_speed < speed_min){
                break
            }
            time_list.push_back(time_list[-1]+ dt)
            Dusty_plasma_crystal.next_frame()    
        }
    }

    // write out data
    ofstream ofs("saved.txt");
    for(int i = 0; i < sizeof(Dusty_plasma_crystal.Dust_grain_list); i++){
	    ofs << Dusty_plasma_crystal.Dust_grain_list[i]; //store the object to file
    }
	ofs.close();

    return 0;
}

    

