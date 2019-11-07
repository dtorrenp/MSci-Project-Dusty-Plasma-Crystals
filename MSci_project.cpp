#define _USE_MATH_DEFINES

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream> 

#include <algorithm>
#include <vector>
#include <string>

using namespace std;

const double n_e0 = 1e15;//electron and ion densities in bulk plasma
const double n_i0 = 1e15;

const double g_z = 9.81;//gravity
const double e_charge = -1.6*1e-19;
const double i_charge = 1.6*1e-19;
const double grain_R = 7*1e-6;
const double m_i = 1.67*1e-27;
const double m_e = 9.11*1e-31;
const double m_D = ((4/3)*M_PI*pow(grain_R,3))*(1.49*1e3);//mass of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
const double epsilon_0 = 8.85*1e-12;
const double k_b = 1.38*1e-23;
const double mu = (m_i/m_e);//normalization used for convienience

const double beta = pow(10,-1);//T_i/T_e
const double T_e = 1e2;
const double T_i = T_e*beta;//doubleined a bit wierdly but basically just in terms of beta and T_e
const double v_i = pow((2*k_b*T_i/m_i),0.5);
const double Z = 1;//atomic number??
const double a_0 = 1;//intial guess for halley's method

const double lambda_de = pow(((epsilon_0*k_b*T_e)/(n_e0*(pow(e_charge,2)))),0.5);
const double lambda_di = pow(((epsilon_0*k_b*T_i)/(n_i0*(pow(e_charge,2)))),0.5);
const double lambda_D = pow((1/(1/(pow(lambda_de,2)) + 1/(pow(lambda_di,2)))),0.5);//get lambda_D

const double container_height = 11*lambda_D;//drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
const double container_radius = 11*lambda_D;//set radius of contianer ie wall radius
const double z_se = 10*lambda_D;//distance from bottom of container to the sheath edge
const double r_se = 10*lambda_D;//distance from wall to the sheathe edge
const double r_se_inv = container_radius - r_se;

const double phi_wall_z = -100;//volts
const double phi_wall_r = -1;//volts
const double k_z_restore = -2*phi_wall_z/pow(z_se,2);//WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
const double k_r_restore = -2*phi_wall_r/pow(r_se,2);

const double v_B = pow((k_b*T_e/m_i),0.5);

const double alpha = 1e-9;//drag coefficient

const double root = 1e-14;//preciscion of root finding method used to get dust charge

const double dt = 1e-4;//time step in rk4, needs to be small enough to be precise but large enough we can actually move the stuff forward in time
const double dust_grain_max = 5;//dust grain max number
const double frames = 1e5;//number of frames, time taken is not linear as teh longer u run it the more particles it adds hence increases quadratically
const double speed_min = 1e-15;
const bool for_run = true;

//make functions for element-wise multiplication and addition of vectors
vector<double> element_mul(const vector<double>& a,double cst){
    vector<double> c;
    for(int i=0; i<a.size(); i++){
        c.push_back(a[i]*cst);
    }
    return c;
}

vector<double> element_add(const vector<double>& a,const vector<double>& b){
    vector<double> c;
    for(int i=0; i<a.size(); i++){
        c.push_back(a[i] + b[i]);
    }
    return c;
}

class Dust_grain{
    //dust grain class handles all properties relating to the dust grain it

    //PUT IT FIRST COS COMPILER GETS ANGSTY AS CONTAINER CALLS IT
    private:
    double mass;
    double grain_R;
    double charge;
    
    public:

    vector<double>W_vec;
    vector<double>a_c;
    vector<double>time_list_dust;
    bool to_be_deleted;
    vector<double> x_history;
    vector<double> y_history;
    vector<double> z_history;

    Dust_grain(double m , double grain_r , double q, double time_init):mass(m),grain_R(grain_r),charge(q),a_c{0,0,0},time_list_dust{time_init},to_be_deleted(false)
    {
        prod_W_vec();        
    }
    
    void prod_W_vec(){
        //when creating new dust grains gives them random x,y positions so the dust grains dont just stack on 0,0//
        W_vec.push_back(-container_radius/6.0 + (rand()/(RAND_MAX + 1.0))*container_radius/3.0);//create random number centred about 0 with width container_radius/3
        W_vec.push_back(-container_radius/6.0 + (rand()/(RAND_MAX + 1.0))*container_radius/3.0);
        W_vec.push_back(container_height);
        W_vec.push_back((rand()/(RAND_MAX + 1.0))*2.0e-12 - 1.0e-12);
        W_vec.push_back((rand()/(RAND_MAX + 1.0))*2.0e-12 - 1.0e-12);
        W_vec.push_back((rand()/(RAND_MAX + 1.0))*2.0e-12 - 1.0e-12);
        x_history.push_back(W_vec[0]/lambda_D);
        y_history.push_back(W_vec[1]/lambda_D);
        z_history.push_back(W_vec[2]/lambda_D);

        //cout << "W_vec = " << W_vec[0]/lambda_D << "," << W_vec[1]/lambda_D << "," << W_vec[2]/lambda_D << "," << W_vec[3] << ","  << W_vec[4] << ","  << W_vec[5] << endl;
    }

    vector<double> f_der(vector<double> W_vec_f){
        double x_diff_new = W_vec_f[3];
        double y_diff_new = W_vec_f[4];
        double z_diff_new = W_vec_f[5];
        double vx_diff_new;
        double vy_diff_new;
        double vz_diff_new;
        double radial = pow((pow(W_vec_f[0],2) + pow(W_vec_f[1],2)),0.5);
        
        //x,y components
        if (radial < r_se_inv){
            vx_diff_new = - (alpha*W_vec_f[3])/mass + a_c[0];//drag and coloumb force
            vy_diff_new = - (alpha*W_vec_f[4])/mass + a_c[1];
        }
        else{
            double v_i_r = pow(v_B,2) + pow((i_charge*k_r_restore*pow((abs(radial - r_se_inv)),2)/m_i),0.5);

            //cout << "abs:" << radial - r_se_inv << "," << abs(radial - r_se_inv) << endl;

            double a_i_r = (M_PI*pow(grain_R,2)*m_i*n_i0*pow(v_i_r,2))/mass;
            double rad_force_mag = (charge/mass)*k_r_restore*abs(radial - r_se_inv);

            //cout << "heeyyy" << endl;
            //cout << "rad_force_mag sheathe:" << rad_force_mag << endl;
            //cout << "drag:" << (alpha*W_vec_f[3])/mass << endl;
            //cout << "inter particle:" <<  a_c[0] << endl;
            //cout << "ion drag:" << (W_vec_f[0]/radial)*a_i_r << endl;
            

            vx_diff_new = rad_force_mag*(W_vec_f[0]/radial) - (alpha*W_vec_f[3])/mass + a_c[0] + (W_vec_f[0]/radial)*a_i_r;//drag, sheathe and coloumb force and ion drag force
            cout << "x accel" << rad_force_mag*(W_vec_f[0]/radial) << "," << - (alpha*W_vec_f[3])/mass << "," << (W_vec_f[0]/radial)*a_i_r << endl;

            vy_diff_new = rad_force_mag*(W_vec_f[1]/radial) - (alpha*W_vec_f[4])/mass + a_c[1] + (W_vec_f[1]/radial)*a_i_r;
        }

        if (W_vec_f[2] > z_se){
            vz_diff_new = - g_z - (alpha*W_vec_f[5])/mass + a_c[2];//drag, gravity, coloumb force and ion drag force
        }
        else{
            double v_i_z = pow((pow(v_B,2) + (i_charge*k_z_restore*pow((W_vec_f[2] - z_se),2))/m_i -2.0*g_z*(W_vec_f[2] - z_se)),0.5);
            vz_diff_new = (charge/mass)*k_z_restore*(W_vec_f[2] - z_se) - g_z - (alpha*W_vec_f[5])/mass  + a_c[2] + (M_PI*pow(grain_R,2)*m_i*n_i0*pow(v_i_z,2))/mass;//drag, sheathe, gravity, coloumb force and ion drag force;
        }

        vector<double> f = {x_diff_new, y_diff_new, z_diff_new, vx_diff_new, vy_diff_new, vz_diff_new};
        return f;
    }
    
    //rk4 is love, rk4 is life
    
    vector<double> k_1(){
        vector<double> res = f_der(W_vec);
        vector<double> k_1_v = element_mul(res,dt);
        return k_1_v;
    }
    
    vector<double> k_2(vector<double> k_1){
        vector<double> W_vec_k2 = element_add(W_vec,element_mul(k_1,1/2) );
        vector<double> res = f_der(W_vec_k2);
        vector<double> k_2_v = element_mul(res,dt);
        return k_2_v;
    }
    
    vector<double> k_3(vector<double> k_2){
        vector<double> W_vec_k3 = element_add(W_vec,element_mul(k_2,1/2));
        vector<double> res = f_der(W_vec_k3);
        vector<double> k_3_v = element_mul(res,dt);
        return k_3_v;
    }
    
    vector<double> k_4(vector<double> k_3){
        vector<double> W_vec_k4 = element_add(W_vec,k_3);
        vector<double> res = f_der(W_vec_k4);
        vector<double> k_4_v = element_mul(res,dt);
        return k_4_v;
    }
    
    void step(){
        vector<double> k1 = k_1();
        vector<double> k2 = k_2(k1);
        vector<double> k3 = k_3(k2);
        vector<double> k4 = k_4(k3);

        cout << "hi" << endl;
        //cout << "k1 = " <<k1[0]/lambda_D << "," << k1[1]/lambda_D << "," << k1[2]/lambda_D << "," << k1[3] << ","  << k1[4] << ","  << k1[5] << endl;
        //cout << "k2 = " <<k2[0]/lambda_D << "," << k1[1]/lambda_D << "," << k1[2]/lambda_D << "," << k1[3] << ","  << k1[4] << ","  << k1[5] << endl;
        //cout << "k3 = " <<k3[0]/lambda_D << "," << k1[1]/lambda_D << "," << k1[2]/lambda_D << "," << k1[3] << ","  << k1[4] << ","  << k1[5] << endl;
        //cout << "k4 = " <<k4[0]/lambda_D << "," << k1[1]/lambda_D << "," << k1[2]/lambda_D << "," << k1[3] << ","  << k1[4] << ","  << k1[5] << endl;

        W_vec = element_add(W_vec,element_mul(element_add(element_add(    element_add(k1,element_mul(k2,2)) ,element_mul(k3,2)      ),k4),1.0/6.0));

        //cout << "W_vec = " << W_vec[0]/lambda_D << "," << W_vec[1]/lambda_D << "," << W_vec[2]/lambda_D << "," << W_vec[3] << ","  << W_vec[4] << ","  << W_vec[5] << endl;
    }
};

class Dust_Container{

    private:

    double container_radius;
    double container_height;
    double m_D;
    double grain_R;
    double phi_grain;

    public:
	vector<double> time_list;
    vector<Dust_grain> Dust_grain_list;

	Dust_Container(double container_Radius, double container_Height, double M_D, double Grain_R) :time_list{0},container_radius(container_Radius), container_height(container_Height), m_D(M_D), grain_R(Grain_R)
    {
            phi_grain = OML_find_surface_potential();
            Dust_grain_list = {Dust_grain(m_D , grain_R , produce_q_D(), 0)};
            //double comb_list = list(it.combinations(self.Dust_grain_list,2))
    } 

    //Find the dust grain sufrace potential using OML model
    double OML_find_surface_potential(){
        //basically we have an equation that has the normalized surface potetnail in it and we rootfind for it, using the lambert W function form though tbh i agree with coppins this isnt actually necesary
        double k = ((pow((mu*beta),0.5))/Z)*exp(beta/Z);
        double W_0_H = find_W_H(a_0,k,root);
        double n_double = W_0_H - (beta/Z);
        ///n double is the normalized surface potential
        return (k_b*T_e*n_double)/(e_charge);// unormalize it
    }
    double f_x(double x_n, double z){
        return x_n*exp(x_n) - z;
    }
    double f_x_first_derv(double x_n){
        return exp(x_n)*(1 + x_n);
    }
    double f_x_second_derv(double x_n){
        return exp(x_n)*(2 + x_n);
    }
    double calc_x_plus(double x_n,double z){
        double f_x_0 = f_x(x_n,z);
        double f_x_1 = f_x_first_derv(x_n);
        double f_x_2 = f_x_second_derv(x_n);
        double x_plus = x_n - ((2*f_x_0*f_x_1)/(2*(pow(f_x_1,2))-f_x_0*f_x_2));
        return x_plus;
    }
    double find_W_H(double a_0,double z,double root_prec){
        //root finding method using "halley's method" like newtons method but third order
        double a_n = a_0;
        double a_plus = calc_x_plus(a_0,z);//do i need to give type if product of a thing?
    
        while(abs(a_plus - a_n) > root_prec){
            a_n = a_plus;
            a_plus = calc_x_plus(a_n,z);
        }
        double W_0 = (a_n + a_plus)/2;
        return W_0;
    }
    double produce_q_D(){
        // get dust charge using the surface potential, note the use of the grain radius//
        return 4*M_PI*epsilon_0*grain_R*phi_grain;
    }
    /*
    double calc_average_speed(){
        double speed = 0
        for( int i = 0; i < sizeof(in Dust_grain_list); i+=1){
            double a = i.W_vec[36]
            double b = np.linalg.norm(a)
            double speed += b
        } 
        return speed/sizeof(Dust_grain_list)
    } 
    */ 

    /*
    double inter_particle_ca(){
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

		time_list.push_back(time_list.back() + dt);

        for (int i = 0; i < Dust_grain_list.size(); i+=1){
            //advance via the rk4 and add the predicted postiosn to the momentary list
            Dust_grain_list[i].step();
            Dust_grain_list[i].time_list_dust.push_back(time_list.back()+ dt);
            Dust_grain_list[i].a_c = {0,0,0};//NOT SURE ABOUT THIS
            Dust_grain_list[i].x_history.push_back(Dust_grain_list[i].W_vec[0]/lambda_D);
            Dust_grain_list[i].y_history.push_back( Dust_grain_list[i].W_vec[1]/lambda_D);
            Dust_grain_list[i].z_history.push_back(Dust_grain_list[i].W_vec[2]/lambda_D);

            //cout << Dust_grain_list[i].W_vec[0]/lambda_D << "," << Dust_grain_list[i].W_vec[1]/lambda_D << "," << Dust_grain_list[i].W_vec[2]/lambda_D << endl;

            if ( (i == (Dust_grain_list.size() - 1)) && (Dust_grain_list[i].W_vec[2] < z_se) && (Dust_grain_list.size() < dust_grain_max) ){
                //""" if the last dust grain added has reached the lower sheathe electrode add another dust grain unless we have reached the dust grain number cap"""
                Dust_grain_list.push_back(Dust_grain(m_D , grain_R , produce_q_D(),time_list.back()));
                //comb_list = list(it.combinations(Dust_grain_list,2))//FIXXXXX
            }
        }
    }
};
int main(){
    //create the dusty plasma container
    Dust_Container Dusty_plasma_crystal = Dust_Container(container_radius, container_height, m_D, grain_R);
    if (for_run){
        for(int i = 0; i < frames; i++){
            //""" do the loop advancing by dt each time"""
            Dusty_plasma_crystal.next_frame();
        }
    }
    /*
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
    */

    //

    ofstream outputFile;
    outputFile.open("test_csv.csv");
    outputFile << "X" << "," << "Y" << "," << "Z" << "," << "VX" << "," << "VY" << "," << "VZ" << endl;
    for(int i = 0; i < Dusty_plasma_crystal.Dust_grain_list.size(); i++){
	    outputFile << Dusty_plasma_crystal.Dust_grain_list[i].W_vec[0]/lambda_D << "," << Dusty_plasma_crystal.Dust_grain_list[i].W_vec[1]/lambda_D << ","  << Dusty_plasma_crystal.Dust_grain_list[i].W_vec[2]/lambda_D << ","  << Dusty_plasma_crystal.Dust_grain_list[i].W_vec[3] << ","  << Dusty_plasma_crystal.Dust_grain_list[i].W_vec[4] << ","  << Dusty_plasma_crystal.Dust_grain_list[i].W_vec[5] << endl;
    }
    outputFile.close();
	
    return 0;
}

    

