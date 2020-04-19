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
#include <random>
#include <chrono>
#include <map>

//CRITICAL VALUES
const int dust_grain_max_input = 20; //dust grain max number
const double dt_a = 1.0e-4;
const double time_limit = 1;
const double frame_req = 5;
const double V_DC = -49.9655;  //V_RF = 50 , rounded to 4 dp, taken from the python

//CONSTANTS TO FUCK ABOUT WITH
const double n_e0 = 1.0e15; //electron number density (in bulk plasma)
const double n_i0 = 1.0e15; //ion number density (in bulk plasma)
const double n_n0 = 1.0e24; //electron number density (in bulk plasma)
const double Z = 1; //ion atomic number
const double Z_n = 18; //neutrals atomic number (Ar)
const double grain_R = 7*1e-6; //dust grain radius
const double dust_grain_density = 1.49*1e3; //dust density
const double phi_wall_r = -1000.0; //radial wall potential [volts]
const double init_speed = 1e-5;

//CONSTANTS DEPENDANT ON ACTUAL PHYSICS
const double g_z = 9.81;//gravity
const double e_charge = 1.6*1e-19;//magnitude of e charge
const double i_charge = 1.6*1e-19;//magnitude of i charge DOES THIS NEED TO CHANGE WHEN USED IN THE ION DRAG?????
const double m_i = 1.67*1e-27;
const double m_e = 9.11*1e-31;
const double m_n = Z_n*m_i;
const double m_D = ((4.0/3.0)*M_PI*pow(grain_R,3))*dust_grain_density;//m_D of the dust grain given by volume*density where the density is space dust NEED TO GET BETTER VALUE
const double epsilon_0 = 8.85*1e-12;
const double k_b = 1.38*1e-23;
const double mu = (m_i/m_e);//normalization used for convienience
const double T_e = 2.0*(1.6*1e-19)/k_b;
const double T_i = 0.03*(1.6*1e-19)/k_b;
const double T_n = 0.03*(1.6*1e-19)/k_b;
const double beta = T_i/T_e;
const double lambda_de = pow(((epsilon_0*k_b*T_e)/(n_e0*(pow(e_charge,2)))),0.5);
const double lambda_di = pow(((epsilon_0*k_b*T_i)/(n_i0*(pow(e_charge,2)))),0.5);
const double lambda_D = pow((1/(1/(pow(lambda_de,2)) + 1/(pow(lambda_di,2)))),0.5);//get lambda_D
const double drop_height = 9.8*lambda_D;//drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
const double container_radius = 25.0*lambda_D;//set radius of contianer ie wall radius
const double coulomb_limit = 5;

//related to finding the charge and other set up stuff of the dust grains
const double Y_DC = -50.0;//NITTER
const double dz_norm = 1/1000.0;
const double dz = lambda_D*dz_norm;
const double upper_lim_z = 15.0*lambda_D;
const int for_loop_max_int = upper_lim_z/dz;
const double root = dz/10;//preciscion of root finding method used to get dust charge
const double Y_epsilon = -1e-2;//how close we define the sheath edge
const double a_0_Sheath = -2.5;//intial guess for halley's method
const double r_se = 25.0*lambda_D; //distance of radial 'sheath' from wall of container
const double v_B = pow((3*k_b*T_e/m_i),0.5);//bohm velocity of sheath
const double v_Tn = pow((3*k_b*T_n/m_n),0.5);//thermal termperature of the neutrals
const double alpha_n = (8/3)*pow(2*M_PI,0.5)*pow(grain_R,2)*m_n*n_n0*v_Tn;//taken from shuckla
const double therm_coeff = sqrt(2*k_b*T_n*alpha_n/dt_a);//for brownian motion
const double wake_safety_factor = grain_R*0.5;


const double k_r_restore = -2.0*phi_wall_r/pow(r_se,2);


//make functions for element-wise multiplication and addition of vectors
std::vector<double> element_mul(const std::vector<double>& a,double cst){
    std::vector<double> c;
    for(int i=0; i<a.size(); i++){
        c.push_back(a[i]*cst);
    }
    return c;
}

std::vector<double> element_add_cst(const std::vector<double>& a,double cst){
    std::vector<double> c;
    for(int i=0; i<a.size(); i++){
        c.push_back(a[i] + cst);
    }
    return c;
}


std::vector<std::vector<double>> special_element_mul(const std::vector<std::vector<double>>& a,double cst){
    std::vector<std::vector<double>> c;
    for(int i=0; i<a.size(); i++){
        for(int v=0; v<a[0].size(); v++){
            c[i].push_back(a[i][v]*cst);
        }
    }
    return c;
}

std::vector<double> element_add(const std::vector<double>& a,const std::vector<double>& b){
    std::vector<double> c;
    for(int i=0; i<a.size(); i++){
        c.push_back(a[i] + b[i]);
    }
    return c;
}

double v_abs(const std::vector<double>& a){
    double c = 0.0;
    for(int i=0; i<a.size(); i++){
        c += pow(a[i],2);
    }
    return pow(c,0.5);
}

class Dust_grain{
    //dust grain class handles all properties relating to the dust grain it
    private:
    std::default_random_engine generator;
    std::normal_distribution<double> dist;
    public:
    std::vector<double>W_vec;
    std::vector<double>a_c;
    std::vector<double> x_history;
    std::vector<double> y_history;
    std::vector<double> z_history;
    std::vector<double> vx_history;
    std::vector<double> vy_history;
    std::vector<double> vz_history;
    double charge;
    double Sheath_E;
    double z_se;
    double v_i;
    double ion_drag;

    Dust_grain():a_c{0,0,0},
    generator(std::default_random_engine(clock())),dist(std::normal_distribution<double>(0.0,sqrt((k_b*T_i)/m_D)))
    {
        prod_W_vec();      
    }

    double calc_speed(){
        std::vector<double> vel (W_vec.end() - 3,W_vec.end());        
        return v_abs(vel);
    }
    
    void prod_W_vec(){
        //when creating new dust grains gives them random x,y positions so the dust grains dont just stack on 0,0
        //std::cout << container_dust_dist_creation << "," <<(-container_dust_dist_creation/2.0 + (rand()/(RAND_MAX + 1.0))*container_dust_dist_creation)/lambda_D << std::endl;
        W_vec.push_back(-container_radius /4.0 + (rand()/(RAND_MAX + 1.0))*container_radius/2);//create random number centred about 0 with width container_radius/3
        W_vec.push_back(-container_radius /4.0 + (rand()/(RAND_MAX + 1.0))*container_radius/2);
        W_vec.push_back(drop_height + (rand()/(RAND_MAX + 1.0))*lambda_D/2);
        W_vec.push_back(-init_speed + (rand()/(RAND_MAX + 1.0))*init_speed/2);
        W_vec.push_back(-init_speed + (rand()/(RAND_MAX + 1.0))*init_speed/2);
        W_vec.push_back(-init_speed + (rand()/(RAND_MAX + 1.0))*init_speed/2);
    }
    void init_history_and_viz(){
        x_history.push_back(W_vec[0]/lambda_D);
        y_history.push_back(W_vec[1]/lambda_D);
        z_history.push_back(W_vec[2]/lambda_D);
        vx_history.push_back(W_vec[3]);
        vy_history.push_back(W_vec[4]);
        vz_history.push_back(W_vec[5]);
    }

    std::vector<double> f_der(std::vector<double> W_vec_f){
        std::vector<double> f; 
        f.push_back(W_vec_f[3]);
        f.push_back(W_vec_f[4]);
        f.push_back(W_vec_f[5]);
        //x,y components
        double radial = pow((pow(W_vec_f[0],2) + pow(W_vec_f[1],2)),0.5);
        double rad_acc_mag = (charge/m_D)*k_r_restore*abs(radial);

        f.push_back(rad_acc_mag*(W_vec_f[0]/radial) - (alpha_n*W_vec_f[3])/m_D + a_c[0] + (therm_coeff*dist(generator))/m_D);// - (alpha_i*pow(W_vec_f[3],2))/m_D TAKE OUT ION DRAG//drag, sheathe and coloumb force and ion drag force
        
        std::cout <<"X_a: " <<rad_acc_mag*(W_vec_f[0]/radial) << ","<<  - (alpha_n*W_vec_f[3])/m_D << ","<<a_c[0] <<","<< + (therm_coeff*dist(generator))/m_D << std::endl;
        f.push_back(rad_acc_mag*(W_vec_f[1]/radial) - (alpha_n*W_vec_f[4])/m_D + a_c[1] + (therm_coeff*dist(generator))/m_D);
        //z component
        if (W_vec_f[2] > z_se){
            std::cout <<"Z_a: " << -g_z<< ","<< a_c[2] << ","<<- (alpha_n*W_vec_f[5])/m_D <<","<< + (therm_coeff*dist(generator))/m_D << std::endl;
            f.push_back(-g_z + a_c[2] - (alpha_n*W_vec_f[5])/m_D  + (therm_coeff*dist(generator))/m_D);//drag, gravity, coloumb force and ion drag force
        }
        else{
            std::cout << "Z_a: " <<charge*Sheath_E << ","<< + ion_drag/m_D<< "," <<-g_z<< ","<< a_c[2] << ","<<- (alpha_n*W_vec_f[5])/m_D <<","<< + (therm_coeff*dist(generator))/m_D << std::endl;
            f.push_back(charge*Sheath_E + ion_drag/m_D - g_z  + a_c[2] - (alpha_n*W_vec_f[5])/m_D -  + (therm_coeff*dist(generator))/m_D); 
        };
        std::cout << f[5]<< std::endl;
        return f;
    }

    void step(double dt){
        std::cout << "step" << std::endl;
        std::vector<double> k1 = element_mul(f_der(W_vec),dt);
        std::vector<double> k2 = element_mul(f_der(element_add(W_vec,element_mul(k1,1/2))),dt);
        std::vector<double> k3 = element_mul(f_der(element_add(W_vec,element_mul(k2,1/2))),dt);
        std::vector<double> k4 = element_mul(f_der(element_add(W_vec,k3)),dt);
        W_vec = element_add(W_vec,element_mul(element_add(element_add(element_add(k1,element_mul(k2,2.0)),element_mul(k3,2.0)),k4),1.0/6.0));
    }
};

class Dust_Container{

    private:
    int dust_grain_max;

    public:
	std::vector<double> time_list;
    std::map<double, double> Y_list;
    std::map<double, double> Q_list;
    std::map<double, double> v_i_list;
    std::map<double, double> charge_list;
    std::map<double, double> ion_drag_list;
    std::vector<Dust_grain> Dust_grain_list;
    std::vector<std::pair<Dust_grain&,Dust_grain&>> combs_list;
    double v_squared_sum;
    double temperature;
    double z_se;
    std::vector<double> temperature_history;

	Dust_Container(int Dust_grain_max):dust_grain_max(Dust_grain_max),time_list{0},v_squared_sum{0}
    {
        //std::cout<< "start calc setup" << std::endl;
        calc_sheath();
        //std::cout<< "sheath done" << std::endl;
        calc_ion_vel();
        //std::cout<< "ion vel done" << std::endl;
        calc_dust_grain_charges();
        //std::cout<< "charge done" << std::endl;
        calc_ion_drag();
        //std::cout<< "ion drag done" << std::endl;
        create_dust_grains();
        //std::cout<< "create done" << std::endl;
    } 

    //CALCULATE THE SHEATH POTENTIAL AND ELECTRIC FIELD
    double f_Y_z_init(double Y_z_0){
        //"""RK4 derivative function"""
        return pow(2*(exp(Y_z_0) + pow(1-2*Y_z_0,0.5) - 2),0.5);
    }

    double step_Y_z_init(double Y_z_0){
        //"""RK4 step of size dz for Y_z_init"""
        double k1 = dz_norm*f_Y_z_init(Y_z_0);
        double k2 = dz_norm*f_Y_z_init(Y_z_0 + k1/2);
        double k3 = dz_norm*f_Y_z_init(Y_z_0 + k2/2);
        double k4 = dz_norm*f_Y_z_init(Y_z_0 + k3);
        double Y_z_1 = Y_z_0 + (k1 + 2*k2 + 2*k3 + k4)/6.0;
        return Y_z_1;
    }

    void calc_sheath(){
        double Y = V_DC;
        double edge_switch = 0;
        for(int i = 0; i < for_loop_max_int; i+= 1){
            //std::cout<< "sheath inside"<< ","<< i << std::endl;
            Y_list.insert(std::pair<double,double> (i*dz, Y));
            Q_list.insert(std::pair<double,double> (i*dz, f_Y_z_init(Y)));
            Y = step_Y_z_init(Y);
            if((Y > Y_epsilon) && (edge_switch == 0)){
                z_se = i*dz;
                edge_switch = 1;
            }
        }
    }

    //CALCULATE THE ION VELOCITY
    double v_i_z_calc(double Y, double z){
        if(z < z_se){
            return -pow(pow(v_B,2) - 2*i_charge*Y/m_i,0.5);          
        }
        else{
            return 0;
        }    
    }

    void calc_ion_vel(){
        for(int i = 0; i < for_loop_max_int; i+= 1){
            //std::cout << i << ","<< for_loop_max_int << ","<<Y_list.at(i*dz) << std::endl;
            v_i_list.insert(std::pair<double,double> (i*dz, v_i_z_calc(Y_list.at(i*dz),i*dz)));
        }
    }

    //CALCULATE THE DUST GRAIN CHARGE

    double f_x_Sheath(double x_n, double A, double B){
        return A*exp(B + x_n) - (2/(2*B - 1))*x_n - 1;
    }
    double f_x_first_derv_Sheath(double x_n, double A, double B){
        return A*exp(B + x_n) - (2/(2*B - 1));
    }
    double f_x_second_derv_Sheath(double x_n, double A, double B){
        return A*exp(B + x_n);
    }
    double calc_x_plus_Sheath(double x_n, double A, double B){
        double f_x_0 = f_x_Sheath(x_n,A,B);
        double f_x_1 = f_x_first_derv_Sheath(x_n,A,B);
        double f_x_2 = f_x_second_derv_Sheath(x_n,A,B);
        double x_plus = x_n - ((2*f_x_0*f_x_1)/(2.0*(pow(f_x_1,2))-f_x_0*f_x_2));
        return x_plus;
    }
    double find_phi_D_norm_OML_Sheath(double a_init,double z,double root_prec,double Phi_DC){
        double A = sqrt((8*k_b*T_e)/(pow(v_B,2)*M_PI*m_e));
        double a_n = a_init;
        double a_plus = calc_x_plus_Sheath(a_init,A,Phi_DC);

        while(std::abs(a_plus - a_n) > root_prec){
            a_n = a_plus;
            a_plus = calc_x_plus_Sheath(a_n,A,Phi_DC);
        }
        double W_0 = (a_n + a_plus)/2;
        return W_0;
    }
    double OML_charge(double z,double Phi_DC){
        //basically we have an equation that has the normalized surface potetnail in it and we rootfind for it, using the lambert W function form though tbh i agree with coppins this isnt actually necesary
        double Phi_D_norm = find_phi_D_norm_OML_Sheath(a_0_Sheath,z,root,Phi_DC);
        double phi_grain =  (k_b*T_e*Phi_D_norm)/(e_charge);// unormalize it
        return 4.0*M_PI*epsilon_0*grain_R*phi_grain;
    }
    void calc_dust_grain_charges(){
        for(int i = 0; i < for_loop_max_int; i+= 1){
            //std::cout << i << ","<< for_loop_max_int << ","<< Y_list.at(i*dz) << std::endl;
            charge_list.insert(std::pair<double,double> (i*dz, OML_charge(i*dz, Y_list.at(i*dz))));
        }
        //std::cout << charge_list.size() << std::endl;
    }

    //CALCULATE THE ION DRAG
    double ion_drag(double z,double charge_z,double v_i){//IT SHOULD DOUBEL BACK SURELY COS OF THE CHARGE FLIP?
        if(z > z_se){
            return 0;
        }
        else{
            //THE COLLECTION TERM DOMIANTES DRAMATICALLY
            double z_ion = (charge_z*i_charge/(4*M_PI*epsilon_0))*(1/(grain_R*k_b*T_e));
            double v_ti = -pow(k_b*T_i/m_i,0.5);
            double u = v_i/v_ti;
            double R = (charge_z*i_charge)/((4*M_PI*epsilon_0)*k_b*T_i*(1+pow(u,2)));
            double beta_bar = R/lambda_D;
            double LAMBDA = (beta_bar + 1)/(beta_bar + (grain_R/lambda_D));
            double F_i = -pow(2*M_PI,0.5)*pow(grain_R,2)*n_e0*m_i*pow(v_ti,2)*( pow(M_PI/2,0.5)*erf(pow(u/2,0.5))*(1+pow(u,2)+(1-pow(u,-2))*(1+2*z_ion*beta) + 4*pow(z_ion,2)*pow(beta,2)*pow(u,-2)*log(LAMBDA)) + pow(u,-1)*(1 + 2*z_ion*beta+pow(u,2) -  4*pow(z_ion,2)*pow(beta,2)*log(LAMBDA))*exp(-pow(u,2)/2));
            return F_i;
        }
    }

    void calc_ion_drag(){
        for(int i = 0; i < for_loop_max_int; i+= 1){
            ion_drag_list.insert( std::pair<double,double> (i*dz, ion_drag(i*dz,  charge_list.at(i*dz),   v_i_list.at(i*dz) ) ));
        }
    }
    //DONE

    void calc_temperature(){
        temperature = m_D*(v_squared_sum/Dust_grain_list.size())/(3*k_b);
        temperature_history.push_back(temperature);
    }

    void combs_list_produce(){
        combs_list.clear();
        for(int i = 0; i <  Dust_grain_list.size(); i++){
            for(int k = i+1;k < Dust_grain_list.size(); k++){ 
                std::pair<Dust_grain&,Dust_grain&> pair_comb (Dust_grain_list[i],Dust_grain_list[k]);
                combs_list.push_back(pair_comb);
            }
        }
    }
    void update_parameters_dust(int i){
        Dust_grain_list[i].charge = charge_list.lower_bound( Dust_grain_list[i].W_vec[2])->second;
        Dust_grain_list[i].Sheath_E = Q_list.lower_bound(Dust_grain_list[i].W_vec[2])->second;
        Dust_grain_list[i].v_i = v_i_list.lower_bound(Dust_grain_list[i].W_vec[2])->second;
        Dust_grain_list[i].ion_drag = ion_drag_list.lower_bound(Dust_grain_list[i].W_vec[2])->second;
    }
    
    void create_dust_grains(){
        //std::cout << dz << "," << root << std::endl;
        double r_01_mag;
        while(Dust_grain_list.size() < dust_grain_max){
            Dust_grain_list.push_back(Dust_grain());
            std::vector<double> pos_1 (Dust_grain_list.back().W_vec.begin(),Dust_grain_list.back().W_vec.begin() + 3);
            for(int v = 0; v <  Dust_grain_list.size() - 1; v++){
                std::vector<double> pos_0 (Dust_grain_list[v].W_vec.begin(),Dust_grain_list[v].W_vec.begin() + 3);
                r_01_mag = v_abs(element_add(pos_1, element_mul(pos_0,-1)));
                if (r_01_mag <= 3*grain_R){
                    Dust_grain_list.pop_back();
                    break;
                }
            }
            //std::cout << Dust_grain_list.size() << std::endl;
        } 
        for(int i = 0; i <  Dust_grain_list.size(); i++){
            Dust_grain_list[i].init_history_and_viz();
            update_parameters_dust(i);
            Dust_grain_list[i].z_se = z_se;//BIT OF SHIT CODE
            v_squared_sum += pow(Dust_grain_list[i].calc_speed(),2);
        }
        calc_temperature();
        v_squared_sum = 0;
        combs_list_produce();
    }

    void inter_particle_ca(){ 
        #pragma omp parallel for       
        for (int i = 0; i < combs_list.size(); i++){
            std::vector<double> force_c;
            std::vector<double> r_01;
            double r_01_mag;
            double p_mag;
            std::vector<double> force_c_pos_01{0,0,0};
            std::vector<double> force_c_pos_10{0,0,0};
            std::vector<double> pos_0 (combs_list[i].first.W_vec.begin(),combs_list[i].first.W_vec.begin() + 3);
            std::vector<double> pos_1 (combs_list[i].second.W_vec.begin(),combs_list[i].second.W_vec.begin() + 3);

            r_01 =  element_add(pos_1, element_mul(pos_0,-1));
            r_01_mag = v_abs(r_01);
            if (r_01_mag > coulomb_limit*lambda_D){
                continue;
            }

            std::vector<double>p_01{r_01[0], r_01[1]};
            p_mag = v_abs(p_01);
            force_c = element_mul(r_01,-((combs_list[i].first.charge*combs_list[i].second.charge)/(4*M_PI*epsilon_0))* exp((grain_R/lambda_D) - (r_01_mag/lambda_D)) * (1/(pow(r_01_mag,3)) + 1/(lambda_D*(pow(r_01_mag,2)))));
            
            if((pos_1[2] < (z_se - wake_safety_factor)) && (pos_1[2] > pos_0[2])){
                double M = combs_list[i].second.v_i/v_B;
                double z_plus = std::abs(r_01[2]) + p_mag*pow((pow(M,2)-1),0.5);
                double z_minus = std::abs(r_01[2]) - p_mag*pow((pow(M,2)-1),0.5);
                double A = cos((z_plus/lambda_D)/pow((pow(M,2)-1),0.5) - M_PI/4);
                double B = cos((z_minus/lambda_D)/pow((pow(M,2)-1),0.5) + M_PI/4);
                double C = sin((z_plus/lambda_D)/pow((pow(M,2)-1),0.5) - M_PI/4);
                double D = sin((z_minus/lambda_D)/pow((pow(M,2)-1),0.5) + M_PI/4);
                double F_p_far = -1*combs_list[i].first.charge*(2*combs_list[i].second.charge/(1 - pow(M,-2))) * (pow(lambda_D/(2*M_PI*p_mag),0.5))*(   (1/p_mag)*(-0.5)*((1/z_plus)*(A - 1/pow(2,0.5)) + (1/z_minus)*(B - 1/pow(2,0.5))) -pow(z_plus,-2)*pow((pow(M,2)-1),0.5)*(A - 1/pow(2,0.5)) - (1/z_plus)*(C/lambda_D) + pow(z_minus,-2)*pow((pow(M,2)-1),0.5)*(B - 1/pow(2,0.5)) + (1/z_minus)*(D/lambda_D) );
                double F_z_far = -1*combs_list[i].first.charge*(2*combs_list[i].second.charge/(1 - pow(M,-2))) * (pow(lambda_D/(2*M_PI*p_mag),0.5))*(  -pow(z_plus,-2)*(A - 1/pow(2,0.5)) - (1/z_plus)*(C*(1/lambda_D)/pow((pow(M,2)-1),0.5)) - pow(z_minus,-2)*(B - 1/pow(2,0.5)) - (1/z_minus)*(D*(1/lambda_D)/pow((pow(M,2)-1),0.5))      );
                force_c_pos_01 = {(p_01[0]/p_mag)*F_p_far, (p_01[1]/p_mag)*F_p_far, F_z_far};
            };
            if((pos_0[2] < (z_se - wake_safety_factor)) && (pos_0[2] > pos_1[2])){
                double M = combs_list[i].first.v_i/v_B;
                double z_plus = std::abs(r_01[2]) + p_mag*pow((pow(M,2)-1),0.5);
                double z_minus = std::abs(r_01[2]) - p_mag*pow((pow(M,2)-1),0.5);
                double A = cos((z_plus/lambda_D)/pow((pow(M,2)-1),0.5) - M_PI/4);
                double B = cos((z_minus/lambda_D)/pow((pow(M,2)-1),0.5) + M_PI/4);
                double C = sin((z_plus/lambda_D)/pow((pow(M,2)-1),0.5) - M_PI/4);
                double D = sin((z_minus/lambda_D)/pow((pow(M,2)-1),0.5) + M_PI/4);
                double F_p_far = -1*combs_list[i].second.charge*(2*combs_list[i].first.charge/(1 - pow(M,-2))) * (pow(lambda_D/(2*M_PI*p_mag),0.5))*(   (1/p_mag)*(-0.5)*((1/z_plus)*(A - 1/pow(2,0.5)) + (1/z_minus)*(B - 1/pow(2,0.5))) -pow(z_plus,-2)*pow((pow(M,2)-1),0.5)*(A - 1/pow(2,0.5)) - (1/z_plus)*(C/lambda_D) + pow(z_minus,-2)*pow((pow(M,2)-1),0.5)*(B - 1/pow(2,0.5)) + (1/z_minus)*(D/lambda_D) );
                double F_z_far = -1*combs_list[i].second.charge*(2*combs_list[i].first.charge/(1 - pow(M,-2))) * (pow(lambda_D/(2*M_PI*p_mag),0.5))*(  -pow(z_plus,-2)*(A - 1/pow(2,0.5)) - (1/z_plus)*(C*(1/lambda_D)/pow((pow(M,2)-1),0.5)) - pow(z_minus,-2)*(B - 1/pow(2,0.5)) - (1/z_minus)*(D*(1/lambda_D)/pow((pow(M,2)-1),0.5))      );
                force_c_pos_10 = {(-p_01[0]/p_mag)*F_p_far, (-p_01[1]/p_mag)*F_p_far, F_z_far};
            };

            combs_list[i].first.a_c = element_add(combs_list[i].first.a_c,element_add(element_mul(force_c,1/m_D), element_mul(force_c_pos_01,1/m_D)));
            combs_list[i].second.a_c = element_add(combs_list[i].second.a_c,element_add(element_mul(force_c,1/-m_D) , element_mul(force_c_pos_10,1/m_D)));
        }
    }

    void next_frame(double dt){
    //The big daddy of the code, this functions loops over advancing the simulation by a time step dt
        //cout<< Dust_grain_list.size() << endl;
        inter_particle_ca();
		time_list.push_back(time_list.back() + dt);
        
        #pragma omp parallel for
        for (int i = 0; i < Dust_grain_list.size(); i++){
            //advance via the rk4 and add the predicted postiosn to the momentary list
            Dust_grain_list[i].step(dt);
            Dust_grain_list[i].a_c = {0,0,0};//NOT SURE ABOUT THIS
            Dust_grain_list[i].x_history.push_back(Dust_grain_list[i].W_vec[0]/lambda_D);
            Dust_grain_list[i].y_history.push_back( Dust_grain_list[i].W_vec[1]/lambda_D);
            Dust_grain_list[i].z_history.push_back(Dust_grain_list[i].W_vec[2]/lambda_D);
            update_parameters_dust(i);
            v_squared_sum += pow(Dust_grain_list[i].calc_speed() ,2);
        };
        calc_temperature();
        v_squared_sum = 0;
    }
};

void write_csv(std::string filename, std::vector<std::pair<std::string, std::vector<double>>> dataset){
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Send column names to the stream
    for(int j = 0; j < dataset.size(); ++j){
        myFile << dataset.at(j).first;//get the column title
        if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
    }
    myFile << std::endl;//next line
    
    // Send data to the stream
    for(int i = 0; i < dataset.front().second.size(); ++i){//loop through rows
        for(int j = 0; j < dataset.size(); ++j){//loop through columbs
            //cout << (dataset.back().second.size() - dataset.at(j).second.size()) << endl;
            if(i <  dataset.at(j).second.size()){
                myFile << dataset.at(j).second.at(i);
            } 
            if(j != dataset.size() - 1){
                myFile << ","; // No comma at end of line
            } 
        }
        myFile << std::endl;
    }    
    // Close the file
    myFile.close();
}

int main(){
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::vector<double> speed_list;
    std::cout << "Running..." << std::endl;

    Dust_Container Dusty_plasma_crystal = Dust_Container(dust_grain_max_input);

    std::cout << "INJECT DUST" << std::endl;

    while (Dusty_plasma_crystal.time_list.back() < time_limit){
        Dusty_plasma_crystal.next_frame(dt_a);
    }

    std::cout << "Simulation finished" << std::endl;

    /////////////////////

    std::string filename = "HPC_Data/Finite_Dust_grain_max_" + std::to_string(dust_grain_max_input);
    filename += "_Final_Temperature_" + std::to_string(Dusty_plasma_crystal.temperature);
    filename += "_frames_" + std::to_string(Dusty_plasma_crystal.time_list.size());
    filename += ".csv";

    std::vector<std::pair<std::string,std::vector<double>>> vals;

    for(int i = 0 ; i < Dusty_plasma_crystal.Dust_grain_list.size(); i++){
        vals.push_back(make_pair("X_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].x_history));
        vals.push_back(make_pair("Y_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].y_history));
        vals.push_back(make_pair("Z_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].z_history));
        vals.push_back(make_pair("VX_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].vx_history));
        vals.push_back(make_pair("VY_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].vy_history));
        vals.push_back(make_pair("VZ_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].vz_history));
        speed_list.push_back(Dusty_plasma_crystal.Dust_grain_list[i].calc_speed());
    };
    vals.push_back(make_pair("Time_list", Dusty_plasma_crystal.time_list));
    vals.push_back(make_pair("Temperature_list", Dusty_plasma_crystal.temperature_history));
    vals.push_back(make_pair("Speed_list", speed_list));

    write_csv(filename, vals);

    std::cout << "FILENAME:" << filename << std::endl;

    //////////////

    std::string filename_end = "HPC_Data/Final_Finite_Dust_grain_max_" + std::to_string(dust_grain_max_input);
    filename_end += "_Final_Temperature_" + std::to_string(Dusty_plasma_crystal.temperature);
    filename_end += "_frames_" + std::to_string(Dusty_plasma_crystal.time_list.size());
    filename_end += ".csv";

    std::vector<std::pair<std::string, std::vector<double> >> vals_end;
    std::vector<double> final_x;
    std::vector<double> final_y;
    std::vector<double> final_z;
    std::vector<double> final_vx;
    std::vector<double> final_vy;
    std::vector<double> final_vz;

    for(int i = 0 ; i < Dusty_plasma_crystal.Dust_grain_list.size(); i++){
        final_x.push_back(Dusty_plasma_crystal.Dust_grain_list[i].W_vec[0]);
        final_y.push_back(Dusty_plasma_crystal.Dust_grain_list[i].W_vec[1]);
        final_z.push_back(Dusty_plasma_crystal.Dust_grain_list[i].W_vec[2]);
        final_vx.push_back(Dusty_plasma_crystal.Dust_grain_list[i].W_vec[3]);
        final_vy.push_back(Dusty_plasma_crystal.Dust_grain_list[i].W_vec[4]);
        final_vz.push_back(Dusty_plasma_crystal.Dust_grain_list[i].W_vec[5]);
    };

    vals_end.push_back(make_pair("X",final_x));
    vals_end.push_back(make_pair("Y",final_y));
    vals_end.push_back(make_pair("Z",final_z));
    vals_end.push_back(make_pair("V_X",final_vx));
    vals_end.push_back(make_pair("V_Y",final_vy));
    vals_end.push_back(make_pair("V_Z",final_vz));

    write_csv(filename_end, vals_end);
    std::cout << "FILENAME FINAL:" << filename_end << std::endl;

    std::cout << "done" << std::endl;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

    return 0;
}

