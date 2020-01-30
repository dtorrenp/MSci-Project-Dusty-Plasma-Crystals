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
const int dust_grain_max_input = 100; //dust grain max number
const double dt_a = 1.0e-4;
const double time_limit = 1.0e-1;
const double frame_req = 5;

//CONSTANTS TO FUCK ABOUT WITH
const double n_e0 = 1.0e15; //electron number density (in bulk plasma)
const double n_i0 = 1.0e15; //ion number density (in bulk plasma)
const double n_n0 = 1.0e24; //electron number density (in bulk plasma)
const double Z = 1; //ion atomic number
const double Z_n = 18; //neutrals atomic number (Ar)
const double grain_R = 7*1e-6; //dust grain radius
const double dust_grain_density = 1.49*1e3; //dust density
const double phi_wall_z = -100.0; //vertical wall potential [volts]
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
const double beta = T_i/T_e;
const double lambda_de = pow(((epsilon_0*k_b*T_e)/(n_e0*(pow(e_charge,2)))),0.5);
const double lambda_di = pow(((epsilon_0*k_b*T_i)/(n_i0*(pow(e_charge,2)))),0.5);
const double lambda_D = pow((1/(1/(pow(lambda_de,2)) + 1/(pow(lambda_di,2)))),0.5);//get lambda_D
const double wake_potential_below = 1*lambda_D;
const double wake_charge_multiplier = 1.0;
const double drop_height = 9.8*lambda_D;//drop particles from this height, low so that we dont waste computational time on calculations as its falling and not interacting with sheathe
const double container_length = 10.0*lambda_D; //container radius
const double coulomb_limit = 5;

//related to finding the charge of the dust grains
const double dz = lambda_D/1000;
const double root = dz/10;//preciscion of root finding method used to get dust charge

const double lower_lim_z = 7.0*lambda_D;
const double upper_lim_z = 15.0*lambda_D;
const double a_0_OML = -0.5;
const double a_0_Sheath = -2.5;//intial guess for halley's method

//const double container_dust_dist_creation = sqrt(container_radius/(2*lambda_D));
const double z_se = 10.0*lambda_D; //distance of vertical sheath from bottom of container
const double r_se = 25.0*lambda_D; //distance of radial 'sheath' from wall of container
const double k_z_restore = -2.0*phi_wall_z/pow(z_se,2); //WIERD MINUS SIGN TO ACCOUNT FOR FACT THAT K MUST BE POSITIVE WE THINK BUT NEED TO COME BACK TO THIS
const double v_B = pow((3*k_b*T_e/m_i),0.5);
const double v_Tn = pow((3*k_b*T_i/m_n),0.5);//thermal termperature of the neutrals
const double alpha_n = (4/3)*M_PI*pow(grain_R,2)*m_n*n_n0*v_Tn;//maybe 8/3
const double alpha_i = M_PI*pow(grain_R,2)*m_i*n_i0;
const double therm_coeff = sqrt(2*k_b*T_i*alpha_n);
const double therm_coeff_i = sqrt(2*k_b*T_i*alpha_i);

std::vector<double> element_mul(const std::vector<double>& a,double cst){
    std::vector<double> c;
    for(int i=0; i<a.size(); i++){
        c.push_back(a[i]*cst);
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

void print_vector(const std::vector<double>& a){
    std::cout << "[";
    for(int i=0; i<a.size(); i++){
        std::cout << a[i];
        if(i != a.size() - 1){
                std::cout << ","; // No comma at end of line
            } 
    }
    std::cout << "]"<< std::endl;
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
    double wake_charge;
    double v_i_z;
    double charge;

    Dust_grain():a_c{0,0,0},v_i_z(0),
    generator(std::default_random_engine(clock())),dist(std::normal_distribution<double>(0.0,sqrt((k_b*T_i)/m_D)))
    {
        prod_W_vec();        
    }

    double calc_speed(){
        std::vector<double> vel (W_vec.end() - 3,W_vec.end());       
        return v_abs(vel);
    }
    
    void prod_W_vec(){
        //when creating new dust grains gives them random x,y positions so the dust grains dont just stack on 0,0//
        W_vec.push_back((rand()/(RAND_MAX + 1.0))*container_length - container_length/2);//create random number centred about 0 with width container_radius/3
        W_vec.push_back((rand()/(RAND_MAX + 1.0))*container_length - container_length/2);
        W_vec.push_back(drop_height + (rand()/(RAND_MAX + 1.0))*lambda_D/2);
        W_vec.push_back(-init_speed + (rand()/(RAND_MAX + 1.0))*init_speed/2);
        W_vec.push_back(-init_speed + (rand()/(RAND_MAX + 1.0))*init_speed/2);
        W_vec.push_back(-init_speed + (rand()/(RAND_MAX + 1.0))*init_speed/2);
        x_history.push_back(W_vec[0]/lambda_D);
        y_history.push_back(W_vec[1]/lambda_D);
        z_history.push_back(W_vec[2]/lambda_D);
    }

    std::vector<double> f_der(std::vector<double> W_vec_f){
        std::vector<double> f; 
        f.push_back(W_vec_f[3]);
        f.push_back(W_vec_f[4]);
        f.push_back(W_vec_f[5]);
        //x,y components
        f.push_back( - (alpha_n*W_vec_f[3])/m_D + a_c[0] + (therm_coeff*dist(generator))/m_D + (therm_coeff_i*dist(generator))/m_D);//drag, sheathe and coloumb force and ion drag force
        f.push_back( - (alpha_n*W_vec_f[4])/m_D + a_c[1] + (therm_coeff*dist(generator))/m_D + (therm_coeff_i*dist(generator))/m_D);
        if (W_vec_f[2] > z_se){
            f.push_back(- g_z - (alpha_n*W_vec_f[5])/m_D + a_c[2] + (therm_coeff*dist(generator))/m_D);//drag, gravity, coloumb force and ion drag force
        }
        else{
            v_i_z = pow((pow(v_B,2) + (i_charge*k_z_restore*pow((W_vec_f[2] - z_se),2))/m_i -2.0*g_z*(W_vec_f[2] - z_se)),0.5);
            f.push_back((charge/m_D)*k_z_restore*(W_vec_f[2] - z_se) - g_z + a_c[2] - (alpha_n*W_vec_f[5])/m_D - (alpha_i*pow(v_i_z - W_vec_f[5],2))/m_D + (therm_coeff*dist(generator))/m_D + (therm_coeff_i*dist(generator))/m_D);//drag, sheathe, gravity, coloumb force and ion drag force;
        };
        return f;
    }

    void step(double dt){
        std::vector<double> k1 = element_mul(f_der(W_vec),dt);
        std::vector<double> k2 = element_mul(f_der(element_add(W_vec,element_mul(k1,1/2))),dt);
        std::vector<double> k3 = element_mul(f_der(element_add(W_vec,element_mul(k2,1/2))),dt);
        std::vector<double> k4 = element_mul(f_der(element_add(W_vec,k3)),dt);
        std::vector<double> W_vec_temp = element_add(W_vec,element_mul(element_add(element_add(element_add(k1,element_mul(k2,2.0)),element_mul(k3,2.0)),k4),1.0/6.0));
        for(int i = 0; i < 2; i++){
            if(std::abs(W_vec_temp[i]) > container_length/2){
                if(W_vec_temp[i]> 0){
                    W_vec[i] = W_vec_temp[i] - container_length;
                }
                else{
                    W_vec[i] = W_vec_temp[i] + container_length;
                };
            }
            else{
                W_vec[i] = W_vec_temp[i];
            }
            W_vec[i+3] = W_vec_temp[i+3];
        };
        W_vec[2] = W_vec_temp[2];
        W_vec[5] = W_vec_temp[5];
        //std::cout << W_vec[0] << "," << W_vec[1] << "," << W_vec[2] << std::endl;
    }
};

class Dust_Container{

    private:
    double q_D;
    int dust_grain_max;

    public:
	std::vector<double> time_list;
    std::map<double, double> charge_list;
    std::vector<Dust_grain> Dust_grain_list;
    std::vector<std::pair<Dust_grain&,Dust_grain&>> combs_list;
    double v_squared_sum;
    double temperature;
    std::vector<double> temperature_history;

	Dust_Container(int Dust_grain_max):dust_grain_max(Dust_grain_max),time_list{0},v_squared_sum{0}
    {
        calc_dust_grain_charges(dz, lower_lim_z,upper_lim_z);
        create_dust_grains();
    } 

    double f_x_OML(double x_n){
        return sqrt(beta)*(1-x_n/beta) -exp(x_n);
    }
    double f_x_first_derv_OML(double x_n){
        return -sqrt(beta)*(1/beta) -exp(x_n);
    }
    double f_x_second_derv_OML(double x_n){
        return -exp(x_n);
    }
    double calc_x_plus_OML(double x_n){
        double f_x_0 = f_x_OML(x_n);
        double f_x_1 = f_x_first_derv_OML(x_n);
        double f_x_2 = f_x_second_derv_OML(x_n);
        double x_plus = x_n - ((2*f_x_0*f_x_1)/(2.0*(pow(f_x_1,2))-f_x_0*f_x_2));
        return x_plus;
    }
    double  find_phi_D_norm_OML(double a_init,double root_prec){
        double a_n = a_init;
        double a_plus = calc_x_plus_OML(a_init);//do i need to give type if product of a thing?

        while(abs(a_plus - a_n) > root_prec){
            a_n = a_plus;
            a_plus = calc_x_plus_OML(a_n);
        }
        double W_0 = (a_n + a_plus)/2;
        return W_0;
    }

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
    double  find_phi_D_norm_OML_Sheath(double a_init,double z,double root_prec){
        double A = sqrt((8*k_b*T_e)/(pow(v_B,2)*M_PI*m_e));
        double B = ((2*e_charge*phi_wall_z)/(k_b*T_e))*pow((z/z_se - 1),2);
        double a_n = a_init;
        double a_plus = calc_x_plus_Sheath(a_init,A,B);

        while(abs(a_plus - a_n) > root_prec){
            a_n = a_plus;
            a_plus = calc_x_plus_Sheath(a_n,A,B);
        }
        double W_0 = (a_n + a_plus)/2;
        return W_0;
    }
    double OML_charge(double z){
        //basically we have an equation that has the normalized surface potetnail in it and we rootfind for it, using the lambert W function form though tbh i agree with coppins this isnt actually necesary
        double Phi_D_norm;
        if (z > z_se){
            Phi_D_norm = find_phi_D_norm_OML(a_0_OML,root);
        }
        else{
            Phi_D_norm = find_phi_D_norm_OML_Sheath(a_0_Sheath,z,root);
        }
        double phi_grain =  (k_b*T_e*Phi_D_norm)/(e_charge);// unormalize it
        return 4.0*M_PI*epsilon_0*grain_R*phi_grain;
    }
    void calc_dust_grain_charges(double d_z, double low_z, double high_z){
        for(double i = low_z; i < high_z; i+= d_z){
            charge_list.insert( std::pair<double,double> (i, OML_charge(i)) );
        }
        std::cout << charge_list.size() << std::endl;
    }

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
    
    void create_dust_grains(){
        std::cout << dz << "," << root << std::endl;
        double r_01_mag;
        while(Dust_grain_list.size() < dust_grain_max){
            Dust_grain_list.push_back(Dust_grain());
            std::vector<double> pos_1 (Dust_grain_list.back().W_vec.begin(),Dust_grain_list.back().W_vec.begin() + 3);
            for(int v = 0; v <  Dust_grain_list.size() - 1; v++){
                std::vector<double> pos_0 (Dust_grain_list[v].W_vec.begin(),Dust_grain_list[v].W_vec.begin() + 3);
                r_01_mag = v_abs(element_add(pos_1, element_mul(pos_0,-1)));
                if (r_01_mag <= 2*grain_R){
                    Dust_grain_list.pop_back();
                    break;
                }
            }
            std::cout << Dust_grain_list.size() << std::endl;
        } 
        for(int i = 0; i <  Dust_grain_list.size(); i++){
            Dust_grain_list[i].charge = charge_list.lower_bound( Dust_grain_list[i].W_vec[2])->second;
            Dust_grain_list[i].wake_charge = abs(Dust_grain_list[i].charge)*wake_charge_multiplier;
            v_squared_sum += pow(Dust_grain_list[i].calc_speed(),2);
        }
        calc_temperature();
        v_squared_sum = 0;
        combs_list_produce();
    }

    void inter_particle_ca(){   
        #pragma omp parallel for
        for (int i = 0; i < combs_list.size(); i++){
            double wake_charge; 
            std::vector<double> force_c;
            std::vector<double> r_01;
            std::vector<double> r_01_pos{0,0,0};
            std::vector<double> r_10_pos{0,0,0};
            double r_01_mag;
            double r_01_pos_mag;
            double r_10_pos_mag;
            std::vector<double> force_c_pos_01{0,0,0};
            std::vector<double> force_c_pos_10{0,0,0};
            std::vector<double> pos_0 (combs_list[i].first.W_vec.begin(),combs_list[i].first.W_vec.begin() + 3);
            std::vector<double> pos_1 (combs_list[i].second.W_vec.begin(),combs_list[i].second.W_vec.begin() + 3);
            
            //std::cout << "heyy" << std::endl;

            r_01 =  element_add(pos_1, element_mul(pos_0,-1));
            
            if(std::abs(r_01[0]) > container_length/2){

                if(pos_1[0] > 0){
                    r_01[0] = pos_1[0] - pos_0[0] - container_length;
                }
                else{
                    r_01[0] = pos_1[0] - pos_0[0] + container_length;
                }
            }

            if(std::abs(r_01[1]) > container_length/2){

                if(pos_1[1] > 0){
                    r_01[1] = pos_1[1] - pos_0[1] - container_length;
                }
                else{
                    r_01[1] = pos_1[1] - pos_0[1] + container_length;
                }
            }

            r_01_mag = v_abs(r_01);

            if (r_01_mag > coulomb_limit*lambda_D){
                continue;//NOT SURE IF THIS IS NECESSARY
            }

            force_c = element_mul(r_01,-((combs_list[i].first.charge*combs_list[i].second.charge)/(4*M_PI*epsilon_0))* exp((grain_R/lambda_D) - (r_01_mag/lambda_D)) * (1/(pow(r_01_mag,3)) + 1/(lambda_D*(pow(r_01_mag,2)))));

            if(pos_1[2] < z_se){
                r_01_pos[0] = r_01[0];
                r_01_pos[1] = r_01[1];
                r_01_pos[2] = r_01[2] - wake_potential_below;
                r_01_pos_mag = v_abs(r_01_pos);
                wake_charge = std::abs(combs_list[i].first.v_i_z/v_B)*combs_list[i].second.wake_charge;
                force_c_pos_01 = element_mul(r_01_pos,-((combs_list[i].first.charge*wake_charge)/(4*M_PI*epsilon_0))* exp((grain_R/lambda_D) - (r_01_pos_mag/lambda_D)) * (1/(pow(r_01_pos_mag,3)) + 1/(lambda_D*(pow(r_01_pos_mag,2)))));
            };

            if(pos_0[2] < z_se){
                r_10_pos[0] = -r_01[0];
                r_10_pos[1] = -r_01[1];
                r_10_pos[2] = -r_01[2] - wake_potential_below;
                r_10_pos_mag = v_abs(r_10_pos);
                wake_charge = std::abs(combs_list[i].first.v_i_z/v_B)*combs_list[i].first.wake_charge;
                force_c_pos_10 = element_mul(r_10_pos,-((combs_list[i].second.charge*wake_charge)/(4*M_PI*epsilon_0))*exp((grain_R/lambda_D) - (r_10_pos_mag/lambda_D)) * (1/(pow(r_10_pos_mag,3)) + 1/(lambda_D*(pow(r_10_pos_mag,2)))));
            };
            combs_list[i].first.a_c = element_add(combs_list[i].first.a_c,element_add(element_mul(force_c,1/m_D), element_mul(force_c_pos_01,1/m_D)));
            combs_list[i].second.a_c = element_add(combs_list[i].second.a_c,element_add(element_mul(force_c,1/-m_D) , element_mul(force_c_pos_10,1/m_D)));  
        }
    }

    void next_frame(double dt){
    //The big daddy of the code, this functions loops over advancing the simulation by a time step dt
        //std::cout << "frame" << std::endl;
        inter_particle_ca();
		time_list.push_back(time_list.back() + dt);

        #pragma omp parallel for
        for (int i = 0; i < Dust_grain_list.size(); i++){
            Dust_grain_list[i].step(dt);
            Dust_grain_list[i].a_c = {0,0,0};
            Dust_grain_list[i].x_history.push_back(Dust_grain_list[i].W_vec[0]/lambda_D);
            Dust_grain_list[i].y_history.push_back( Dust_grain_list[i].W_vec[1]/lambda_D);
            Dust_grain_list[i].z_history.push_back(Dust_grain_list[i].W_vec[2]/lambda_D);
            v_squared_sum += pow(Dust_grain_list[i].calc_speed(),2);

            Dust_grain_list[i].charge = charge_list.lower_bound(Dust_grain_list[i].W_vec[2])->second;
            Dust_grain_list[i].wake_charge = abs(Dust_grain_list[i].charge)*wake_charge_multiplier;
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

    while (Dusty_plasma_crystal.time_list.back() < time_limit){
        Dusty_plasma_crystal.next_frame(dt_a);
    }

    std::cout << "Simulation finished" << std::endl;

    /////////////////////

    std::string filename = "HPC_Data/Periodic_Dust_grain_max_" + std::to_string(dust_grain_max_input);
    filename += "_Final_Temperature_" + std::to_string(Dusty_plasma_crystal.temperature);
    filename += "_frames_" + std::to_string(Dusty_plasma_crystal.time_list.size());
    filename += ".csv";
    std::vector<std::pair<std::string,std::vector<double>>> vals;

    for(int i = 0 ; i < Dusty_plasma_crystal.Dust_grain_list.size(); i++){
        vals.push_back(make_pair("X_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].x_history));
        vals.push_back(make_pair("Y_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].y_history));
        vals.push_back(make_pair("Z_" + std::to_string(i), Dusty_plasma_crystal.Dust_grain_list[i].z_history));
        speed_list.push_back(Dusty_plasma_crystal.Dust_grain_list[i].calc_speed());
    };
    vals.push_back(make_pair("Time_list", Dusty_plasma_crystal.time_list));
    vals.push_back(make_pair("Temperature_list", Dusty_plasma_crystal.temperature_history));
    vals.push_back(make_pair("Speed_list", speed_list));

    write_csv(filename, vals);
    std::cout << "FILENAME:" << filename << std::endl;


    //////////////

    std::string filename_end = "HPC_Data/Final_Periodic_Dust_grain_max_" + std::to_string(dust_grain_max_input);
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
    std::cout << "FILENAME_END:" << filename_end << std::endl;

    /////////////

    std::cout << "done" << std::endl;

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;

    return 0;
}
