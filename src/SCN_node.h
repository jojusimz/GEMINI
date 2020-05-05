
#ifndef SCN_NODE_H
#define SCN_NODE_H

#include <string>
#include <utility>      // std::pair
#include <sstream>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <iterator>
#include <functional>
#include <cmath>
#include <iomanip>
#include <limits>   // double precisio


using namespace std;

class SCN_node
{
public:
    SCN_node(int scn_type,double dl, double tfactor, float er,float sigma_ey, int id, int xx, int yy, int zz,bool special_node_flag,bool is_PEC, bool PML_check,int PML_thickness=0,double Refn_fctr=1,int conduct_profile=0, int Lx=-1, int Ly=-1,int Lz=-1);
    SCN_node(const SCN_node& copy_SCN_node);
    SCN_node();
    SCN_node & operator=(const SCN_node& copy_SCN_node);
    void excitation (int port_id=0,double V_i=0);            //sets port_id and incidence voltage
    friend ostream& operator <<(std::ostream& out,SCN_node& node);            // output streaming
    void set_neighbour(SCN_node *node);
    void print_neighbours();
    void scattering();
    void scattering_PML();
    void scattering_PML2();
    void scattering_EPML();
    void scattering_PML_df(bool check);
    void scattering_source_node(double v_inc, double Rs, double & Vy_prime);
    double compute_v_series_k(int k);
    double compute_v_dummy_kj(double cndtvty_i,double cndtvty_j,double cndtvty_k, double i_k_pml,double Ik, int axis, double v_dummy_k_n_1);
    double compute_i_k_pml(int k,double v_series_k);
    double compute_i_shunt_pml_a(int port_id, double H, double cndtvty);
    double compute_i_shunt_pml(int j);
    double compute_Vj_pml(int j);
    double compute_alpha(int j);
    double compute_Ik(int k);
    double compute_Vj(int j);
    void compute_source_voltage_y(const double v_inc,int source_res,double &Vz,double &Vz_prime);
    double get_reflected(int port_id);
    void  set_reflected( int port_id,double Vr);
    void set_all_to_one();
    void print_node_voltages();
    double get_incidence( int port_id);
    int get_iD();
    double get_cndtvty_x();
    int get_coord(int i);
    double get_cndtvty_y();
    double get_cndtvty_z();
    float get_alpha_x(){return alpha_x;}
    float get_alpha_y(){return alpha_y;}
    float get_alpha_z(){return alpha_z;}
    void  set_incidence(int port_id,double Vr);
    bool check_special_node();
    void set_special_SCN_node(bool flag);
    bool  is_PML_SCN();
    bool  is_PEC();
    void set_PEC();
    void set_PML (int Lz, double Refn_fctr, int c_profile , int scn_typ, int num_PML);
    void set_soft_source_node();
    double get_G(int i);
    double get_R(int i);
    double get_Y(int i);
    double get_Z(int i);
    vector<double> get_SCN_inc_ports();
    virtual ~SCN_node();
private:
    int SCN_id;
    int scn_type;
    bool PEC;
    double SCN_dl;
    int conduct_prof;
    bool soft_source_node;
    vector < double > sigma_e;              // sigma_e[0]: sigma_ex ; sigma_e[1] : sigma_ey;
    vector < double> sigma_h;
    bool is_PML_node;
    vector < double > G;                    //we assume length in all dimension is equal dx=dy=dz
    vector < double > R;
    vector < double > Y;
    vector < double > Z;
    double Refn_fctr;
    int PML_Lx;
    int PML_Ly;
    int PML_Lz;
    vector< pair <double,double> > SCN_ports;     // SCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.
    // PML DF component
    double cndtvty_x;
    double cndtvty_y;
    double cndtvty_z;

    float alpha_x;
    float alpha_y;
    float alpha_z;
    vector < double >  V_inc_n_2;
    vector < double >  V_inc_n_1;
    vector < double >  V_j_n_1;
    vector < double >  V_j_pml_n_1;
    vector < double >  V_j_tot_n_1;
    vector < double >  I_j_sh_pml_n_1;
    vector < double >  I_j_sh_pml_a_n_1;
    vector < double >  I_j_sh_tot_a_n_1;
    vector < double > I_k_n_1;
    vector < double > I_k_pml_n_1;
    vector < double > I_k_tot_n_1;
    vector < double > V_k_series_n_1;
    vector < double > V_kj_dummy_n_1;
    vector < double > V_kj_dummy_tot_n_1;
    //vector<SCN_node*> SCN_neighbours;
    vector<int> SCN_coord;
    bool special_node;
    bool special_node_set_flag;
};


#endif
