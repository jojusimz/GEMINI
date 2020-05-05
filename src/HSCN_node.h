

#ifndef HSCN_NODE_H
#define HSCN_NODE_H

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
#include <limits>   //   double precisio


using namespace std;

//int return_coordinates(  double width,   double height,   double length,int npml,   double dl,   double x,   double y,   double z);

class HSCN_node
{
public:
    //CONSTRUCTOR
    HSCN_node(int sim_identifier,float dl,  double dt, int id, int xx, int yy, int zz,bool special_node_flag,int medium_type,bool is_PEC, bool PML_check,
              int PML_thickness=0,  double Refn_fctr=1,int conduct_profile=0, int Lx=-1, int Ly=-1,int Lz=-1);
    HSCN_node(const HSCN_node& copy_HSCN_node);
    HSCN_node();
    HSCN_node & operator=(const HSCN_node& copy_HSCN_node);

    //EXCITATION
    void excitation (int port_id=0,  double V_i=0);            //sets port_id and incidence voltage

    //SCATTERING
    void scattering_standard(vector<float> impedance_y, vector<float> conductance_G, vector<float> Rm);
    void scattering(vector<float> impedance_y, vector<float> conductance_G, vector<float> Rm);

    //COMPUTE VOLTAGES
    void compute_Vj(vector<float> impedance_y, vector <float> conductance_g ,  double & Vx,   double &Vy ,   double &Vz,
                       double & Vx_n_1,   double &Vy_n_1 ,   double &Vz_n_1,   double & Vx_n_2,   double &Vy_n_2 ,   double &Vz_n_2);

    void compute_standard_Vj( vector<float> impedance_y, vector <float> conductance_g ,   double & Vx,   double &Vy ,  double &Vz);

    void compute_complex_Vj(vector<float> impedance_y, vector <float> conductance_g ,   double & Vx,   double &Vy ,   double &Vz);

    void compute_V_xy_yx(float y_xy, float Rmz,   double & Vyx,   double &Vxy);

    void compute_V_xz_zx(float y_xz, float Rmy,   double & Vxz,   double &Vzx);

    void compute_V_yz_zy(float y_zy, float Rmx,   double &Vyz ,   double &Vzy);

    void compute_standard_V_xy_yx(float y_xy, float Rmz,   double & Vyx,   double &Vxy);

    void compute_standard_V_xz_zx(float y_xz, float Rmy,   double & Vxz,   double &Vzx);

    void compute_standard_V_yz_zy(float y_zy, float Rmx,   double &Vyz ,   double &Vzy);

    void update_voltages( );

    //PRINT
    void print_node_voltages();
    void print_neighbours();
    friend ostream& operator <<(std::ostream& out,HSCN_node& node);            // output streaming

    //IS A
    bool  is_PML_HSCN();
    bool  is_PEC();
    bool check_special_node();

    //GETS
      double get_incidence( int port_id);
    int get_iD();
    int get_coord(int i);
      double get_cndtvty_x();
      double get_cndtvty_y();
      double get_cndtvty_z();
      double get_Vn(int i);
      double get_reflected(int port_id);
    int  get_medium_type();

    //SETS
     void set_PEC();
     void set_special_HSCN_node(bool flag);
     void  set_incidence(int port_id,  double Vr);
     void  set_reflected( int port_id,  double Vr);
     void set_all_to_one();
     void set_neighbour(HSCN_node *node);
     void set_material_type(int media_type);

    virtual ~HSCN_node();
private:

    int HSCN_node_sim_identifier;
    int media_identifier;
    bool is_PML_node;
    unsigned int conduct_prof;
    bool PEC;
    float HSCN_dl;
    int HSCN_id;
    int PML_Lx;
    int PML_Ly;
    int PML_Lz;
    vector< pair <  double,  double> > HSCN_ports;     // HSCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.
      double cndtvty_x;
      double cndtvty_y;
      double cndtvty_z;
    vector <   double >  V_inc_n_2;
    vector <   double >  V_inc_n_1;
    vector <   double >  V_j_tot_n;
    vector <   double >  V_j_tot_n_1;
    vector <   double >  V_j_tot_n_2;
    vector <   double >  V_IZ_tot_n;
    vector <   double >  V_IZ_tot_n_1;
    vector <   double >  V_IZ_tot_n_2;

    vector<int> HSCN_coord;
    bool special_node;
    bool special_node_set_flag;
};

#endif
