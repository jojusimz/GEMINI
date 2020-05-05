
#ifndef mesh_handler_SCN_H
#define mesh_handler_SCN_H

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



#include "SCN_node.h"
#include "gemini_utility.h"

using namespace std;


class mesh_handler_SCN
{
public:
    mesh_handler_SCN(int simtype,double width, double height, double length, double dl,double tfactor,float er,float sigma_ey, int scn_type ,int n_PMl=0,int conduct_prof=0,double Refn_factor=1);
    mesh_handler_SCN(double width, double height, double length,int simtype, double dl,double tfactor,float er,float sigma_ey,  int scn_type, int n_PMl=0,int conduct_prof=0,double Refn_factor=1);
    mesh_handler_SCN(double width, double height,int simtype, double dl,double tfactor,float er,float sigma_ey,  int scn_type, int n_PMl=0,int conduct_prof=0,double Refn_factor=1);
    mesh_handler_SCN(const mesh_handler_SCN& copy_mesh_handler_SCN);
    mesh_handler_SCN & operator=(const mesh_handler_SCN & Other);
    void Re_mesh_mesh_handler_SCN(double width, double height, double length, double dl);              // creates a new rectangular waveguide and re-meshes using the input parameters.
    void set_neighbours();           // finds the neighbours for each node and points them to each other.
    void print_all_neighbours();
    void TLM_simulation1(int dt_total, double tfactor, int output_node, int excited_node, double Zz1, double Zz2, double Zx1, double Zx2, double Zy1, double Zy2, f_gaussian func , double excitation_type=-1,int srce_res=0, bool scn_xy=false) ;     //simulates wave propagation using structured TLM
    void TLM_simulation_2D_SCN(int dt_total,double tfactor, int output_node,int excited_node,double Zz1, double Zz2, double Zy1, double Zy2,double Zx1, double Zx2, f_gaussian func, double z);
    void TLM_scattering();
    void TLM_connection(double Zz1, double Zz2, double Zy1, double Zy2, double Zx1, double Zx2);
    void TLM_connection_sd(double Zz1, double Zz2, double Zy1, double Zy2, double Zx1, double Zx2);
    void TLM_connection_pml_df();
    void TLM_connection_2D_SCN(double Zz1, double Zz2, double Zy1, double Zy2, double Zx1, double Zx2);
    void TLM_connection_2D_SCN_wg(double Zz1, double Zz2, double Zy1, double Zy2, double Zx1, double Zx2);
    void TLM_connection_optimized(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 );
    void TLM_excitation_te10(int node , double v_inc);
    void TLM_excitation_single_node(int node , double v_inc, int axis=2, int srce_res=0);
    void TLM_excitation_plane_node(int node , double v_inc);
    void TLM_excitation_by_imposing_inc( int node, vector <double> scn_inc );
    double Ex_output_at_node(int node_id);    //selects which node to read an output from
    double Ey_output_at_node(int node_id);
    double Ez_output_at_node(int node_id);
    double Hx_output_at_node(int node_id);    //selects which node to read an output from
    double Hy_output_at_node(int node_id);
    double Hz_output_at_node(int node_id);
    double compute_field_average_in_z_plane(int z_plane);
    vector<double> compute_field_along_line( int x_plane_id,int y_plane, int z_plane_id);
    int get_true_centre_node();
    void print_Ey_output();
    void print_Ez_output();
    void print_SCN_neighbour(int node_id);
    void print_WG_nodes();
    void print_WG_plane_par(int z_plane, int par=1);
    void breakpoint_(int i);
    void set_special_nodes();
    void print_connect_bins();
    void create_connect_bins();
    void place_dipole_antenna(double dipole_length, int x=-1, int y=-1, int z=-1 );
    void insert_square_loop(double sq_length, int x=-1, int y=-1, int z=-1 );
    void insert_FSS_square(double thcknss , float pos);
    void insert_FSS_jerusalem_cross(double thcknss , float pos);
    void insert_PML(int num_PML, double loc, double Reflctn_f, int c_profile, int scn_typ);

    void insert_iris(int thickness , float z_plane);
    void insert_perfect_cube(int cube_l,int nx0,int ny0, int nz0);
    bool print_G_PML_parameters(int L_R,int PML_layer);
    void print_WG_node(int excited_node);
    void write_output_file (const string& nodefileEX,const string& nodefileEY, const string& nodefileEZ, const string& nodefileHX,const string& nodefileHY,const string& nodefileHZ);                                           //writes to node file
    void write_output_file_for_xy_plane (const string& nodefileEY) ;
    void write_output_file_for_yz_plane (const string& nodefileEY) ;
    void write_output_file_for_zx_plane (const string& nodefileEY) ;
    void write_output_file_for_line (const string& nodefileEY) ;
    void write_incident_file_for_xy_plane (const string& nodefileEY) ;

    void print_this_neighbour(int excited_node);
    void set_WG_neighbour();
    void analytical_fc();
    void emptyy();
    vector<double> compute_far_field(const vector<int> &plane_xy1);
    virtual ~mesh_handler_SCN();
private:
    int scn_type;
    double WG_width;
    double WG_height;
    double WG_length;
    double WG_dl;                // discretization length
    double WG_PML_length;
    int sim_type;
    int Nx;                     //number of nodes on the x , y and z axis
    int Ny;
    int Nz;
    int Nzz;
    int Nxx;
    int Nyy;
    int npml;
    int centre_node;              //node of excitation for short dipole - default set to -1 if waveguide
    vector < double > Is;
    vector <double > Vs;
    vector< double > WG_Ex;
    vector< double > WG_Ey;
    vector< double > WG_Ez;
    vector< double > WG_Hx;
    vector< double > WG_Hy;
    vector< double > WG_Hz;
    vector< double > average_Ey;
    vector <vector < double > > WG_Ey_plane_xy;
    vector <vector < double > > WG_Ey_plane_yz;
    vector <vector < double > > WG_Ey_plane_zx;
    vector <vector < double > > V_incident_xy;
    vector <vector < double > > Ey_field_along_line;
    vector <double> Ey_outputs;
    int Ntotal;
    int Ntotal_;
    vector<SCN_node> WG_nodes;   // represents the waveguide nodes; this is set in the constructor depending on the simulation parameters.
    vector< vector <unsigned int> > WG_node_neighbours;
    vector< vector <SCN_node*> > connect_bins;
    bool check_meshed;           //checks if the waveguide has been meshed / filled with cubes
    bool neighbours_set;
};

#endif
