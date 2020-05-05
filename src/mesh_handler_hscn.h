

#ifndef mesh_handler_HSCN_H
#define mesh_handler_HSCN_H

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
#include <limits>   //  double precisio
#include "HSCN_node.h"
#include "mesh_handler_SCN.h"

using namespace std;

/*
int return_coordinates( double width,  double height,  double length,int npml,  double dl,  double x,  double y,  double z);
//vector<vector< double> > parse2DCsvFile(int &simtype, string inputFileName, vector< string >& e_filemames, vector< string >& h_filenames);

class f_gaussian
{
public:
    f_gaussian( double f0_, double bw,  double dt_) : f0(f0_),bandwidth(bw),dt(dt_){  }
    inline  double operator () (int ti) const { return ( (dt*ti-(dt))*bandwidth* exp( -pow( bandwidth*(dt*ti-(dt*10) ),2))*sin(2*3.14159265*f0*dt*ti)  ); }

private:
     double f0;
     double bandwidth;
     double dt;
};
*/

class mesh_handler_HSCN
{
public:
    mesh_handler_HSCN(int simtype, double width,  double height,  double length,  double dl,float ery, float urx,float sigma_ey, float sigma_mx, int n_PMl=0,int conduct_prof=0, double Refn_factor=1);
    mesh_handler_HSCN( double width,  double height,  double length,int simtype,  double dl,float ery, float urx,float sigma_ey, float sigma_mx,int n_PMl=0,int conduct_prof=0, double Refn_factor=1);
    mesh_handler_HSCN( double width,  double height,int simtype,  double dl,float ery, float urx,float sigma_ey, float sigma_mx, int n_PMl=0,int conduct_prof=0, double Refn_factor=1);
    mesh_handler_HSCN(const mesh_handler_HSCN& copy_mesh_handler_HSCN);
    mesh_handler_HSCN & operator=(const mesh_handler_HSCN & Other);
    void Re_mesh_mesh_handler_HSCN( double width,  double height,  double length,  double dl);              // creates a new rectangular waveguide and re-meshes using the input parameters.
    void set_neighbours();           // finds the neighbours for each node and points them to each other.
    void print_all_neighbours();
    void TLM_simulation1(int dt_total,int output_node, int excited_node,  double Zz1,  double Zz2,  double Zx1,  double Zx2,  double Zy1,  double Zy2, f_gaussian func ,  double excitation_type=-1) ;     //simulates wave propagation using structured TLM
    void TLM_simulation_2D_HSCN(int dt_total, int output_node,int excited_node, double Zz1,  double Zz2,  double Zy1,  double Zy2, double Zx1,  double Zx2, f_gaussian func,  double z);
    void TLM_scattering();
    void TLM_connection( double Zz1,  double Zz2,  double Zy1,  double Zy2,  double Zx1,  double Zx2);
    void TLM_connection_sd( double Zz1,  double Zz2,  double Zy1,  double Zy2,  double Zx1,  double Zx2);
    void TLM_connection_pml_df();
    void TLM_connection_2D_HSCN( double Zz1,  double Zz2,  double Zy1,  double Zy2,  double Zx1,  double Zx2);
    void TLM_connection_2D_HSCN_mesh_handler_HSCN( double Zz1,  double Zz2,  double Zy1,  double Zy2,  double Zx1,  double Zx2);
    void create_connect_bins();
    void TLM_scatter_optimized();
    void print_connect_bins();
    void TLM_connection_optimized(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 );
    void TLM_excitation_te10(int node ,  double v_inc);
    void TLM_excitation_single_node(int node ,  double v_inc, int axis=2);
    void TLM_excitation_plane_node(int node ,  double v_inc, int exctn_plane=1);
    void E_output_at_node(int output_node);
    void H_output_at_node(int output_node);
     double Ex_output_at_node(int node_id,  double Vx);    //selects which node to read an output from
     double Ey_output_at_node(int node_id,  double Vy);
     double Ez_output_at_node(int node_id,  double Vz);
     double Hx_output_at_node(int node_id,  double Ix=0);    //selects which node to read an output from
     double Hy_output_at_node(int node_id,  double Iy=0);
     double Hz_output_at_node(int node_id,  double Iz=0);
    int get_true_centre_node() {return (centre_node);}
    int get_node_coord(int node_iD, int coord_x_y_z);
    int get_Nxx() {return (Nxx); }
    int get_Nx() {return (Nx);}
    int get_Nyy() {return (Nyy); }
    int get_Ny() {return (Ny);}
    int get_Nzz() {return (Nzz); }
    int get_Nz() {return (Nz);}
    void print_Ey_output();
    void print_Ez_output();
    void print_HSCN_neighbour(int node_id);
    void print_mesh_handler_HSCN_nodes();
    void print_mesh_handler_HSCN_plane_par(int z_plane, int par=1);
    void breakpoint_(int i);
    void set_special_nodes();
    void place_dipole_antenna( double dipole_length, int x=-1, int y=-1, int z=-1 );
    void insert_iris(int thickness , float z_plane);
    void insert_perfect_cube(int cube_l,int nx0,int ny0, int nz0);
    void insert_dielectric_sheet(double sheet_thickness, double loc);
    bool print_G_PML_parameters(int L_R,int PML_layer);
    void print_mesh_handler_HSCN_node(int excited_node);
    void write_output_file (const string& nodefileEX,const string& nodefileEY, const string& nodefileEZ, const string& nodefileHX,const string& nodefileHY,const string& nodefileHZ);                                           //writes to node file
    void write_output_file_for_xy_plane (const string& nodefileEY) ;
    void write_output_file_for_yz_plane (const string& nodefileEY) ;
    void write_output_file_for_zx_plane (const string& nodefileEY) ;
    void print_this_neighbour(int excited_node);
    void set_mesh_handler_HSCN_neighbour();
    void analytical_fc();
    vector< double> compute_far_field(const vector<int> &plane_xy1);
    virtual ~mesh_handler_HSCN();
private:
     double mesh_handler_HSCN_width;
     double mesh_handler_HSCN_height;
     double mesh_handler_HSCN_length;
     double mesh_handler_HSCN_dl;                // discretization length
     double mesh_handler_HSCN_PML_length;
     double timestep;
    vector< float > line_Y;
    vector< float > line_G;
    vector< float > line_R;
    vector< float > line_Y2;
    vector< float > line_G2;
    vector< float > line_R2;
    int media_identifier;
    int sim_type;
    int Nx;                     //number of nodes on the x , y and z axis
    int Ny;
    int Nz;
    int Nzz;
    int Nxx;
    int Nyy;
    int Ntotal;
    int Ntotal_;
    int npml;
    int centre_node;              //node of excitation for short dipole - default set to -1 if waveguide
    vector<  double > mesh_handler_HSCN_Ex;
    vector<  double > mesh_handler_HSCN_Ey;
    vector<  double > mesh_handler_HSCN_Ez;
    vector<  double > mesh_handler_HSCN_Hx;
    vector<  double > mesh_handler_HSCN_Hy;
    vector<  double > mesh_handler_HSCN_Hz;
    vector <vector <  double > > mesh_handler_HSCN_Ey_plane_xy;
    vector <vector <  double > > mesh_handler_HSCN_Ey_plane_yz;
    vector <vector <  double > > mesh_handler_HSCN_Ey_plane_zx;
    vector<HSCN_node> mesh_handler_HSCN_nodes;   // represents the waveguide nodes; this is set in the constructor depending on the simulation parameters.
    vector< vector <unsigned int> > mesh_handler_HSCN_node_neighbours;
    vector< vector <HSCN_node*> > connect_bins;
    vector< vector <unsigned int> > vctr_medium_iD;
    vector< vector <unsigned int> > vctr_PML_medium_iD;
    bool check_meshed;           //checks if the waveguide has been meshed / filled with cubes
    bool neighbours_set;
};

#endif
