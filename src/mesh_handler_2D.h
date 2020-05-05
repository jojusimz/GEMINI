#ifndef mesh_handler_2D_H
#define mesh_handler_2D_H

#include "shunt_node.h"

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
#include <limits>   // double precision
#include "mesh_handler_SCN.h"
int return_coordinates_2D(double width, double height,int npmlx,int npmly, double dl, double x, double y, double z);

class naca_f
{
public:
    naca_f(float t_) : t(t_){}
    inline int operator () (double x, float dl) const
    {

       double f =  5*t*( 0.2969*sqrt(x)   -   0.1260*x  -  0.3516*pow(x,2)   +   0.2843*pow(x,3)  -  0.1015*pow(x,4)  );

       int ret = int( (1e3*f/dl) + 0.5);
      // cout<<" f here : "<<ret<<endl;
       return (ret);

    }

private:
    float t;


};
class mesh_handler_2D
{
    public:
        mesh_handler_2D(bool freespace,float er, float sigma_e,float h, float w,float dl_,int PMLx=0, int PMLy=0, double Refn_fctr=1, int cndct_prof=0);

        mesh_handler_2D(const mesh_handler_2D& copy_mesh_handler_2D);
        mesh_handler_2D(mesh_handler_2D& geom, int PML_x ,int PML_y);
        mesh_handler_2D & operator=(const mesh_handler_2D & Other);
        void Simulate_2D(int ttltimestep,vector<int> &bdry_cdn,int excit_iD,int excit_lngth,int obsrv_iD,f_gaussian func, int exct_typ=1, int d =0, int src_res = -1);
        void Geometry_excitation_Ez (float v_inc, int excit_iD,int source_res=-1);  // use a functor here to pass excitation function.
        double Geometry_observation_Ez (int obsrv_iD);
        double Geometry_observation_Hx (int obsrv_iD);
        double far_field_RCS(int bdry, double r);
        void insert_structure_Naca_x_axis(int start_node,int end_node,naca_f func);
        void insert_structure_dipole_y_axis(int mid_point,int length);
        void insert_iris_2D(int thickness , float pos);
        void write_output_file_2D (const string& nodefileEZ,const string& nodefileHX,const string& nodefileIs,const string& nodefileVs);                                         //writes to node file
        void write_output_file_for_xy_plane_shuntnode(const string& nodefileEZ) ;
        void print_nodes(int propty=1);
        void print_all_PEC();
        void make_PEC(int node_id);
        void TLM_scatter();
        void TLM_connect(vector<int> &bdry_cdn);
        virtual ~mesh_handler_2D();
    private:
        vector < shunt_node > prob_domain2D;
        vector < int > bdry_cdn;
        vector <vector < double > > Ez_plane_xy;
        vector<double> Ez_output;
        vector<double> Ez_far_field;
        vector<double> Hx_output;
        vector <double > Is;
        vector <double >Vs;
        float height;
        float width;
        float dl;
        int Ny;
        int Nx;
        int Nxy;
        int npmlx;
        int npmly;
        int Nxx;
        int Nyy;
        int Nxxyy;
        int nodetype;
        //Simulation related member variables

};

#endif // mesh_handler_2D_H
