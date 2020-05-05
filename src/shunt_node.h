#ifndef SHUNT_NODE_H
#define SHUNT_NODE_H

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

//#include "SCN_node.h"

using namespace std;

class shunt_node
{
    public:
        //CONSTRUCTORS
        shunt_node (bool freespace_,float err,float sigma_e,int id_, double dt_,int xx, int yy, double dl_, bool isPML_=0, int PML_x=0, int PML_y=0, double sigmax=0,double sigmay=0,int cndct_prof=0);
        shunt_node(const shunt_node& copy_shunt_node);                                   // copy constructor
        shunt_node & operator=(const shunt_node& copy_shunt_node);                       // assignment operator
        //TLM PROCESS
        void shunt_scatter();                                                               // scattering
        void shunt_scatter_with_stubs ();
        void shunt_scatter_with_stubs2 ();
        void shunt_scatter_with_stubs3 ();
        void shunt_excitation (double V_i);                                   // sets port_id and incidence voltage

        //OPERATOR OVERLOAD <<
        friend ostream& operator <<(std::ostream& out,shunt_node& node);                 // output streaming
        //GETS
        double get_reflected(int port_id);
        double get_incidence( int port_id);
        void compute_source_voltages(float ez_inc,int source_res,float &Vz,float &Vz_prime);
        float get_cnductvty_x();
        float get_cnductvty_y();
        inline bool is_PEC(){return (this->PEC);}
        inline bool is_PML(){return (this->isPML);}
        inline void make_PEC_node(){this->PEC = true;}
        int get_iD();
        double get_dt();
        int get_coord(int i);
        double get_Ez();
        double get_Hx();
        double get_Hy();
        //SETS
        void  set_reflected( int port_id,double Vr);
        void  set_incidence(int port_id,double Vr);
        //MISC
        bool  is_PML_shunt();
        void print_node_voltages();
        void compute_total_i();
        //UPDATE
        void update_voltages_shunt_node();

        //DESTRUCTOR
        virtual ~shunt_node();

    private:

        vector< pair <double,double> > shunt_port;
        vector< int > shunt_coords;
        int id;
        double dt;
        double dl;
        double v_total;
        double i_pml_sum;
        bool PEC;
       // double i_sum;

       //STUBS
        bool freespace;
        float er_;
        float Ge;
        double Yst;

        //PML
        bool isPML;                                             //
        int PML_Lx;
        int PML_Ly;
        float cnductvty_x;                                        //conductivity in x
        float cnductvty_y;                                        //conductivity in y
        vector< pair < double,double> > v_i_ndt__1;             // stored voltage and current from previous time step
        vector < double > v_ndt__2;
        pair <double,double>  v_i_total_ndt__1;                 // stored total voltage and current from previous time step.
        double v_total_ndt__2;

};

#endif // SHUNT_NODE_H
