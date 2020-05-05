

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

#include "HSCN_node.h"

using namespace std;

const   double pi = 3.14159265;
const   double  u0 = 12566370614e-16;
const   double e0 = 88541878176e-22;
const   double c = 299792458;
const   double Z0 = sqrt(u0/e0);//
const float ur(1);

HSCN_node::HSCN_node( int sim_identifier, float dl,   double dt, int id, int xx, int yy, int zz,bool special_node_flag,int medium_type,bool is_PEC, bool PML_check,int PML_thickness,   double Refn_fctr, int conduct_profile, int Lx, int Ly, int Lz)
{

    HSCN_node_sim_identifier = sim_identifier;

    HSCN_id=id;
    HSCN_dl =dl;
    HSCN_coord.push_back(xx);
    HSCN_coord.push_back(yy);
    HSCN_coord.push_back(zz);
    conduct_prof = conduct_profile;
    is_PML_node=PML_check;
    PEC = is_PEC;
    media_identifier = medium_type;
    special_node_set_flag =special_node_flag;

    //SETTING Time step dt
    HSCN_ports.assign(15,make_pair(0, 0) );

//SETTING THE CONDUCTIVITY
      double sigma = -(e0*0.5*c *log(Refn_fctr))/(1e-3*HSCN_dl*pow(PML_thickness,conduct_profile+1));
      double sigma_ex(0),sigma_ey(0),sigma_ez(0);
    if (Lx != -1) sigma_ex = sigma*(pow((Lx+1),conduct_prof+1)- pow(Lx,conduct_prof+1));
    if (Ly != -1) sigma_ey = sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1));
    if (Lz != -1) sigma_ez = sigma*(pow((Lz+1),conduct_prof+1)- pow(Lz,conduct_prof+1));

    cndtvty_x = sigma_ex*dt/e0;
    cndtvty_y = sigma_ey*dt/e0;
    cndtvty_z = sigma_ez*dt/e0;

// DIGITAL FILTER COMPONENTS

    V_j_tot_n.assign(3,0);
    V_IZ_tot_n.assign(6,0);

    if(PML_check==true)
    {

        V_inc_n_1.assign(15,0);
        V_inc_n_2.assign(15,0);

        //V_j_tot_n.assign(3,0);
        V_j_tot_n_1.assign(3,0);             //V_j.first = shunt voltage , and V_j.second = stored in previous time step
        V_j_tot_n_2.assign(3,0);

        V_IZ_tot_n.assign(6,0);
        V_IZ_tot_n_1.assign(6,0);             //V_j.first = shunt voltage , and V_j.second = stored in previous time step
        V_IZ_tot_n_2.assign(6,0);

    }

}

HSCN_node::HSCN_node (const HSCN_node& copy_HSCN_node)
{
        HSCN_node_sim_identifier = copy_HSCN_node.HSCN_node_sim_identifier;
        media_identifier = copy_HSCN_node.media_identifier;
        is_PML_node = copy_HSCN_node.is_PML_node;
        conduct_prof = copy_HSCN_node.conduct_prof;
        PEC = copy_HSCN_node.PEC;
        HSCN_dl = copy_HSCN_node.HSCN_dl;
        HSCN_id = copy_HSCN_node.HSCN_id;
        HSCN_ports = copy_HSCN_node.HSCN_ports;     // HSCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.

        cndtvty_x = copy_HSCN_node.cndtvty_x;
        cndtvty_y = copy_HSCN_node.cndtvty_y;
        cndtvty_z = copy_HSCN_node.cndtvty_z;

        V_inc_n_1 = copy_HSCN_node.V_inc_n_1;
        V_inc_n_2 = copy_HSCN_node.V_inc_n_2;

        V_j_tot_n = copy_HSCN_node.V_j_tot_n;
        V_j_tot_n_1 = copy_HSCN_node.V_j_tot_n_1;
        V_j_tot_n_2 = copy_HSCN_node.V_j_tot_n_2;

        V_IZ_tot_n = copy_HSCN_node.V_IZ_tot_n;
        V_IZ_tot_n_1 = copy_HSCN_node.V_IZ_tot_n_1;
        V_IZ_tot_n_2 = copy_HSCN_node.V_IZ_tot_n_2;

        HSCN_coord = copy_HSCN_node.HSCN_coord;
        special_node = copy_HSCN_node.special_node;
        special_node_set_flag = copy_HSCN_node.special_node_set_flag;
}

HSCN_node::HSCN_node()
{

}

HSCN_node &HSCN_node::operator=(const HSCN_node& copy_HSCN_node)
{
    if( this != &copy_HSCN_node)
        {
        HSCN_node_sim_identifier = copy_HSCN_node.HSCN_node_sim_identifier;
        media_identifier = copy_HSCN_node.media_identifier;
        is_PML_node = copy_HSCN_node.is_PML_node;
        conduct_prof = copy_HSCN_node.conduct_prof;
        PEC = copy_HSCN_node.PEC;
        HSCN_dl = copy_HSCN_node.HSCN_dl;
        HSCN_id = copy_HSCN_node.HSCN_id;
        HSCN_ports = copy_HSCN_node.HSCN_ports;     // HSCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.

        cndtvty_x = copy_HSCN_node.cndtvty_x;
        cndtvty_y = copy_HSCN_node.cndtvty_y;
        cndtvty_z = copy_HSCN_node.cndtvty_z;

        V_inc_n_1 = copy_HSCN_node.V_inc_n_1;
        V_inc_n_2 = copy_HSCN_node.V_inc_n_2;

        V_j_tot_n = copy_HSCN_node.V_j_tot_n;
        V_j_tot_n_1 = copy_HSCN_node.V_j_tot_n_1;
        V_j_tot_n_2 = copy_HSCN_node.V_j_tot_n_2;

        V_IZ_tot_n = copy_HSCN_node.V_IZ_tot_n;
        V_IZ_tot_n_1 = copy_HSCN_node.V_IZ_tot_n_1;
        V_IZ_tot_n_2 = copy_HSCN_node.V_IZ_tot_n_2;

        HSCN_coord = copy_HSCN_node.HSCN_coord;
        special_node = copy_HSCN_node.special_node;
        special_node_set_flag = copy_HSCN_node.special_node_set_flag;
        }
    return (*this);
}

 HSCN_node::~HSCN_node()
{
    //destructor
}

int HSCN_node::get_iD()
{
    return this->HSCN_id;
}

int HSCN_node::get_coord(int i)
{
    return(this->HSCN_coord[i]);
}
ostream& operator <<(ostream& out,HSCN_node& node)
{
    out<<" [ "<< node.HSCN_id<<" | "<< node.HSCN_coord[0] <<" "<< node.HSCN_coord[1] <<" "<<node.HSCN_coord[2]<< " ] "<<endl;
    return out;
}

void HSCN_node::excitation (int port_id,  double V_i)
{
    HSCN_ports[port_id].first=V_i;
}

void HSCN_node::set_neighbour(HSCN_node *node)
{

}

void HSCN_node::set_material_type(int media_type)
{
    media_identifier = media_type;
}

bool HSCN_node::check_special_node()
{
    return this->special_node_set_flag;
}

void HSCN_node::set_special_HSCN_node(bool flag)
{
    this->special_node_set_flag = flag;
}

void HSCN_node::print_neighbours()
{
}

bool HSCN_node::is_PML_HSCN()
{
    return (is_PML_node);
}


void HSCN_node::set_PEC()
{
    if(!is_PML_node) this->PEC = true;
}

bool HSCN_node::is_PEC()
{
    return (PEC);
}

   double HSCN_node::get_reflected(int port_id)
{
    return  HSCN_ports[port_id-1].second;
}

 void HSCN_node:: set_reflected( int port_id,  double Vr)
{
    HSCN_ports[port_id-1].second = Vr;
}


   double HSCN_node:: get_incidence( int port_id)
{
    return  HSCN_ports[port_id-1].first;
}

   double HSCN_node::get_Vn(int i)
 {
     return  V_j_tot_n[i-1];
 }

 void HSCN_node:: set_incidence(int port_id,  double Vi)
{
    HSCN_ports[port_id-1].first = Vi;
}

   double HSCN_node:: get_cndtvty_x()
{
    return  this->cndtvty_x;
}

  double HSCN_node:: get_cndtvty_y()
{
    return  this->cndtvty_y;
}

  double HSCN_node:: get_cndtvty_z()
{
    return  this->cndtvty_z;
}

int HSCN_node::get_medium_type()
{
    return (media_identifier);
}


void HSCN_node::set_all_to_one()
{
    //HSCN_ports.assign(24,make_pair(1,0));

    set_incidence(1,1);
    set_incidence(2,2);
    set_incidence(3,3);
    set_incidence(4,4);
    set_incidence(5,5);
    set_incidence(6,6);
    set_incidence(7,7);
    set_incidence(8,8);
    set_incidence(9,9);
    set_incidence(10,10);
    set_incidence(11,11);
    set_incidence(12,12);
    set_incidence(13,0);
    set_incidence(14,0);
    set_incidence(15,0);

}

void HSCN_node::print_node_voltages()
{

    cout<<endl;
    for(int i = 1; i<= 15 ; i++)
    {
        cout<<i<<" V_inc = "<<std::setprecision(numeric_limits<  double>::digits10+1)<<this->HSCN_ports[i-1].first << "    V_r = " << this->HSCN_ports[i-1].second <<endl;
    }

}

void HSCN_node::scattering_standard(vector<float> impedance_y, vector<float> conductance_g, vector<float> Rm) // needs looking into
{
    int ii;

      double Vx (0), Vy(0) , Vz(0);
      double Vxz(0), Vzx(0), Vyz(0), Vzy(0), Vyx(0), Vxy(0);
    double   y_xy = impedance_y[0] ,y_zy =impedance_y[1] , y_xz = impedance_y[2];
    float yox = impedance_y[3] , yoy = impedance_y[4] , yoz = impedance_y[5] ;
    float gex = conductance_g[0] , gey = conductance_g[1] , gez = conductance_g[2];

    double Rmx = Rm[0], Rmy = Rm[1] , Rmz = Rm[2];

    //this->set_all_to_one();

    //compute_standard_Vj( impedance_y, conductance_G , Vx, Vy ,Vz );

    //compute_V_xy_yx(y_xy, Rmz, Vyx, Vxy);

    double Yx_v = 1/ (y_xy + y_xz + (Z0*gex +yox)*0.5);
    double Yy_v = 1/ (y_xy + y_zy + (Z0*gey +yoy)*0.5);
    double Yz_v = 1/ (y_xz + y_zy + (Z0*gez +yoz)*0.5);

    Vx = Yx_v*( y_xy*(HSCN_ports[12-1].first + HSCN_ports[1-1].first)  + y_xz*(HSCN_ports[2-1].first + HSCN_ports[9-1].first)  + yox*HSCN_ports[13-1].first );

    Vy = Yy_v*( y_xy*(HSCN_ports[3-1].first +  HSCN_ports[11-1].first) + y_zy*(HSCN_ports[4-1].first + HSCN_ports[8-1].first)  + yoy*HSCN_ports[14-1].first );

    Vz = Yz_v*( y_zy*(HSCN_ports[7-1].first +  HSCN_ports[5-1].first)  + y_xz*(HSCN_ports[6-1].first + HSCN_ports[10-1].first) + yoz*HSCN_ports[15-1].first);


    float Zxy = Z0/y_xy;
    double Yz = 1/(2*Zxy + Rmz*0.5);
    double Iz = (HSCN_ports[1-1].first - HSCN_ports[12-1].first + HSCN_ports[11-1].first - HSCN_ports[3-1].first)*Yz;
    Vxy = Zxy*Iz;
    Vyx = Vxy;

    //compute_V_xz_zx(y_xz, Rmy, Vxz, Vzx);

    float Zxz = Z0/y_xz;
    double Yy = 1/(2*Zxz + Rmy*0.5);
    double Iy = (HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first)*Yy;
    Vxz = Zxz*Iy;
    Vzx = Vxz;

    //compute_V_yz_zy(y_zy, Rmx, Vyz, Vzy);

    double Zzy = Z0/y_zy;
    double Yx = 1/(2*Zzy + Rmx*0.5);
    double Ix = (HSCN_ports[7-1].first - HSCN_ports[5-1].first + HSCN_ports[4-1].first - HSCN_ports[8-1].first)*Yx;
    Vzy = Zzy*Ix;
    Vyz = Vzy;


    V_j_tot_n[0] = Vx; //total voltages
    V_j_tot_n[1] = Vy;
    V_j_tot_n[2] = Vz;
/*
    V_IZ_tot_n[0] = Vxy;
    V_IZ_tot_n[1] = Vyx;
    V_IZ_tot_n[2] = Vxz;
    V_IZ_tot_n[3] = Vzx;
    V_IZ_tot_n[4] = Vzy;
    V_IZ_tot_n[5] = Vyz;
*/
      double Vr1 = Vx - Vyx - get_incidence(12);
      double Vr2 = Vx + Vzx - get_incidence(9);
      double Vr3 = Vy + Vxy - get_incidence(11);
      double Vr4 = Vy - Vzy - get_incidence(8);
      double Vr5 = Vz + Vyz - get_incidence(7);
      double Vr6 = Vz - Vxz - get_incidence(10);
      double Vr7 = Vz - Vyz - get_incidence(5);
      double Vr8 = Vy + Vzy - get_incidence(4);
      double Vr9 = Vx - Vzx - get_incidence(2);
      double Vr10 = Vz + Vxz - get_incidence(6);
      double Vr11 = Vy - Vxy - get_incidence(3);
      double Vr12 = Vx + Vyx - get_incidence(1);
      double Vr13 = Vx - get_incidence(13);
      double Vr14 = Vy - get_incidence(14);
      double Vr15 = Vz - get_incidence(15);

    set_reflected(1, Vr1 );
    set_reflected(2, Vr2 );
    set_reflected(3, Vr3 );
    set_reflected(4, Vr4 );
    set_reflected(5, Vr5 );
    set_reflected(6, Vr6 );
    set_reflected(7, Vr7 );
    set_reflected(8, Vr8 );
    set_reflected(9, Vr9 );
    set_reflected(10,Vr10 );
    set_reflected(11, Vr11 );
    set_reflected(12, Vr12 );
    set_reflected(13, Vr13 );
    set_reflected(14, Vr14 );
    set_reflected(15, Vr15 );

//}  //STUB CONNECTION
    this->set_incidence(13,Vr13);            // capacitive stubs become  open circuits
    this->set_incidence(14,Vr14);
    this->set_incidence(15,Vr15);

}


void HSCN_node::scattering(vector<float> impedance_y, vector<float> conductance_G, vector<float> Rm) // needs looking into
{
    /*if( this->is_PML_HSCN()== false )
    {
        //cout<<"here"<<endl;
        this->scattering_standard(impedance_y, conductance_G, Rm);
    }

    else
    {
*/
    int ii;
    //if(HSCN_id == 1){
      double Vx (0), Vy(0) , Vz(0);
      double Vxz(0), Vzx(0), Vyz(0), Vzy(0), Vyx(0), Vxy(0);
      double   y_xy = impedance_y[0] ,y_zy =impedance_y[1] , y_xz = impedance_y[2];
      double Rmx = Rm[0], Rmy = Rm[1] , Rmz = Rm[2];

    //this->set_all_to_one();

    compute_complex_Vj( impedance_y, conductance_G , Vx, Vy ,Vz );

    compute_V_xy_yx(y_xy, Rmz, Vyx, Vxy);

    compute_V_xz_zx(y_xz, Rmy, Vxz, Vzx);
    compute_V_yz_zy(y_zy, Rmx, Vyz, Vzy);

    V_j_tot_n[0] = Vx;
    V_j_tot_n[1] = Vy;
    V_j_tot_n[2] = Vz;

    V_IZ_tot_n[0] = Vxy;
    V_IZ_tot_n[1] = Vyx;
    V_IZ_tot_n[2] = Vxz;
    V_IZ_tot_n[3] = Vzx;
    V_IZ_tot_n[4] = Vzy;
    V_IZ_tot_n[5] = Vyz;

      double Vr1 = Vx - Vyx - get_incidence(12);
      double Vr2 = Vx + Vzx - get_incidence(9);
      double Vr3 = Vy + Vxy - get_incidence(11);
      double Vr4 = Vy - Vzy - get_incidence(8);
      double Vr5 = Vz + Vyz - get_incidence(7);
      double Vr6 = Vz - Vxz - get_incidence(10);
      double Vr7 = Vz - Vyz - get_incidence(5);
      double Vr8 = Vy + Vzy - get_incidence(4);
      double Vr9 = Vx - Vzx - get_incidence(2);
      double Vr10 = Vz + Vxz - get_incidence(6);
      double Vr11 = Vy - Vxy - get_incidence(3);
      double Vr12 = Vx + Vyx - get_incidence(1);
      double Vr13 = Vx - get_incidence(13);
      double Vr14 = Vy - get_incidence(14);
      double Vr15 = Vz - get_incidence(15);


       double Iz = (HSCN_ports[1-1].first - HSCN_ports[12-1].first + HSCN_ports[11-1].first - HSCN_ports[3-1].first)/(2*Z0);
       double Ix = (HSCN_ports[7-1].first - HSCN_ports[5-1].first + HSCN_ports[4-1].first - HSCN_ports[8-1].first)/(2*Z0);
       double Iy = (HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first)/(2*Z0);

    set_reflected(1, exp(-cndtvty_y)*Vr1 );
    set_reflected(2, exp(-cndtvty_z)*Vr2 );
    set_reflected(3, exp(-cndtvty_x)*Vr3 );
    set_reflected(4, exp(-cndtvty_z)*Vr4 );
    set_reflected(5, exp(-cndtvty_y)*Vr5 );
    set_reflected(6, exp(-cndtvty_x)*Vr6 );
    set_reflected(7, exp(-cndtvty_y)*Vr7 );
    set_reflected(8, exp(-cndtvty_z)*Vr8 );
    set_reflected(9, exp(-cndtvty_z)*Vr9 );
    set_reflected(10, exp(-cndtvty_x)*Vr10 );
    set_reflected(11, exp(-cndtvty_x)*Vr11 );
    set_reflected(12, exp(-cndtvty_y)*Vr12 );
    set_reflected(13, Vr13 );
    set_reflected(14, Vr14 );
    set_reflected(15, Vr15 );

    update_voltages();
//}  //STUB CONNECTION
    this->set_incidence(13,Vr13);            // capacitive stubs become  open circuits
    this->set_incidence(14,Vr14);
    this->set_incidence(15,Vr15);


}

void HSCN_node::compute_Vj( vector<float> impedance_y, vector <float> conductance_g ,   double & Vx,   double &Vy ,
                             double &Vz,   double & Vx_n_1,   double &Vy_n_1 ,   double &Vz_n_1,   double & Vx_n_2,   double &Vy_n_2 ,   double &Vz_n_2)
{

    float yxy = impedance_y[0] , yyz = impedance_y[1] , yzx =impedance_y[2] ;
    float yox = impedance_y[3] , yoy = impedance_y[4] , yoz = impedance_y[5] ;

    float gex = conductance_g[0] , gey = conductance_g[1] , gez = conductance_g[2];

      double Yx = 1/ (yxy + yzx + (Z0*gex +yox)*0.5);
      double Yy = 1/ (yxy + yyz + (Z0*gey +yoy)*0.5);
      double Yz = 1/ (yzx + yyz + (Z0*gez +yoz)*0.5);


    int ii=0;

    Vx = Yx*( yxy*(HSCN_ports[12-1].first + HSCN_ports[1-1].first)  + yzx*(HSCN_ports[2-1].first + HSCN_ports[9-1].first)  + yox*HSCN_ports[13-1].first );

    Vy = Yy*( yxy*(HSCN_ports[3-1].first +  HSCN_ports[11-1].first) + yyz*(HSCN_ports[4-1].first + HSCN_ports[8-1].first)  + yoy*HSCN_ports[14-1].first );

    Vz = Yz*( yyz*(HSCN_ports[7-1].first +  HSCN_ports[5-1].first)  + yzx*(HSCN_ports[6-1].first + HSCN_ports[10-1].first) + yoz*HSCN_ports[15-1].first);



    Vx_n_1 = Yx*( yxy*(V_inc_n_1[12-1] + V_inc_n_1[1-1])  + yzx*(V_inc_n_1[2-1] + V_inc_n_1[9-1])  + yox*V_inc_n_1[13-1]);

    Vy_n_1 = Yy*( yxy*(V_inc_n_1[3-1]  + V_inc_n_1[11-1]) + yyz*(V_inc_n_1[8-1] + V_inc_n_1[4-1])  + yoy*V_inc_n_1[14-1]);

    Vz_n_1 = Yz*( yyz*(V_inc_n_1[7-1]  + V_inc_n_1[5-1] ) + yzx*(V_inc_n_1[6-1] + V_inc_n_1[10-1]) + yoz*V_inc_n_1[15-1]);



    Vx_n_2 =  Yx*( yxy*(V_inc_n_2[12-1] + V_inc_n_2[1-1])  + yzx*(V_inc_n_2[2-1] + V_inc_n_2[9-1])  + yox*V_inc_n_2[13-1]);

    Vy_n_2 =  Yy*( yxy*(V_inc_n_2[3-1]  + V_inc_n_2[11-1]) + yyz*(V_inc_n_2[8-1] + V_inc_n_2[4-1])  + yoy*V_inc_n_2[14-1]);

    Vz_n_2 =  Yz*( yyz*(V_inc_n_2[7-1]  + V_inc_n_2[5-1] ) + yzx*(V_inc_n_2[6-1] + V_inc_n_2[10-1]) + yoz*V_inc_n_2[15-1]);

/*
    cout<<endl;
    cout<<" Yx "<< Yx << " Yy" << Yy << " Yz" << Yz<<endl;
    cout<<" V_n "<< Vx << " Vy" << Vy << " Vz" << Vz<<endl;
    cout<<" V_n_1 "<< Vx_n_1 << " Vy" << Vy_n_1 << " Vz" << Vz_n_1<<endl;
    cout<<" V_n_2 "<< Vx_n_2 << " Vy" << Vy_n_2 << " Vz" << Vz_n_2<<endl;
    cout<<endl;

    cin>>ii;
*/
}


void HSCN_node::compute_standard_Vj( vector<float> impedance_y, vector <float> conductance_g ,   double & Vx,   double &Vy ,
                             double &Vz)
{

    float yxy = impedance_y[0] , yyz = impedance_y[1] , yzx =impedance_y[2] ;
    float yox = impedance_y[3] , yoy = impedance_y[4] , yoz = impedance_y[5] ;

    float gex = conductance_g[0] , gey = conductance_g[1] , gez = conductance_g[2];

      double Yx = 1/ (yxy + yzx + (Z0*gex +yox)*0.5);
      double Yy = 1/ (yxy + yyz + (Z0*gey +yoy)*0.5);
      double Yz = 1/ (yzx + yyz + (Z0*gez +yoz)*0.5);

    int ii=0;

    Vx = Yx*( yxy*(HSCN_ports[12-1].first + HSCN_ports[1-1].first)  + yzx*(HSCN_ports[2-1].first + HSCN_ports[9-1].first)  + yox*HSCN_ports[13-1].first );

    Vy = Yy*( yxy*(HSCN_ports[3-1].first +  HSCN_ports[11-1].first) + yyz*(HSCN_ports[4-1].first + HSCN_ports[8-1].first)  + yoy*HSCN_ports[14-1].first );

    Vz = Yz*( yyz*(HSCN_ports[7-1].first +  HSCN_ports[5-1].first)  + yzx*(HSCN_ports[6-1].first + HSCN_ports[10-1].first) + yoz*HSCN_ports[15-1].first);

}


void HSCN_node::compute_complex_Vj( vector<float> impedance_y, vector <float> conductance_g ,   double & Vx,   double &Vy ,   double &Vz)
{
   // VARIABLES DECLARATIONS + DEFINITIONS
      double Vx_n(0), Vy_n(0), Vz_n(0), Vx_n_1(0), Vy_n_1(0), Vz_n_1(0), Vx_n_2(0), Vy_n_2(0), Vz_n_2(0);

      double yxy = impedance_y[0] , yyz = impedance_y[1] , yzx = impedance_y[2] ;
      double yox = impedance_y[3] , yoy = impedance_y[4] , yoz = impedance_y[5] ;

      double gex = conductance_g[0] , gey = conductance_g[1] , gez = conductance_g[2];

      double Yx = (yxy + yzx + (Z0*gex +yox)*0.5) ;
      double Yy = (yxy + yyz + (Z0*gey +yoy)*0.5);
      double Yz = (yzx + yyz + (Z0*gez +yoz)*0.5);

    // COMPUTE AND RETURN BY REFERENCE THE STANDARD VOLTAGES
    int ii=0;

    compute_Vj( impedance_y, conductance_g,  Vx_n, Vy_n, Vz_n, Vx_n_1, Vy_n_1, Vz_n_1, Vx_n_2, Vy_n_2, Vz_n_2);

    //COMPUTING THE PML COEFFICIENTS FOR X
      double Gx = (Z0*gex +yox)*0.5;
      double bx = cndtvty_z*cndtvty_y*0.5;
      double ax = ( cndtvty_z + cndtvty_y + bx )*0.5;
      double cx = -( ( cndtvty_z + cndtvty_y - bx )*0.5);
      double dx = 2*Yx - Gx*bx ;
      double ex = ( 2*Yx - yxy*cndtvty_z - yzx*cndtvty_y - Gx*( cndtvty_z + cndtvty_y - bx ) )*0.5;
      double alpha_mhscn_x = 2/( 2*Yx + yxy*cndtvty_z  + yzx*cndtvty_y + Gx*( cndtvty_z + cndtvty_y + bx ) );

    //COMPUTING THE PML COEFFICIENTS FOR Y
      double Gy = (Z0*gey +yoy)*0.5;
      double by = cndtvty_z*cndtvty_x*0.5;
      double ay = ( cndtvty_z + cndtvty_x + by )*0.5;
      double cy = -( ( cndtvty_z + cndtvty_x - by )*0.5);
      double dy = 2*Yy - Gy*by ;
      double ey = ( 2*Yy - yxy*cndtvty_z - yyz*cndtvty_x - Gy*( cndtvty_z + cndtvty_x - by ) )*0.5;
      double alpha_mhscn_y = 2/( 2*Yy + yxy*cndtvty_z + yyz*cndtvty_x + Gy*( cndtvty_z + cndtvty_x + by ) );


    //COMPUTING THE PML COEFFICIENTS FOR Z
      double Gz = (Z0*gez +yoz)*0.5;
      double bz = cndtvty_y*cndtvty_x*0.5;
      double az = ( cndtvty_y + cndtvty_x + bz )*0.5;
      double cz = -( ( cndtvty_y + cndtvty_x - bz )*0.5);
      double dz = 2*Yz - Gz*bz ;
      double ez = ( 2*Yz - yyz*cndtvty_x - yzx*cndtvty_y - Gz*( cndtvty_x + cndtvty_y - bz ) )*0.5;
      double alpha_mhscn_z = 2/( 2*Yz + yyz*cndtvty_x + yzx*cndtvty_y + Gz*( cndtvty_x + cndtvty_y + bz ) );

    //cout<< " ??"<<endl;
    //cout<<Gz <<" "<<bz << " " << az << "" << cz <<" " <<dz << " " << ez << " " << alpha_mhscn_z<<endl;
    //cin>>ii;

    //COMPLEX VX
      double V_x_part1 = Yx *(Vx_n - 2*Vx_n_1 + Vx_n_2);
      double V_x_part2 = yxy*cndtvty_z*0.5* ( HSCN_ports[12-1].first + HSCN_ports[1-1].first -  V_inc_n_2[12-1] - V_inc_n_2[1-1]);
      double V_x_part3 = yzx*cndtvty_y*0.5* ( HSCN_ports[9-1].first  + HSCN_ports[2-1].first -  V_inc_n_2[9-1]  - V_inc_n_2[2-1]);
      double V_x_part4 = yox*( ax*HSCN_ports[13-1].first + bx*V_inc_n_1[13-1] + cx*V_inc_n_2[13-1]);
      double V_x_part5 = dx*V_j_tot_n_1[0] - ex*V_j_tot_n_2[0];

    Vx = alpha_mhscn_x*( V_x_part1 + V_x_part2 + V_x_part3 + V_x_part4 + V_x_part5 );
    //Vx = alpha_mhscn_x*Yx*Vx_n;

    //COMPLEX VY
      double V_y_part1 = Yy *(Vy_n);
      double V_y_part1_ =  - 2*Yy *Vy_n_1 + Yy *Vy_n_2;
      double V_y_part2 = yxy*cndtvty_z*0.5 *( HSCN_ports[11-1].first + HSCN_ports[3-1].first -  V_inc_n_2[11-1] - V_inc_n_2[3-1]);
      double V_y_part3 = yyz*cndtvty_x*0.5 *( HSCN_ports[4-1].first  + HSCN_ports[8-1].first -  V_inc_n_2[4-1]  - V_inc_n_2[8-1]);
      double V_y_part4 = yoy*( ay*HSCN_ports[14-1].first + by*V_inc_n_1[14-1] + cy*V_inc_n_2[14-1]);
      double V_y_part5 = dy*V_j_tot_n_1[1] - ey*V_j_tot_n_2[1];

    Vy = alpha_mhscn_y*( V_y_part1 + V_y_part1_ + V_y_part2 + V_y_part3 + V_y_part4 + V_y_part5 );
    //Vy = alpha_mhscn_y*Yy*Vy_n;
/*
    if( HSCN_id == 532)
    {
        cout<<" This should be zero "<< V_y_part2 << " " << V_y_part3 << " " << V_y_part4 << " " <<endl<<endl;
        cout<<" This should be zero too " << V_y_part1_ + V_y_part5 << " " <<endl;
        cout<< " " <<endl;
    }
*/
    //COMPLEX VZ
      double V_z_part1 = Yz *(Vz_n - 2*Vz_n_1 + Vz_n_2);
      double V_z_part2 = yzx*cndtvty_y*0.5 *( HSCN_ports[6-1].first + HSCN_ports[10-1].first -  V_inc_n_2[10-1] - V_inc_n_2[6-1]);
      double V_z_part3 = yyz*cndtvty_x*0.5* ( HSCN_ports[5-1].first  + HSCN_ports[7-1].first -  V_inc_n_2[5-1]  - V_inc_n_2[7-1]);
      double V_z_part4 = yoz*( az*HSCN_ports[15-1].first + bz*V_inc_n_1[15-1] + cz*V_inc_n_2[15-1]);
      double V_z_part5 = dz*V_j_tot_n_1[2] - ez*V_j_tot_n_2[2];

    Vz = alpha_mhscn_z*( V_z_part1 + V_z_part2 + V_z_part3 + V_z_part4 + V_z_part5 );
    //Vz = alpha_mhscn_z*Yz*Vz_n;

/*if( HSCN_id == 5320000000)
    {
        cout<<" difference one "<< -2*Yx*Vx_n_1 + dx*V_j_tot_n_1[0]<<endl<<endl;
        cout<<" difference two "<< Yx*Vx_n_2 - ex*V_j_tot_n_2[0]<<endl<<endl;
        cout<< " " <<endl;
    }
*/
    /*cout<<"zz"<<endl;
    cout<<Vz<<endl;
    cout<<alpha_mhscn_z<< " "<< V_z_part1 << " " << V_z_part2 << " " << V_z_part3 << " " << V_z_part4 << V_z_part5<<" "<<endl;
    cout<< " Y" << Yy << " "<< Yx << " " << Yz <<endl;
    cout<< " yoy " <<yoy << " " << yox<< " " <<yoz <<endl;
    cout<< " Vz " << Vz<<endl;
    //int ii=0;
    cin>>ii;*/

}

void HSCN_node::compute_V_xy_yx(float y_xy, float Rmz,   double & Vyx,   double &Vxy)
{

    float Zxy = Z0/y_xy;
      double Zz = 2*Zxy + Rmz*0.5;
      double Yz = 1/Zz;
      double axy = (4*Zxy + Rmz)*Yz;
      double bxy = cndtvty_x*cndtvty_y*0.25*Yz*Rmz;
      double cxy = (2*Zxy + Rmz)*(cndtvty_x + cndtvty_y)*0.25*Yz;
      double gxy = bxy - axy;
      double hxy = 0.5*( axy + bxy - 2*cxy);
      double b_MHSCN_z = 2/(axy + bxy +2*cxy);
      double Iz = (HSCN_ports[1-1].first - HSCN_ports[12-1].first + HSCN_ports[11-1].first - HSCN_ports[3-1].first)*Yz;
      double Iz_n_1 = ( V_inc_n_1[1-1] - V_inc_n_1[12-1] + V_inc_n_1[11-1] - V_inc_n_1[3-1])*Yz;
      double Iz_n_2 = ( V_inc_n_2[1-1] - V_inc_n_2[12-1] + V_inc_n_2[11-1] - V_inc_n_2[3-1])*Yz;

    Vxy = b_MHSCN_z* Zxy * (Iz );
      double Vxy_ = b_MHSCN_z* ( Zxy * ( - 2*Iz_n_1 + Iz_n_2 + cndtvty_y*0.5*(Iz - Iz_n_2)) - gxy*V_IZ_tot_n_1[0] - hxy*V_IZ_tot_n_2[0] );
    Vxy = Vxy_ + Vxy;

    Vyx = b_MHSCN_z* Zxy * (Iz);
      double Vyx_ =  b_MHSCN_z* ( Zxy * (- 2*Iz_n_1 + Iz_n_2 + cndtvty_x*0.5*(Iz - Iz_n_2)) - gxy*V_IZ_tot_n_1[1] - hxy*V_IZ_tot_n_2[1] );
    Vyx = Vyx_ + Vyx;

     if( HSCN_id == 53200000)
    {
        cout<<" difference three "<< 0.5*(HSCN_ports[1-1].first - HSCN_ports[12-1].first + HSCN_ports[11-1].first - HSCN_ports[3-1].first) - Vyx<<endl<<endl;
        cout<<" difference three "<< 0.5*(HSCN_ports[1-1].first - HSCN_ports[12-1].first + HSCN_ports[11-1].first - HSCN_ports[3-1].first) - Vxy<<endl<<endl;
        cout<< " " <<endl<<endl;
    }

    //Vxy = 0.5*(HSCN_ports[1-1].first - HSCN_ports[12-1].first + HSCN_ports[11-1].first - HSCN_ports[3-1].first);
   // Vyx = Vxy;

     /*int ii =0;

    cout<< " Zxy " <<Zxy << " " <<Zz << " " << axy << " " << bxy << " " << cxy << " " << gxy << " "<< hxy << " "<<b_MHSCN_z <<endl;
    cout<< " Iz " << Iz << " " << Iz_n_1 << " " << Iz_n_2;
    cout << " V " << Vxy << Vyx <<endl;

    cin>>ii;*/

}

void HSCN_node::compute_V_xz_zx(float y_xz, float Rmy,   double & Vxz,   double &Vzx)
{

    float Zxz = Z0/y_xz;
      double Zy = 2*Zxz + Rmy*0.5;

      double Yy = 1/Zy;
      double axz = (4*Zxz + Rmy)*Yy;
      double bxz = cndtvty_x*cndtvty_z*0.25*Yy*Rmy;
      double cxz = (2*Zxz + Rmy)*(cndtvty_x + cndtvty_z)*0.25*Yy;
      double gxz = bxz - axz;
      double hxz = .5*(axz + bxz -2*cxz);
      double b_MHSCN_y = 2/( axz + bxz + 2*cxz);
      double Iy = (HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first)*Yy;
      double Iy_n_1 = ( V_inc_n_1[6-1] - V_inc_n_1[10-1] + V_inc_n_1[9-1] - V_inc_n_1[2-1])*Yy;
      double Iy_n_2 =  ( V_inc_n_2[6-1] - V_inc_n_2[10-1] + V_inc_n_2[9-1] - V_inc_n_2[2-1])*Yy;

      double Vxz_ = b_MHSCN_y* ( Zxz * (- 2*Iy_n_1 + Iy_n_2 + cndtvty_z*0.5*(Iy - Iy_n_2)) - gxz*V_IZ_tot_n_1[2] - hxz*V_IZ_tot_n_2[2] );
    Vxz = b_MHSCN_y* Zxz * Iy;
    Vxz = Vxz_ + Vxz;

      double Vzx_ = b_MHSCN_y* ( Zxz * (- 2*Iy_n_1 + Iy_n_2 + cndtvty_x*0.5*(Iy - Iy_n_2)) - gxz*V_IZ_tot_n_1[3] - hxz*V_IZ_tot_n_2[3] );
    Vzx = b_MHSCN_y* Zxz * Iy;
    Vzx = Vzx_ + Vzx;

    if( HSCN_id == 532000000)
    {
        cout<<" difference two "<< 0.5*(HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first) - Vzx<<endl<<endl;
        cout<<" difference two "<< 0.5*(HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first) - Vxz<<endl<<endl;
        cout<< " " <<endl<<endl;
    }
    //Vxz = 0.5*(HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first);
   // Vzx = 0.5*(HSCN_ports[6-1].first - HSCN_ports[10-1].first + HSCN_ports[9-1].first - HSCN_ports[2-1].first);

   /*int ii =0;

    cout<< " Zxz " <<Zxz << " " <<Zy << " " << axz << " " << bxz << " " << cxz << " " << gxz << " "<< hxz << " "<<b_MHSCN_y <<endl;
    cout<< " Iy " << Iy << " " << Iy_n_1 << " " << Iy_n_2;
    cout << " V " << Vxz << " "<< Vzx <<endl;

    cin>>ii;
    */
}



void HSCN_node::compute_V_yz_zy(float y_zy, float Rmx,   double &Vyz ,   double &Vzy)
{
    double Zzy = Z0/y_zy;
    Rmx = Rmx/Zzy;

    double azy = (4 + Rmx );
    double bzy = cndtvty_z*cndtvty_y*0.25*Rmx;
    double czy = (2 + Rmx)*(cndtvty_z + cndtvty_y)*0.25;
    double gzy = (bzy - azy);
    double hzy = (azy + bzy -2*czy )*0.5;
    double b_MHSCN_x = 2/( azy + bzy + 2*czy);
    double Ix = (HSCN_ports[7-1].first - HSCN_ports[5-1].first + HSCN_ports[4-1].first - HSCN_ports[8-1].first);
    double Ix_n_1 = ( V_inc_n_1[7-1] - V_inc_n_1[5-1] + V_inc_n_1[4-1] - V_inc_n_1[8-1]);
    double Ix_n_2 = ( V_inc_n_2[7-1] - V_inc_n_2[5-1] + V_inc_n_2[4-1] - V_inc_n_2[8-1]);

    Vzy  = b_MHSCN_x*Ix;
    double Vzy_ =  b_MHSCN_x*( - 2*Ix_n_1 + Ix_n_2 + cndtvty_y*0.5*(Ix - Ix_n_2)) - b_MHSCN_x*gzy*V_IZ_tot_n_1[4] - b_MHSCN_x*hzy*V_IZ_tot_n_2[4] ;
    Vzy =  Vzy+Vzy_ ; // adding a small number to a big number eliminates rounding errors.

    Vyz  = b_MHSCN_x*Ix;
    double  Vyz_ =  b_MHSCN_x*( - 2*Ix_n_1 + Ix_n_2 + cndtvty_z*0.5*(Ix - Ix_n_2)) - b_MHSCN_x*gzy*V_IZ_tot_n_1[5] - b_MHSCN_x*hzy*V_IZ_tot_n_2[5] ;
    Vyz =  Vyz_+Vyz ;//

    //Vzy = 0.5*(HSCN_ports[7-1].first - HSCN_ports[5-1].first + HSCN_ports[4-1].first - HSCN_ports[8-1].first);
    //Vyz = 0.5*(HSCN_ports[7-1].first - HSCN_ports[5-1].first + HSCN_ports[4-1].first - HSCN_ports[8-1].first);
/*
     if( HSCN_id == 532)
    {

        cout<<" difference two- which should be zero  "<<  Vzy_ <<endl<<endl;
        cout<< " interesting n-1 "<< b_MHSCN_x*( - 2*Ix_n_1)/(b_MHSCN_x*gzy*V_IZ_tot_n_1[4]) <<endl<<endl;
        cout<< " interesting n-2 "<<b_MHSCN_x*Ix_n_2 / (b_MHSCN_x*hzy*V_IZ_tot_n_2[4] )<<endl<<endl<<endl;

    }
*/

    /*int ii =0;

    cout<< " Zzy " <<Zzy << " " <<Zx << " " << azy << " " << bzy << " " << czy << " " << gzy << " "<< hzy << " "<<b_MHSCN_x <<endl;
    cout<< " Ix " << Ix << " " << Ix_n_1 << " " << Ix_n_2;
    cout << " V " << Vzy << Vyz <<endl;

    cin>>ii;*/

}

void HSCN_node::update_voltages( )
{
    V_j_tot_n_2 = V_j_tot_n_1;
    V_j_tot_n_1 = V_j_tot_n;

    V_IZ_tot_n_2 = V_IZ_tot_n_1;
    V_IZ_tot_n_1 = V_IZ_tot_n;

    V_inc_n_2 = V_inc_n_1;

    this->V_inc_n_1[0] = this->get_incidence(1);
    this->V_inc_n_1[1] = this->get_incidence(2);
    this->V_inc_n_1[2] = this->get_incidence(3);
    this->V_inc_n_1[3] = this->get_incidence(4);
    this->V_inc_n_1[4] = this->get_incidence(5);
    this->V_inc_n_1[5] = this->get_incidence(6);
    this->V_inc_n_1[6] = this->get_incidence(7);
    this->V_inc_n_1[7] = this->get_incidence(8);
    this->V_inc_n_1[8] = this->get_incidence(9);
    this->V_inc_n_1[9] = this->get_incidence(10);
    this->V_inc_n_1[10] = this->get_incidence(11);
    this->V_inc_n_1[11] = this->get_incidence(12);
    this->V_inc_n_1[12] = this->get_incidence(13);
    this->V_inc_n_1[13] = this->get_incidence(14);
    this->V_inc_n_1[14] = this->get_incidence(15);

}



