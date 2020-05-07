
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

using namespace std;

const double pi = 3.14159265;
const double  u0 = 12566370614e-16;
const double e0 = 88541878176e-22;
const double c = 299792458.;
const double Z0 = sqrt(u0/e0);//
const float ur(1);

SCN_node::SCN_node(int scn_typ, double dl, double tfactor, float er,float sigma_y, int id, int xx, int yy, int zz,bool special_node_flag,
                      bool is_PEC, bool PML_check,int PML_thickness,double Refn_fctr,int conduct_profile, int Lx, int Ly, int Lz)
{

    soft_source_node = false;
    SCN_id=id;
    SCN_dl =dl;
    SCN_coord.push_back(xx);
    SCN_coord.push_back(yy);
    SCN_coord.push_back(zz);
    conduct_prof = conduct_profile;
    is_PML_node=PML_check;
    PEC = is_PEC;
    PML_Lx = Lx;
    PML_Ly = Ly;
    PML_Lz = Lz;
    special_node_set_flag =special_node_flag; //default setting
    scn_type = scn_typ;

    //cout<<" instantiation:"<<"ID "<< id << " "<<"PEC" << PEC <<" PML" << PML_check<<endl;

    //SETTING Time step dt
    if( scn_type >=2) tfactor =1;
    double dt = (dl*1e-3/(2*c)) *tfactor;

    if( PML_check ==true ){

    SCN_ports.assign(24,make_pair(0, 0) );

    double sigma = -(e0*0.5*c *log(Refn_fctr))/(1e-3*SCN_dl*pow(PML_thickness,conduct_profile+1));

    double sigma_ex(0),sigma_ey(0),sigma_ez(0);
    if (Lx != -1) sigma_ex = sigma*(pow((Lx+1),conduct_prof+1)- pow(Lx,conduct_prof+1));
    if (Ly != -1) sigma_ey = sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1));
    if (Lz != -1) sigma_ez = sigma*(pow((Lz+1),conduct_prof+1)- pow(Lz,conduct_prof+1));

    alpha_x = 1; alpha_y =1 ;alpha_z=1;
    if (Lx != -1)
    {
        if (Lx<=4)alpha_x = 0.6*(Lx+1)+1;//5; //sigma*(pow((Lx+1),conduct_prof+1)- pow(Lx,conduct_prof+1));
        else alpha_x = 0.6*(5)+1;
         //alpha_x =1;
        //cout<<alpha_x<<endl;
    }

    if (Ly != -1)
    {
        if (Ly<=4)alpha_y = 0.6*(Ly+1)+1;//5; //sigma_ey = sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1));
        else alpha_y = 0.6*(5)+1;
         //alpha_y =1;
        //cout<<alpha_y<<endl;
    }

    if (Lz != -1)
    {
       if (Lz<=4)alpha_z = .6*(Lz+1)+1;//5; // sigma*(pow((Lz+1),conduct_prof+1)- pow(Lz,conduct_prof+1))
       else alpha_z = .6*(5)+1;
        //alpha_z =1;
        //cout<<"z "<<alpha_z<<endl;
    }

    //Parameters along the X axis
    sigma_e.push_back(sigma_ex);
    sigma_h.push_back(sigma_e[0]*u0/e0);
    G.push_back(sigma_e[0]*Z0*1e-3*dl);
    R.push_back(sigma_h[0]*dl*1e-3/Z0);

    //Parameters along the Y axis
    //sigma_e.push_back(sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1)));
    sigma_e.push_back(sigma_ey);
    sigma_h.push_back(sigma_e[1]*u0/e0);
    G.push_back(sigma_e[1]*Z0*1e-3*dl);
    R.push_back(sigma_h[1]*dl*1e-3/Z0);

    //Parameters along the Z axis
    //sigma_e.push_back(sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1)));
    sigma_e.push_back(sigma_ez);
    sigma_h.push_back(sigma_e[2]*u0/e0);
    G.push_back(sigma_e[2]*Z0*1e-3*dl);
    R.push_back(sigma_h[2]*dl*1e-3/Z0);

    if( scn_type == 2) // PML digital filter
    {
            G.assign(3,0);
            R.assign(3,0);
            cndtvty_x = sigma_e[0]*dt/e0;
            cndtvty_y = sigma_e[1]*dt/e0;
            cndtvty_z = sigma_e[2]*dt/e0;
          //  cndtvty_z = 1.6798;
            V_inc_n_1.assign(12,0);
            V_j_n_1.assign(3,0);             //V_j.first = shunt voltage , and V_j.second = stored in previous time step
            V_j_pml_n_1.assign(3,0);
            V_j_tot_n_1.assign(3,0);
            I_j_sh_pml_n_1.assign(3,0);
            I_j_sh_pml_a_n_1.assign(12,0);
            I_j_sh_tot_a_n_1.assign(12,0);
            I_k_n_1.assign(3,0);
            I_k_pml_n_1.assign(3,0);
            I_k_tot_n_1.assign(3,0);
            V_k_series_n_1.assign(3,0);
            V_kj_dummy_n_1.assign(6,0);
            V_kj_dummy_tot_n_1.assign(6,0);

    }

    }
    else
    {

    SCN_ports.assign(24,make_pair(0, 0) );
    sigma_e.assign(3,0);
    sigma_h.assign(3,0);
    G.assign(3,sigma_y*Z0*1e-3*dl);
    G[1]=0;
    //G[1]=1000*356*1e-3*dl;//1;

    R.assign(3,0);
    cndtvty_x = 0;
    cndtvty_y = 0;
    cndtvty_z = 0;
    alpha_x = 1;
    alpha_y = 1 ;
    alpha_z = 1;

    V_inc_n_1.assign(12,0);
    V_j_n_1.assign(3,0);             //V_j.first = shunt voltage , and V_j.second = stored in previous time step
    V_j_pml_n_1.assign(3,0);
    V_j_tot_n_1.assign(3,0);
    I_j_sh_pml_n_1.assign(3,0);
    I_j_sh_pml_a_n_1.assign(12,0);
    I_j_sh_tot_a_n_1.assign(12,0);
    I_k_n_1.assign(3,0);
    I_k_pml_n_1.assign(3,0);
    I_k_tot_n_1.assign(3,0);
    V_k_series_n_1.assign(3,0);
    V_kj_dummy_n_1.assign(6,0);
    V_kj_dummy_tot_n_1.assign(6,0);

    }

    //SETTING THE STUB ADMITTANCE AND IMPEDANCE
    double y_ =  (2*er*dl*1e-3)/(c*dt) - 4 ; // note: this is assuming the grid dimensions are all equal to dl; an assumption made throughout this code.
    Y.assign(3,y_);

    double z_ = (2*ur*dl*1e-3)/(c*dt)   - 4 ;
    Z.assign(3,z_);


}

SCN_node::SCN_node (const SCN_node& copy_SCN_node)
{
    SCN_id = copy_SCN_node.SCN_id;
    PEC = copy_SCN_node.PEC;
    SCN_dl = copy_SCN_node.SCN_dl;
    conduct_prof = copy_SCN_node.conduct_prof;
    sigma_e = copy_SCN_node.sigma_e;            // sigma_e[0]: sigma_ex ; sigma_e[1] : sigma_ey;
    sigma_h = copy_SCN_node.sigma_h;
    is_PML_node = copy_SCN_node.is_PML_node;
    G = copy_SCN_node.G;                    //we assume length in all dimension is equal dx=dy=dz
    R = copy_SCN_node.R;
    Y = copy_SCN_node.Y;
    Z = copy_SCN_node.Z;
    Refn_fctr = copy_SCN_node.Refn_fctr;
    PML_Lx = copy_SCN_node.PML_Lx;
    PML_Ly = copy_SCN_node.PML_Ly;
    PML_Lz = copy_SCN_node.PML_Lz;
    SCN_ports = copy_SCN_node.SCN_ports;     // SCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.

    soft_source_node = copy_SCN_node.soft_source_node;

    //SCN_neighbours = copy_SCN_node.SCN_neighbours;
    SCN_coord = copy_SCN_node.SCN_coord;
    special_node = copy_SCN_node.special_node;
    special_node_set_flag = copy_SCN_node.special_node_set_flag;

    cndtvty_x = copy_SCN_node.cndtvty_x;
    cndtvty_y = copy_SCN_node.cndtvty_y;
    cndtvty_z = copy_SCN_node.cndtvty_z;
    alpha_x = copy_SCN_node.alpha_x;
    alpha_y = copy_SCN_node.alpha_y;
    alpha_z = copy_SCN_node.alpha_z;

    V_inc_n_1 = copy_SCN_node.V_inc_n_1;
    V_j_n_1   = copy_SCN_node.V_j_n_1;
    V_j_pml_n_1 = copy_SCN_node.V_j_pml_n_1;
    V_j_tot_n_1 = copy_SCN_node.V_j_tot_n_1;
    I_j_sh_pml_n_1 = copy_SCN_node.I_j_sh_pml_n_1;
    I_j_sh_pml_a_n_1 = copy_SCN_node.I_j_sh_pml_a_n_1;
    I_j_sh_tot_a_n_1 = copy_SCN_node.I_j_sh_tot_a_n_1;
    I_k_n_1 = copy_SCN_node.I_k_n_1;
    I_k_pml_n_1 = copy_SCN_node.I_k_pml_n_1;
    I_k_tot_n_1 = copy_SCN_node.I_k_tot_n_1;
    V_k_series_n_1 = copy_SCN_node.V_k_series_n_1;
    V_kj_dummy_n_1 = copy_SCN_node.V_kj_dummy_n_1;
    V_kj_dummy_tot_n_1=copy_SCN_node.V_kj_dummy_tot_n_1;
}

SCN_node::SCN_node()
{
    SCN_id = 0;
    PEC = 1;
    SCN_dl = 0;
    conduct_prof = 0;
    sigma_e.assign(3,1);                       // sigma_e[0]: sigma_ex ; sigma_e[1] : sigma_ey;
    sigma_h.assign(3,1);
    is_PML_node = 0;
    G.assign(3,1);                              //we assume length in all dimension is equal dx=dy=dz
    R.assign(3,1);
    Y.assign(3,1);
    Z.assign(3,1);
    Refn_fctr = 0;
    PML_Lx = 0;
    PML_Ly = 0;
    PML_Lz = 0;
    SCN_ports.assign( 24,make_pair(1,1));               // SCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.
    //SCN_neighbours = copy_SCN_node.SCN_neighbours;
    SCN_coord.assign(3,0);
    special_node = 0;
    special_node_set_flag = 0;

}
SCN_node &SCN_node::operator=(const SCN_node& copy_SCN_node)
{
    if( this != &copy_SCN_node){
    SCN_id = copy_SCN_node.SCN_id;
    PEC = copy_SCN_node.PEC;
    SCN_dl = copy_SCN_node.SCN_dl;
    conduct_prof = copy_SCN_node.conduct_prof;
    sigma_e = copy_SCN_node.sigma_e;            // sigma_e[0]: sigma_ex ; sigma_e[1] : sigma_ey;
    sigma_h = copy_SCN_node.sigma_h;
    is_PML_node = copy_SCN_node.is_PML_node;
    G = copy_SCN_node.G;                    //we assume length in all dimension is equal dx=dy=dz
    R = copy_SCN_node.R;
    Y = copy_SCN_node.Y;
    Z = copy_SCN_node.Z;
    Refn_fctr = copy_SCN_node.Refn_fctr;
    PML_Lx = copy_SCN_node.PML_Lx;
    PML_Ly = copy_SCN_node.PML_Ly;
    PML_Lz = copy_SCN_node.PML_Lz;
    SCN_ports = copy_SCN_node.SCN_ports;     // SCN_ports where each port has an incidence and a reflected. the first is incidence and the second is reflected.

    soft_source_node = copy_SCN_node.soft_source_node;

    //SCN_neighbours = copy_SCN_node.SCN_neighbours;
    SCN_coord = copy_SCN_node.SCN_coord;
    special_node = copy_SCN_node.special_node;
    special_node_set_flag = copy_SCN_node.special_node_set_flag;

    cndtvty_x = copy_SCN_node.cndtvty_x;
    cndtvty_y = copy_SCN_node.cndtvty_y;
    cndtvty_z = copy_SCN_node.cndtvty_z;
    alpha_x = copy_SCN_node.alpha_x;
    alpha_y = copy_SCN_node.alpha_y;
    alpha_z = copy_SCN_node.alpha_z;
    V_inc_n_1 = copy_SCN_node.V_inc_n_1;
    V_j_n_1   = copy_SCN_node.V_j_n_1;
    V_j_pml_n_1 = copy_SCN_node.V_j_pml_n_1;
    V_j_tot_n_1 = copy_SCN_node.V_j_tot_n_1;
    I_j_sh_pml_n_1 = copy_SCN_node.I_j_sh_pml_n_1;
    I_j_sh_pml_a_n_1 = copy_SCN_node.I_j_sh_pml_a_n_1;
    I_j_sh_tot_a_n_1 = copy_SCN_node.I_j_sh_tot_a_n_1;
    I_k_n_1 = copy_SCN_node.I_k_n_1;
    I_k_pml_n_1 = copy_SCN_node.I_k_pml_n_1;
    I_k_tot_n_1 = copy_SCN_node.I_k_tot_n_1;
    V_k_series_n_1 = copy_SCN_node.V_k_series_n_1;
    V_kj_dummy_n_1 = copy_SCN_node.V_kj_dummy_n_1;
    V_kj_dummy_tot_n_1=copy_SCN_node.V_kj_dummy_tot_n_1;
    }
    return (*this);
}

 SCN_node::~SCN_node()
{
    //destructor
}

int SCN_node::get_iD()
{
    return this->SCN_id;
}

double SCN_node::get_G(int i)
{
    return this->G[i];
}

double SCN_node::get_R(int i)
{
    return this->R[i];

}

double SCN_node::get_Y(int i)
{
    return this->Y[i];

}

double SCN_node::get_Z(int i)
{
    return this->Z[i];

}
int SCN_node::get_coord(int i)
{
    return(this->SCN_coord[i]);

}
ostream& operator <<(ostream& out,SCN_node& node)
{
    out<<" [ "<< node.SCN_id<<" | "<< node.SCN_coord[0] <<" "<< node.SCN_coord[1] <<" "<<node.SCN_coord[2]<< " ] "<<endl;
    return out;
}

void SCN_node::excitation (int port_id,double V_i)
{
    SCN_ports[port_id].first=V_i;
    //cout<<SCN_ports[port_id].first;
}

void SCN_node::set_neighbour(SCN_node *node)
{
    //this->SCN_neighbours.push_back(node);
    //cout<<*(this->SCN_neighbours[this->SCN_neighbours.size()-1]);
    //cout<<"here"<<*node<<endl;

}

bool SCN_node::check_special_node()
{
    return this->special_node_set_flag;
}

void SCN_node::set_special_SCN_node(bool flag)
{
    this->special_node_set_flag = flag;
}

void SCN_node::print_neighbours()
{
    /*int sizee = this->SCN_neighbours.size();
    for(int i=0; i<sizee; i++)
    {
        cout<< " Node:" << this->SCN_id <<endl;
        cout<< " Neighbours:  " <<endl;
            cout<< *(SCN_neighbours[i])<<endl;
    }*/
}

bool SCN_node::is_PML_SCN()
{
    return (is_PML_node);
}

void SCN_node::set_PEC()
{
    if(!is_PML_node) this->PEC = true;
}

void SCN_node::set_soft_source_node()
{
     this->soft_source_node = true;
}

bool SCN_node::is_PEC()
{
    return (PEC);
}

 double SCN_node::get_reflected(int port_id)
{
    return  SCN_ports[port_id-1].second;
}

 void SCN_node:: set_reflected( int port_id,double Vr)
{
    SCN_ports[port_id-1].second = Vr;
}

void SCN_node::set_PML (int Lz, double Refn_fctr, int c_profile , int scn_typ, int num_PML)
{
    double dl = this-> SCN_dl;

    double sigma = -(e0*0.5*c *log(Refn_fctr))/(1e-3*SCN_dl*pow(num_PML,c_profile+1));

    double sigma_ex(0),sigma_ey(0),sigma_ez(0);

    scn_type = scn_typ;

    is_PML_node = true;

    if( scn_type == 1) // PML digital filter
    {

        sigma_ez = sigma*(pow((Lz+1),c_profile+1)- pow(Lz,c_profile+1));
        //Parameters along the Z axis
        //sigma_e.push_back(sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1)));
        //cout<< "sigma " << sigma_ez <<endl;

        sigma_e[2] = (sigma_ez);
        sigma_h[2] = (sigma_e[2]*u0/e0);
        G[2]=(sigma_e[2]*Z0*1e-3*dl);
        R[2]= (sigma_h[2]*dl*1e-3/Z0);

    }

    else if( scn_type == 2) // PML digital filter
    {
        sigma_ez = sigma*(pow((Lz+1),c_profile+1)- pow(Lz,c_profile+1));

        //cout<< Lz << "  "<< pow((Lz+1),c_profile+1) << " - "<< pow(Lz,c_profile+1)<<endl;
       // cout<< " n = " << c_profile << pow( 2,3)<<endl;
        //cin>> sigma;
        //cout<< "sigma 2..." << sigma <<" ; "<< sigma_ez <<endl;
        double dt = (dl*1e-3/(2*c));
        //Parameters along the Z axis
        //sigma_e.push_back(sigma*(pow((Ly+1),conduct_prof+1)- pow(Ly,conduct_prof+1)));
        sigma_e[2] = (sigma_ez);
        cndtvty_z = sigma_e[2]*dt/e0;

    }




}

 double SCN_node:: get_incidence( int port_id)
{
    return  SCN_ports[port_id-1].first;
}

 void SCN_node:: set_incidence(int port_id,double Vi)
{
    SCN_ports[port_id-1].first = Vi;
}

 double SCN_node:: get_cndtvty_x()
{
    return  this->cndtvty_x;
}

double SCN_node:: get_cndtvty_y()
{
    return  this->cndtvty_y;
}

double SCN_node:: get_cndtvty_z()
{
    return  this->cndtvty_z;
}

vector<double> SCN_node::get_SCN_inc_ports()
{
    vector<double>temp;
    for ( int i =0; i<12;i++)
        temp.push_back(SCN_ports[i].first);

    return temp;
}

void SCN_node::set_all_to_one()
{
    //SCN_ports.assign(24,make_pair(1,0));

    set_incidence(1,1);
    set_incidence(2,2);
    set_incidence(3,3);
    set_incidence(4,4);
    //set_incidence(3,5);
    //set_incidence(4,-0.000687975398236191);
    //set_incidence(5,0.000230325302996784);
    set_incidence(5,5);
    set_incidence(6,6);
    //set_incidence(7,0.000230325302996784);
    set_incidence(7,7);
    set_incidence(8,8);
    set_incidence(9,9);
    set_incidence(10,10);
    set_incidence(11,11);
    set_incidence(12,12);



    set_incidence(13,1);
    set_incidence(14,2);
    set_incidence(15,2);
    set_incidence(16,4);
    set_incidence(17,4);
    set_incidence(18,6);
    set_incidence(19,6);
    set_incidence(20,8);
    set_incidence(21,8);
    set_incidence(22,10);
    set_incidence(23,10);
    set_incidence(24,12);
   /* cout<<" Y "<<endl;
    cout<<Y[1]<<endl;

    cout<<" Z "<<endl;
    cout<<Z[1]<<endl; */

}

void SCN_node::scattering_source_node(double v_inc , double Rs , double & Vy_prime)
{
    vector<double> vr;
    vr.assign(13,0);

    double Ix = this->compute_Ik(1), Iy  = this->compute_Ik(2) , Iz = this->compute_Ik(3);
    double Vx = this->compute_Vj(1), Vz  = this->compute_Vj(3);

    double Y0 = 1/Z0;
    double vi3 = this->get_incidence(3);
    double vi11 = this->get_incidence(11);
    double vi4 = this->get_incidence(4);
    double vi8 = this->get_incidence(8);

    Vy_prime = ( 2*Y0*Rs*( vi3 + vi11 + vi4 + vi8  )  + v_inc )/(4*Rs*Y0 + 1);

   // cout<<"scattering part "<< endl<< vi3 << " "<< vi11 << " "
    //<<vi4 << " " <<vi8 << "....."<<v_inc<<endl;

    //REFLECTED PULSES!!
    vr[3] =  Vy_prime + Z0*Iz - get_incidence(11);  //3
    vr[11] = Vy_prime - Z0*Iz - get_incidence(3);  //11
    vr[5] =  Vz + Z0*Ix - get_incidence(7);  //5
    vr[7] =  Vz - Z0*Ix - get_incidence(5);  //7

    vr[2] = Vx + Z0*Iy - get_incidence(9);  //2
    vr[9] = Vx - Z0*Iy - get_incidence(2);  //9
    vr[4] = Vy_prime - Z0*Ix - get_incidence(8);  //4
    vr[8] = Vy_prime + Z0*Ix - get_incidence(4);  //8

    vr[6]  = Vz - Z0*Iy - get_incidence(10);  //6
    vr[10] = Vz + Z0*Iy - get_incidence(6); //10
    vr[1] =  Vx - Z0*Iz - get_incidence(12); //1
    vr[12] = Vx + Z0*Iz - get_incidence(1); //12

 //scatter + connect combined here...
    this->set_reflected(3,   vr[3] );  //3  , zy
    this->set_reflected(11, vr[11] );  //11  , zy
    this->set_reflected(5,  vr[5]  );   //5 , xz
    this->set_reflected(7,  vr[7]  );    //7 , xz

    this->set_reflected(2,  vr[2] );    //2 , yx
    this->set_reflected(9,  vr[9] );   //9 ,yx
    this->set_reflected(4,  vr[4] );   //4 , xy
    this->set_reflected(8,  vr[8] );   //8 ,xy

    this->set_reflected(6,  vr[6]  );    //6 , yz
    this->set_reflected(10, vr[10] ); //10 , yz
    this->set_reflected(1,  vr[1]  ); //1 , zx
    this->set_reflected(12, vr[12] ); //12 , zx

}

void SCN_node::print_node_voltages()
{
    cout<<endl;
    for(int i = 1; i<= 24 ; i++)
    {
        cout<<i<<" V_inc = "<<std::setprecision(numeric_limits<double>::digits10+1)
        <<this->SCN_ports[i-1].first << "    V_r = " << this->SCN_ports[i-1].second <<endl;

        if( scn_type == 2)
        {
             //cin>>i;
             //cout<<i<<" V_inc_n_1; = "<<std::setprecision(numeric_limits<double>::digits10+1)<<this->V_inc_n_1[i-1]<<endl;

             //cout<< "   I_j_sh_pml_a_n_1 = " << this->I_j_sh_pml_a_n_1[i-1] <<endl;
        }
        cout<<endl;
    }
/*
    cout<<" V_j_n_1 :" << V_j_n_1[0]<< " "<< V_j_n_1[1] << " "<< V_j_n_1[2] <<endl;
    cout<<" V_j_pml_n_1 : " <<V_j_pml_n_1[0]<< " "<< V_j_pml_n_1[1] << " "<< V_j_pml_n_1[2] <<endl;
    cout<<" V_j_tot_n_1 : "<<V_j_tot_n_1[0]<< " "<< V_j_tot_n_1[1] << " "<< V_j_tot_n_1[2] <<endl;
    cout<<" I_j_sh_pml_n_1 :"<<I_j_sh_pml_n_1[0]<< " "<< I_j_sh_pml_n_1[1] << " "<< I_j_sh_pml_n_1[2] <<endl;
    cout<<" I_k_n_1 :" << I_k_n_1[0]<< " "<< I_k_n_1[1] << " "<< I_k_n_1[2] <<endl;
    cout<<" I_k_pml_n_1 :" << I_k_pml_n_1 [0]<< " "<< I_k_pml_n_1 [1] << " "<< I_k_pml_n_1 [2] <<endl;
    cout<<" I_k_tot_n_1 :"<< I_k_tot_n_1 [0]<< " "<< I_k_tot_n_1 [1] << " "<< I_k_tot_n_1 [2] <<endl;
    cout<<" V_k_series_n_1 :"<< V_k_series_n_1 [0]<< " "<< V_k_series_n_1 [1] << " "<< V_k_series_n_1 [2] <<endl;
    cout<<" V_kj_dummy_n_1 :"<< V_kj_dummy_n_1 [0]<< " "<< V_kj_dummy_n_1 [1] << " "<< V_kj_dummy_n_1 [2] << " "<< V_kj_dummy_n_1 [3]<< " "<< V_kj_dummy_n_1 [4] << " "<< V_kj_dummy_n_1 [5] <<endl;
    cout<<" V_kj_dummy_tot_n_1 : "<< V_kj_dummy_tot_n_1 [0] << "  " << V_kj_dummy_tot_n_1 [1] << "  " << V_kj_dummy_tot_n_1 [2] << "  " << V_kj_dummy_tot_n_1 [3] << "  "
                                    << V_kj_dummy_tot_n_1 [4] << "  " << V_kj_dummy_tot_n_1 [5] << endl;
    */
    cout<<endl<<endl;
}

void SCN_node::scattering()
{

    if(this->PEC) return ;
    if(this->soft_source_node) return;

    //this->print_node_voltages();
    //cout<<"......."<<endl;
    //this->set_all_to_one();
/*
    double Yx(this->Y[0]),Zx(this->Z[0]);       // note: this is assuming the grid dimensions are all equal to dl
    double Yy(this->Y[1]),Zy(this->Z[1]);
    double Yz(this->Y[2]),Zz(this->Z[2]);

    double Gxy (this->G[1]), Gxz(this->G[2]) , Gyx(this->G[0]), Gyz(this->G[2]), Gzx(this->G[0]), Gzy(this->G[1]);
    double Rxy (this->R[1]), Rxz(this->R[2]) , Ryx(this->R[0]), Ryz(this->R[2]), Rzx(this->R[0]), Rzy(this->R[1]);


    double Exy =  2*(get_incidence(1)+get_incidence(12) + (Yx+2)*get_incidence(13)-2*get_incidence(14) );
    Exy = Exy/ ( Yx+Gxy+4);

    double Exz =  2*(get_incidence(2)+get_incidence(9) + (Yx+2)*get_incidence(14)-2*get_incidence(13) );
    Exz = Exz/ ( Yx+Gxz+4);

    double Eyx =  2*(get_incidence(3)+get_incidence(11) + (Yy+2)*get_incidence(15)-2*get_incidence(16) );
    Eyx = Eyx/ ( Yy+Gyx+4);

    double Eyz =  2*(get_incidence(4)+ get_incidence(8) + (Yy+2)*get_incidence(16)-2*get_incidence(15) );
    Eyz = Eyz/ ( Yy+Gyz+4);

    double Ezx =  2*(get_incidence(6)+ get_incidence(10) + (Yz+2)*get_incidence(17)-2*get_incidence(18) );
    Ezx = Ezx/ ( Yz + Gzx+4);

    double Ezy =  2*(get_incidence(5)+ get_incidence(7) + (Yz+2)*get_incidence(18)-2*get_incidence(17) );
    Ezy = Ezy/ ( Yz + Gzy+4);

    double Hxy =  2*(get_incidence(5)- get_incidence(7) + (1+2/Zx)*get_incidence(19)-2*get_incidence(20)/Zx );
    Hxy = Hxy/ ( Zx+Rxy+4);

    double Hxz =  2*(get_incidence(8)- get_incidence(4) + (1+2/Zx)*get_incidence(20)-2*get_incidence(19)/Zx );
    Hxz = Hxz/ ( Zx+Rxz+4);

    double Hyx =  2*(get_incidence(10)- get_incidence(6) + (1+2/Zy)*get_incidence(21)-2*get_incidence(22)/Zy );
    Hyx = Hyx/ ( Zy+Ryx+4);

    double Hyz =  2*(get_incidence(2)- get_incidence(9) + (1+2/Zy)*get_incidence(22)-2*get_incidence(21)/Zy);
    Hyz = Hyz/ ( Zy+Ryz+4);

    double Hzx =  2*(get_incidence(3)- get_incidence(11) + (1+2/Zz)*get_incidence(23)-2*get_incidence(24)/Zz);
    Hzx = Hzx/ ( Zz+Rzx+4);

    double Hzy =  2*(get_incidence(12)- get_incidence(1) + (1+2/Zz)*get_incidence(24)-2*get_incidence(23)/Zz);
    Hzy = Hzy/ ( Zz+Rzy+4);


    double Vr1 = Exy + Exz + Hzx + Hzy - get_incidence(12);
    set_reflected(1,Vr1);

    double Vr2 = Exy + Exz - Hyx - Hyz - get_incidence(9);
    set_reflected(2,Vr2);

    double Vr3 = Eyx + Eyz - Hzx - Hzy - get_incidence(11);
    set_reflected(3,Vr3);

    double Vr4 = Eyx + Eyz + Hxy + Hxz - get_incidence(8);
    set_reflected(4,Vr4);

    double Vr5 = Ezx + Ezy - Hxy - Hxz - get_incidence(7);
    set_reflected(5,Vr5);

    double Vr6 = Ezx + Ezy + Hyx + Hyz - get_incidence(10);
    set_reflected(6,Vr6);

    double Vr7 = Ezx + Ezy + Hxy + Hxz - get_incidence(5);
    set_reflected(7,Vr7);

    double Vr8 = Eyx + Eyz - Hxy - Hxz - get_incidence(4);
    set_reflected(8,Vr8);

    double Vr9 = Exy + Exz + Hyx + Hyz - get_incidence(2);
    set_reflected(9,Vr9);

    double Vr10 = Ezx + Ezy - Hyx - Hyz - get_incidence(6);
    set_reflected(10,Vr10);

    double Vr11 = Eyx + Eyz + Hzx + Hzy - get_incidence(3);
    set_reflected(11,Vr11);

    double Vr12 = Exy + Exz - Hzx - Hzy - get_incidence(1);
    set_reflected(12,Vr12);

    double  Vr13 = Exy - get_incidence(13);
    set_reflected(13,Vr13);

    double Vr14 = Exz - get_incidence(14);
    set_reflected(14,Vr14);

    double Vr15 = Eyx - get_incidence(15);
    set_reflected(15,Vr15);

    double Vr16 = Eyz - get_incidence(16);
    set_reflected(16,Vr16);

    double Vr17 = Ezx - get_incidence(17);
    set_reflected(17,Vr17);

    double Vr18 = Ezy - get_incidence(18);
    set_reflected(18,Vr18);

    double Vr19 = -Zx*Hxy + get_incidence(19);
    set_reflected(19,Vr19);

    double Vr20 = -Zx*Hxz + get_incidence(20);
    set_reflected(20,Vr20);

    double Vr21 = -Zy*Hyx + get_incidence(21);
    set_reflected(21,Vr21);

    double Vr22 = -Zy*Hyz + get_incidence(22);
    set_reflected(22,Vr22);

    double Vr23 = -Zz*Hzx + get_incidence(23);
    set_reflected(23,Vr23);

    double Vr24 = -Zz*Hzy + get_incidence(24);
    set_reflected(24,Vr24);

/*
    //STUB CONNECTION
    //STUB CONNECTION IS PERFORMED HERE SINCE IT DOES NOT AFFECT ADJACENT NODES.
    //set_incidence(12,get_reflected(12));
    set_incidence(13,get_reflected(13));
    set_incidence(14,get_reflected(14));
    set_incidence(15,get_reflected(15));
    set_incidence(16,get_reflected(16));
    set_incidence(17,get_reflected(17));
    set_incidence(18,get_reflected(18));
    set_incidence(19,-1*get_reflected(19));
    set_incidence(20,-1*get_reflected(20));
    set_incidence(21,-1*get_reflected(21));
    set_incidence(22,-1*get_reflected(22));
    set_incidence(23,-1*get_reflected(23));
    set_incidence(24,-1*get_reflected(24));
*/

    //this->print_node_voltages();
    //cin>>Exz;


         double Vr1 = get_incidence(2) + get_incidence(3) + get_incidence(9) - get_incidence(11);
        set_reflected(1,0.5*Vr1);

        double Vr2 = get_incidence(1) + get_incidence(6) - get_incidence(10) + get_incidence(12);
        set_reflected(2,0.5*Vr2);

        double Vr3 = get_incidence(1) + get_incidence(4) + get_incidence(8) - get_incidence(12);
        set_reflected(3,0.5*Vr3);

        double Vr4 = get_incidence(3) + get_incidence(5) - get_incidence(7) + get_incidence(11);
        set_reflected(4,0.5*Vr4);

        double Vr5 = get_incidence(4) + get_incidence(6) - get_incidence(8) + get_incidence(10);
        set_reflected(5,0.5*Vr5);

        double Vr6 = get_incidence(2) + get_incidence(5) + get_incidence(7) - get_incidence(9);
        set_reflected(6,0.5*Vr6);

        double Vr7 = -get_incidence(4) + get_incidence(6) + get_incidence(8) + get_incidence(10);
        set_reflected(7,0.5*Vr7);

        double Vr8 = get_incidence(3) - get_incidence(5) + get_incidence(7) + get_incidence(11);
        set_reflected(8,0.5*Vr8);

        double Vr9 = get_incidence(1) - get_incidence(6) + get_incidence(10) + get_incidence(12);
        set_reflected(9,0.5*Vr9);

        double Vr10 =  get_incidence(5) - get_incidence(2)+ get_incidence(7) + get_incidence(9);
        set_reflected(10,0.5*Vr10);

        double Vr11 =  get_incidence(4) + get_incidence(8) + get_incidence(12)- get_incidence(1);
        set_reflected(11,0.5*Vr11);

        double Vr12 = get_incidence(2) - get_incidence(3) + get_incidence(9) + get_incidence(11);
        set_reflected(12,0.5*Vr12);

        //if( this->SCN_id==37001 ) cout<< "scatter of left node Vr11-> Vr3: "<<Vr11 << " ... " <<endl;

        //cout<<"id :  "<<SCN_id<<endl;//if( this->SCN_id==3700 ) cout<< "scatter of left node Vr11-> Vr3: "<<Vr11 << " ... " <<endl;
}

void SCN_node::scattering_PML2()
{
    //Scattering Matrix Parameters

    //this->print_node_voltages();
    //cout<<"......."<<endl;
    // this->set_all_to_one();

    /*int L = this->PML_L;
    int n = 0;   //  conductivity profile
    double sigma = -(e0*0.5*c *log(R))/(1e-3*dl*pow(N,n+1));
    sigma = 0;
    double sigma_e = sigma*(pow((L+1),n+1)- pow(L,n+1));
    double sigma_h = sigma_e*u0/e0;
    double g = sigma_e*Z0*1e-3*dl;                       // this here assumes the permittivity are same in all direction
    double r_ = sigma_h*dl*1e-3/Z0;                      // this here assumes the permeability are same in all direction
*/
    if(this->PEC==1) return;

    if(this->is_PML_SCN()!=true)
    {
        this->scattering();
    }
    else{
    double Yx(this->Y[0]),Zx(this->Z[0]);       // note: this is assuming the grid dimensions are all equal to dl
    double Yy(this->Y[0]),Zy(this->Z[0]);
    double Yz(this->Y[0]),Zz(this->Z[0]);
    //double sigma_ex (0), sigma_ey(0),sigma_ez(0);
    //double sigma_hx(0),sigma_hy(0),sigma_hz(0);

    double Gxy (this->G[1]), Gxz(this->G[2]) , Gyx(this->G[0]), Gyz(this->G[2]), Gzx(this->G[0]), Gzy(this->G[1]);
    double Rxy (this->R[1]), Rxz(this->R[2]) , Ryx(this->R[0]), Ryz(this->R[2]), Rzx(this->R[0]), Rzy(this->R[1]);

    //double Gxy (g), Gxz(g) , Gyx(g), Gyz(g), Gzx(g), Gzy(g);
    //double Rxy (r_), Rxz(r_) , Ryx(r_), Ryz(r_), Rzx(r_), Rzy(r_);

    //COMPUTING VR1
    double a1 = 2/(Yx + Gxy + 4) - 2/(Zz + Rzy + 4);
    double b1 = 2/(Yx + Gxz + 4);
    double c1 = 2/(Yx + Gxy + 4) + 2/(Zz + Rzy +4) - 1;
    double d1 = 2/(Zz + Rzx + 4);
    double e1 = 2*(Yx + 2)/(Yx + Gxy + 4) - 2*( 2/(Yx + Gxz + 4));
    double f1 = (Yx + 2)*2/(Yx + Gxz +4) - 4/(Yx + Gxy +4);
    double g1 = (2*(Zz +2)/(Zz *(Zz + Rzx + 4))) - (2*(2/(Zz + Rzy + 4))/Zz);
    double h1 = (( (Zz + 2)/(Zz) ) * ( 2 /(Zz + Rzy +4))) - 4 / (Zz *(Zz +Rzx +4));
    double Vr1 = a1*get_incidence(1) + b1*get_incidence(2) + d1*get_incidence(3) + b1*get_incidence(9);
           Vr1 = Vr1 - d1*get_incidence(11) + c1*get_incidence(12) + e1*get_incidence(13) + f1*get_incidence(14);
           Vr1 = Vr1 + g1*get_incidence(23) + h1*get_incidence(24);

//COMPUTING VR2
    double a2 = 2/(Yx + Gxz + 4) - 2/(Zy + Ryz + 4);
    double b2 = 2/(Yx + Gxz + 4);                                       //same as a1
    double c2 = 2/(Yx + Gxz + 4) + 2/(Zy + Ryz +4) - 1;
    double d2 = 2/(Zy + Ryx + 4);                                       //same as d1
    double e2 = 2*(Yx + 2)/(Yx + Gxy + 4) - 2*( 2/(Yx + Gxz + 4));      // same as e1
    double f2 = (((Yx + 2)*2) /(Yx + Gxz +4)) - (4/(Yx + Gxy +4));      //same as f1
    double g2 = (2*(Zy +2)/(Zy *(Zy + Ryx + 4))) - (2*(2/(Zy + Ryz + 4))/Zy);
    double h2 = (( (Zy + 2)/(Zy) ) * ( 2 /(Zy + Ryz +4))) - (4 / (Zy *(Zy +Ryx +4)));

    double Vr2 = b2*get_incidence(1) + a2*get_incidence(2) + d2*get_incidence(6) + c2*get_incidence(9);
    Vr2 = Vr2- d2*get_incidence(10) + b2*get_incidence(12) + e2*get_incidence(13) + f2*get_incidence(14);
    Vr2 = Vr2 -g2*get_incidence(21) - h2*get_incidence(22);

  //COMPUTING VR3
    double a3 = 2/(Yy + Gyx + 4) - 2/(Zz + Rzx + 4);
    double b3 = 2/(Yy + Gyz + 4);
    double c3 = 2/(Yy + Gyx + 4) + 2/(Zz + Rzx +4) - 1;
    double d3 = 2/(Zz + Rzy + 4);
    double e3 = 2*(Yy + 2)/(Yy + Gyx + 4) - 2*( 2/(Yy + Gyz + 4));
    double f3 = (((Yy + 2)*2) /(Yy + Gyz +4)) - (4/(Yy + Gyx+4));
    double g3 = (2*(Zz +2)/(Zz *(Zz + Rzx + 4))) - (2*(2/(Zz + Rzy + 4))/Zz);
    double h3 = (( (Zz + 2)/(Zz) ) * ( 2 /(Zz + Rzy +4))) - (4 / (Zz *(Zz +Rzx +4)));

    double Vr3 = d3*get_incidence(1) + a3*get_incidence(3) + b3*get_incidence(4) + b3*get_incidence(8);
    Vr3 = Vr3+ c3*get_incidence(11) - d3*get_incidence(12) + e3*get_incidence(15) + f3*get_incidence(16);
    Vr3 = Vr3 -g3*get_incidence(23) - h3*get_incidence(24);

      //COMPUTING VR4
    double a4 = 2/(Yy + Gyz + 4) - 2/(Zx + Rxz + 4);
    double b4 = 2/(Yy + Gyx + 4);
    double c4 = 2/(Yy + Gyz + 4) + 2/(Zx + Rxz +4) - 1;
    double d4 = 2/(Zx + Rxy + 4);
    double e4 = e3;//2*(Yy + 2)/(Yy + Gyx + 4) - 2*( 2/(Yy + Gyz + 4));
    double f4 = f3;// (((Yy + 2)*2) /(Yy + Gyz +4)) - (4/(Yy + Gyx+4));
    double g4 = (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h4 = (( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr4 = b4*get_incidence(3) + a4*get_incidence(4) + d4*get_incidence(5) - d4*get_incidence(7);
    Vr4 = Vr4+ c4*get_incidence(8) + b4*get_incidence(11) + e4*get_incidence(15) + f4*get_incidence(16);
    Vr4 = Vr4 + g4*get_incidence(19) + h4*get_incidence(20);


      //COMPUTING VR5
    double a5 = 2/(Yz + Gzy + 4) - 2/(Zx + Rxy + 4);
    double b5 = 2/(Yz + Gzx + 4);
    double c5 = 2/(Yz + Gzy + 4) + 2/(Zx + Rxy +4) - 1;
    double d5 = 2/(Zx + Rxz + 4);
    double e5 = 2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f5 = (((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g5 = g4;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h5 = h4; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr5 = d5*get_incidence(4) + a5*get_incidence(5) + b5*get_incidence(6) + c5*get_incidence(7);
    Vr5 = Vr5  - d5*get_incidence(8) + b5*get_incidence(10) + e5*get_incidence(17) +f5*get_incidence(18);
    Vr5 = Vr5 - g5*get_incidence(19) -h5*get_incidence(20);


      //COMPUTING VR6
    double a6 = 2/(Yz + Gzx + 4) - 2/(Zy + Ryx + 4);
    double b6 = 2/(Yz + Gzy + 4);
    double c6 = 2/(Yz + Gzx + 4) + 2/(Zy + Ryx +4) - 1;
    double d6 = 2/(Zy + Ryz + 4);
    double e6 = e5; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f6 = f5;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g6 = g2;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h6 = h2; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr6 = d6*get_incidence(2) + b6*get_incidence(5) + a6*get_incidence(6) + b6*get_incidence(7);
    Vr6 = Vr6 - d6*get_incidence(9) + c6*get_incidence(10) + e6*get_incidence(17) + f6*get_incidence(18);
    Vr6 = Vr6 + g6*get_incidence(21) +h6*get_incidence(22);

          //COMPUTING VR7
    double a7 = 2/(Yz + Gzy + 4) - 2/(Zx + Rxy + 4);
    double b7 = b5;//2/(Yz + Gzy + 4);
    double c7 = 2/(Yz + Gzy + 4) + 2/(Zx + Rxy +4) - 1;
    double d7 = d5; //2/(Zx + Rxz + 4);
    double e7 = e5; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f7 = f5;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g7 = g4;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h7 = h4; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr7 = -d7*get_incidence(4) + c7*get_incidence(5) + b7*get_incidence(6) + a7*get_incidence(7);
    Vr7 = Vr7 + d7*get_incidence(8) + b7*get_incidence(10) + e7*get_incidence(17) + f7*get_incidence(18);
    Vr7 = Vr7 + g7*get_incidence(19) +h7*get_incidence(20);


              //COMPUTING VR8
    double a8 = 2/(Yy + Gyz + 4) - 2/(Zx + Rxz + 4);
    double b8 = b4;   //2/(Yz + Gzy + 4);
    double c8 = 2/(Yy + Gyz + 4) + 2/(Zx + Rxz +4) - 1;
    double d8 = d4; //2/(Zx + Rxz + 4);
    double e8 = e3; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f8 = f3;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g8 = g4;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h8 = h4; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr8 = b8*get_incidence(3) + c8*get_incidence(4) - d8*get_incidence(5) + d8*get_incidence(7);
    Vr8 = Vr8 + a8*get_incidence(8) + b8*get_incidence(11) + e8*get_incidence(15) + f8*get_incidence(16);
    Vr8 = Vr8 - g8*get_incidence(19) -h8*get_incidence(20);

                  //COMPUTING VR9
    double a9 = 2/(Yy + Gxz + 4) - 2/(Zy + Ryz + 4);
    double b9 = b2;//2/(Yz + Gzy + 4);
    double c9 = 2/(Yx + Gxz + 4) + 2/(Zy + Ryz +4) - 1;
    double d9 = d2; //2/(Zx + Rxz + 4);
    double e9 = e2; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f9 = f1;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g9 = g2;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h9 = h2; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr9 = b9*get_incidence(1) + c9*get_incidence(2) - d9*get_incidence(6) + a9*get_incidence(9);
    Vr9 = Vr9 + d9*get_incidence(10) + b9*get_incidence(12) + e9*get_incidence(13) + f9*get_incidence(14);
    Vr9 = Vr9 + g9*get_incidence(21) +h9*get_incidence(22);

    //COMPUTING VR10
    double a10 = 2/(Yz + Gzx + 4) - 2/(Zy + Ryx + 4);
    double b10 = b6;//2/(Yz + Gzy + 4);
    double c10 = 2/(Yz + Gzx + 4) + 2/(Zy + Ryx +4) - 1;
    double d10 = d6; //2/(Zx + Rxz + 4);
    double e10 = e5; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f10 = f5;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g10 = g2;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h10 = h2; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr10 = -d10*get_incidence(2) + b10*get_incidence(5) + c10*get_incidence(6) + b10*get_incidence(7);
    Vr10 = Vr10 + d10*get_incidence(9) + a10*get_incidence(10) + e10*get_incidence(17) + f10*get_incidence(18);
    Vr10 = Vr10 - g10*get_incidence(21) -h10*get_incidence(22);

    //COMPUTING VR11
    double a11 = 0.5*(a3 + c3 + 1 ) + 0.5*(a3 -c3 - 1); //2/(Yz + Gzx + 4) - 2/(Zy + Ryx + 4);
    double b11 = b3;//2/(Yz + Gzy + 4);
    double c11 = 0.5*(a3 + c3 + 1 ) - 0.5*(a3 -c3 - 1) - 1;//a 2/(Yz + Gzx + 4) + 2/(Zy + Ryx +4) - 1;
    double d11 = d3; //2/(Zx + Rxz + 4);
    double e11 = e3; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f11 = f3;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g11 = g1;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h11 = h1; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr11 = -d11*get_incidence(1) + c11*get_incidence(3) + b11*get_incidence(4) + b11*get_incidence(8);
    Vr11 = Vr11 + a11*get_incidence(11) + d11*get_incidence(12) + e11*get_incidence(15) + f11*get_incidence(16);
    Vr11 = Vr11 + g11*get_incidence(23) +h10*get_incidence(24);

                              //COMPUTING VR12
    double a12 = 0.5*(a1 + c1 + 1 ) + 0.5*(a1 -c1 - 1); //2/(Yz + Gzx + 4) - 2/(Zy + Ryx + 4);
    double b12 = b1;//2/(Yz + Gzy + 4);
    double c12 = 0.5*(a1 + c1 + 1 ) - 0.5*(a1 -c1 - 1) - 1;//a 2/(Yz + Gzx + 4) + 2/(Zy + Ryx +4) - 1;
    double d12 = d1; //2/(Zx + Rxz + 4);
    double e12 = e1; //2*(Yz + 2)/(Yz + Gzx + 4) - 2*( 2/(Yz + Gzy + 4));
    double f12 = f1;//(((Yz + 2)*2) /(Yz + Gzy +4)) - (4/(Yz + Gzx+4));
    double g12 = g3;// (2*(Zx +2)/(Zx *(Zx + Rxy + 4))) - (2*(2/(Zx + Rxz + 4))/Zx);
    double h12 = h3; //(( (Zx + 2)/(Zx) ) * ( 2 /(Zx + Rxz +4))) - (4 / (Zx *(Zx +Rxy +4)));

    double Vr12 = c12*get_incidence(1) + b12*get_incidence(2) - d12*get_incidence(3) + b12*get_incidence(9);
    Vr12 = Vr12 + d12*get_incidence(11) + a12*get_incidence(12) + e12*get_incidence(13) + f12*get_incidence(14);
    Vr12 = Vr12 - g12*get_incidence(23) - h12*get_incidence(24);

    //COMPUTING VR13
    double i13 = 0.5*(a1 + c1 + 1);
    double j13 = ( (( Yx + 2 )*2) /(Yx+ Gxy + 4) ) - 1 ;
    double k13 = -4/(Yx + Gxy + 4);

    double Vr13 =  i13*get_incidence(1) + i13*get_incidence(12) + j13*get_incidence(13)+k13*get_incidence(14);

    //COMPUTING VR14
    double i14 = 0.5*( a2 + c2 + 1);
    double m14 = ( -4/(Yx + Gxz + 4));
    double l14 = (Yx + 2)*2/(Yx + Gxz +4) -1;

    double Vr14 =  i14*get_incidence(2) + i14*get_incidence(9) + m14*get_incidence(13) + l14*get_incidence(14);

    //COMPUTING VR15
    double i15 = 0.5*(a3 + c3 + 1);
    double j15 = ( 2*(Yy + 2)/ ( Yy + Gyz + 4 ) ) - 1;
    double k15 = -4 /(Yy + Gyx + 4);

    double Vr15 = i15*get_incidence(3) + i15*get_incidence(11) + j15*get_incidence(15) + k15*get_incidence(16);


    //COMPUTING VR16
    double i16 = 0.5*(a4 + c4 + 1);
    double m16 = -4/(Yy + Gyz + 4);
    double l16 = (((Yy + 2)*2) /(Yy + Gyz +4)) -1 ;

    double Vr16 = i16*get_incidence(4) + i16*get_incidence(8) + m16*get_incidence(15) + l16*get_incidence(16);

    //COMPUTING VR17
    double i17 = 0.5*(a6 + c6 + 1);
    double j17 = 2*( Yz + 2) /(Yz + Gzx + 4) - 1;
    double k17 = -4/(Yz + Gzx +4);

    double Vr17 = i17*get_incidence(6) + i17*get_incidence(10) + j17*get_incidence(17) + k17*get_incidence(18);

    //COMPUTING VR18
    double i18 = 0.5*(a5 + c5 + 1);
    double m18 = -4/(Yz + Gzy + 4);
    double l18 = (((Yz + 2)*2) /(Yz + Gzy +4)) -1 ;

    double Vr18 = i18*get_incidence(5) + i18*get_incidence(7) + m18*get_incidence(17) + l18*get_incidence(18);

    //COMPUTING VR19
    double n19 = -Zx*0.5*(a7 - c7 - 1);
    double o19 = ( (Zx + 2) *2*Zx) / ( Zx *( Zx + Rxy + 4)) -1 ;
    double p19 = 4*Zx/( Zx *( Zx + Rxy + 4));

    double Vr19 = -n19*get_incidence(5) + n19*get_incidence(7) - o19*get_incidence(19) + p19*get_incidence(20);

    //COMPUTING VR20
    double n20 = 2*Zx/(Zx+Rxz +4);
    double r20 = 4/(Zx +Rxz + 4);
    double q20 = (Zx + 2)*(2/(Zx + Rxz +4)) - 1;

    double Vr20 = n20*get_incidence(4) -  n20*get_incidence(8)  + r20*get_incidence(19) - q20*get_incidence(20);

    //COMPUTING VR21
    double n21 =  2*Zy/( Zy + Ryx + 4);
    double o21 =  (( Zy + 2 )/Zy ) *((2*Zy)/ (Zy + Ryx + 4)) - 1;
    double p21 =  4*Zy/(Zy*(Zy + Rzx + 4));

    double Vr21 = n21*get_incidence(6) -  n21*get_incidence(10)  - o21*get_incidence(21) + p21*get_incidence(22);


    //COMPUTING VR22
    double n22 = 2*Zy/( Zy + Ryz + 4);
    double r22 = 4/( Zy + Ryz + 4);
    double q22 = (Zy + 2)*(2/(Zy + Ryz +4)) - 1;

    double Vr22 = -n22*get_incidence(2) +  n22*get_incidence(9)  + r22*get_incidence(21) - q22*get_incidence(22);

    //COMPUTING VR23
    double n23 = 2*Zz/( Zz + Rzx + 4);
    double o23 = (( Zz + 2 )/Zz ) *((2*Zz)/ (Zz + Rzx + 4)) - 1;
    double p23 = (2/Zz)*((2*Zz)/ (Zz + Rzx + 4));

    double Vr23 = -n23*get_incidence(3) +  n23*get_incidence(11)  - o23*get_incidence(23) + p23*get_incidence(24);

    //COMPUTING VR24
    double n24 =  2*Zz/( Zz + Rzy + 4);
    double r24 = 4/( Zz + Rzy + 4);
    double q24 = (Zz + 2)*(2/(Zz + Rzy +4)) - 1;

    double Vr24 = n24*get_incidence(1) -  n24*get_incidence(12)  + r24*get_incidence(23) - q24*get_incidence(24);

    //SET REFLECTED VOLTAGE
    set_reflected(1,Vr1);
    set_reflected(2,Vr2);
    set_reflected(3,Vr3);
    set_reflected(4,Vr4);
    set_reflected(5,Vr5);
    set_reflected(6,Vr6);
    set_reflected(7,Vr7);
    set_reflected(8,Vr8);
    set_reflected(9,Vr9);
    set_reflected(10,Vr10);
    set_reflected(11,Vr11);
    set_reflected(12,Vr12);
    set_reflected(13,Vr13);
    set_reflected(14,Vr14);
    set_reflected(15,Vr15);
    set_reflected(16,Vr16);
    set_reflected(17,Vr17);
    set_reflected(18,Vr18);
    set_reflected(19,Vr19);
    set_reflected(20,Vr20);
    set_reflected(21,Vr21);
    set_reflected(22,Vr22);
    set_reflected(23,Vr23);
    set_reflected(24,Vr24);

    //STUB CONNECTION
    //STUB CONNECTION IS PERFORMED HERE SINCE IT DOES NOT AFFECT ADJACENT NODES.
    //set_incidence(12,get_reflected(12));
   /* set_incidence(13,get_reflected(13));            // capacitive stubs become  open circuits
    set_incidence(14,get_reflected(14));
    set_incidence(15,get_reflected(15));
    set_incidence(16,get_reflected(16));
    set_incidence(17,get_reflected(17));
    set_incidence(18,get_reflected(18));
    set_incidence(19,-1*get_reflected(19));    //inductive stubs become  short circuits
    set_incidence(20,-1*get_reflected(20));
    set_incidence(21,-1*get_reflected(21));
    set_incidence(22,-1*get_reflected(22));
    set_incidence(23,-1*get_reflected(23));
    set_incidence(24,-1*get_reflected(24));

    //this->print_node_voltages();
    //cin>>a1;
*/
    }
}


void SCN_node::scattering_PML()
{

    //this->print_node_voltages();
    //cout<<"....PML..."<<endl;
    //this->set_all_to_one();
/*
    int L = this->PML_L;
    int n = 2;   //  conductivity profile
    double sigma = -(e0*0.5*c *log(R))/(1e-3*dl*pow(N,n+1));
    sigma = 0;
    double sigma_e = sigma*(pow((L+1),n+1)- pow(L,n+1));
    double sigma_h = sigma_e*u0/e0;
    double g = sigma_e*Z0*1e-3*dl;
    double r_ = sigma_h*dl*1e-3/Z0;*/
    if(this->PEC==1) return;

    if(this->is_PML_SCN()!=true)
    {
        //cout<<"here"<<endl;
        this->scattering();
    }
    else
    {



    double Gxy (this->G[1]), Gxz(this->G[2]), Gyx(this->G[0]), Gyz(this->G[2]), Gzx(this->G[0]), Gzy(this->G[1]);
    double Rxy (this->R[1]), Rxz(this->R[2]), Ryx(this->R[0]), Ryz(this->R[2]), Rzx(this->R[0]), Rzy(this->R[1]);

    double Yx(this->Y[0]),Zx(this->Z[0]);       // note: this is assuming the grid dimensions are all equal to dl
    double Yy(this->Y[0]),Zy(this->Z[0]);
    double Yz(this->Y[0]),Zz(this->Z[0]);

    //if(this->SCN_id == 25263) {cout<<Gyz<<" "<<Rxz <<endl; this->print_node_voltages();}

    double Exy =  2*(get_incidence(1)+get_incidence(12) + (Yx+2)*get_incidence(13)-2*get_incidence(14) );
    Exy = Exy/ ( Yx+Gxy+4);

    double Exz =  2*(get_incidence(2)+get_incidence(9) + (Yx+2)*get_incidence(14)-2*get_incidence(13) );
    Exz = Exz/ ( Yx+Gxz+4);

    double Eyx =  2*(get_incidence(3)+get_incidence(11) + (Yy+2)*get_incidence(15)-2*get_incidence(16) );
    Eyx = Eyx/ ( Yy+Gyx+4 );

    double Eyz =  2*(get_incidence(4)+ get_incidence(8) + (Yy+2)*get_incidence(16)-2*get_incidence(15) );
    Eyz = Eyz/ ( Yy+Gyz+4 );

    double Ezx =  2*(get_incidence(6)+ get_incidence(10) + (Yz+2)*get_incidence(17)-2*get_incidence(18) );
    Ezx = Ezx/ ( Yz + Gzx+4);

    double Ezy =  2*(get_incidence(5)+ get_incidence(7) + (Yz+2)*get_incidence(18)-2*get_incidence(17) );
    Ezy = Ezy/ ( Yz + Gzy+4);

    double Hxy =  2*(get_incidence(5)- get_incidence(7) + (1+2/Zx)*get_incidence(19)-2*get_incidence(20)/Zx );
    Hxy = Hxy/ ( Zx+Rxy+4);

    double Hxz =  2*(get_incidence(8)- get_incidence(4) + (1+2/Zx)*get_incidence(20)-2*get_incidence(19)/Zx );
    Hxz = Hxz/ ( Zx+Rxz+4);

    double Hyx =  2*(get_incidence(10)- get_incidence(6) + (1+2/Zy)*get_incidence(21)-2*get_incidence(22)/Zy );
    Hyx = Hyx/ ( Zy+Ryx+4);

    double Hyz =  2*(get_incidence(2)- get_incidence(9) + (1+2/Zy)*get_incidence(22)-2*get_incidence(21)/Zy);
    Hyz = Hyz/ ( Zy+Ryz+4);

    double Hzx =  2*(get_incidence(3)- get_incidence(11) + (1+2/Zz)*get_incidence(23)-2*get_incidence(24)/Zz);
    Hzx = Hzx/ ( Zz+Rzx+4);

    double Hzy =  2*(get_incidence(12)- get_incidence(1) + (1+2/Zz)*get_incidence(24)-2*get_incidence(23)/Zz);
    Hzy = Hzy/ ( Zz+Rzy+4);

    double Vr1 = Exy + Exz + Hzx + Hzy - get_incidence(12);
    set_reflected(1,Vr1);

    double Vr2 = Exy + Exz - Hyx - Hyz - get_incidence(9);
    set_reflected(2,Vr2);

    double Vr3 = Eyx + Eyz - Hzx - Hzy - get_incidence(11);
    set_reflected(3,Vr3);

    double Vr4 = Eyx + Eyz + Hxy + Hxz - get_incidence(8);
    //Vr4 = 1e-15*Vr4;
    set_reflected(4,Vr4);

    double Vr5 = Ezx + Ezy - Hxy - Hxz - get_incidence(7);
    set_reflected(5,Vr5);

    double Vr6 = Ezx + Ezy + Hyx + Hyz - get_incidence(10);
    set_reflected(6,Vr6);

    double Vr7 = Ezx + Ezy + Hxy + Hxz - get_incidence(5);
    set_reflected(7,Vr7);

    double Vr8 = Eyx + Eyz - Hxy - Hxz - get_incidence(4);
    set_reflected(8,Vr8);

    double Vr9 = Exy + Exz + Hyx + Hyz - get_incidence(2);
    set_reflected(9,Vr9);

    double Vr10 = Ezx + Ezy - Hyx - Hyz - get_incidence(6);
    set_reflected(10,Vr10);

    double Vr11 = Eyx + Eyz + Hzx + Hzy - get_incidence(3);
    set_reflected(11,Vr11);

    double Vr12 = Exy + Exz - Hzx - Hzy - get_incidence(1);
    set_reflected(12,Vr12);

    double  Vr13 = Exy - get_incidence(13);
    set_reflected(13,Vr13);

    double Vr14 = Exz - get_incidence(14);
    set_reflected(14,Vr14);

    double Vr15 = Eyx - get_incidence(15);
    set_reflected(15,Vr15);

    double Vr16 = Eyz - get_incidence(16);
    set_reflected(16,Vr16);

    double Vr17 = Ezx - get_incidence(17);
    set_reflected(17,Vr17);

    double Vr18 = Ezy - get_incidence(18);
    set_reflected(18,Vr18);

    double Vr19 = -Zx*Hxy + get_incidence(19);
    set_reflected(19,Vr19);

    double Vr20 = -Zx*Hxz + get_incidence(20);
    set_reflected(20,Vr20);

    double Vr21 = -Zy*Hyx + get_incidence(21);
    set_reflected(21,Vr21);

    double Vr22 = -Zy*Hyz + get_incidence(22);
    set_reflected(22,Vr22);

    double Vr23 = -Zz*Hzx + get_incidence(23);
    set_reflected(23,Vr23);

    double Vr24 = -Zz*Hzy + get_incidence(24);
    set_reflected(24,Vr24);

    //if(this->SCN_id == 25263) {this->print_node_voltages(); cout<< " E ..." << Eyx <<" "<< Eyz << " " << Hxy << " " << Hxz <<endl; }

    //this->print_node_voltages();
   /* cout<<"Ey..."<<Eyx+Eyz<< " "<<endl;
    cout<<"Ex..."<<Exy+Exz<< " "<<endl;
    cout<<"Ez..."<<Ezy+Ezx<< " "<<endl;
    cout<<"Hy..."<<Hyx+Hyz<< " "<<endl;
    cout<<"Hx..."<<Hxy+Hxz<< " "<<endl;
    cout<<"Hz..."<<Hzy+Hzx<< " "<<endl;
*/

    //STUB CONNECTION
    //STUB CONNECTION IS PERFORMED HERE SINCE IT DOES NOT AFFECT ADJACENT NODES.
    set_incidence(13,get_reflected(13));            // capacitive stubs become  open circuits
    set_incidence(14,get_reflected(14));
    set_incidence(15,get_reflected(15));
    set_incidence(16,get_reflected(16));
    set_incidence(17,get_reflected(17));
    set_incidence(18,get_reflected(18));
    set_incidence(19,-1*get_reflected(19));    //inductive stubs become  short circuits
    set_incidence(20,-1*get_reflected(20));
    set_incidence(21,-1*get_reflected(21));
    set_incidence(22,-1*get_reflected(22));
    set_incidence(23,-1*get_reflected(23));
    set_incidence(24,-1*get_reflected(24));

   //cout<<".........."<<endl;

   //cin>>Exz;
    }
}

void SCN_node::scattering_EPML()
{
    if(this->PEC==1) return;
    // cout<<"... EPML ...."<<endl;
    // this->set_all_to_one();
    if(this->is_PML_SCN()!=true)
    {
        this->scattering();
    }

    else {
    double Gxy (this->G[1]), Gxz(this->G[2]) , Gyx(this->G[0]), Gyz(this->G[2]), Gzx(this->G[0]), Gzy(this->G[1]);
    double Rxy (this->R[1]), Rxz(this->R[2]) , Ryx(this->R[0]), Ryz(this->R[2]), Rzx(this->R[0]), Rzy(this->R[1]);

    double Axy = 4/(4+ Gxy) ,Axz = 4/(4+ Gxz) ,Ayx = 4/(4+ Gyx), Ayz = 4/(4+ Gyz) ,Azx = 4/(4+ Gzx), Azy = 4/(4+ Gzy);
    double Amxy = 4/(4+ Rxy) ,Amxz = 4/(4+ Rxz) ,Amyx = 4/(4+ Ryx), Amyz = 4/(4+ Ryz) ,Amzx = 4/(4+ Rzx), Amzy = 4/(4+ Rzy);

    double Yx(0),Zx(0);       // note: this is assuming the grid dimensions are all equal to dl
    double Yy(0),Zy(0);
    double Yz(0),Zz(0);

    Yx = 4*( alpha_x - 0.5);   // directional
    Zx = Yx;

    Yy = 4*( alpha_y - 0.5);
    Zy = Yy;

    Yz = 4*( alpha_z - 0.5);
    Zz = Yz;

    float Cx = 0.5/alpha_x;
    float Dx = Cx;

    float Cy = 0.5/alpha_y;
    float Dy = Cy;

    float Cz = 0.5/alpha_z;
    float Dz = Cz;

    double Exy =  Axy*Cy*(get_incidence(1)+get_incidence(12) + Yy*get_incidence(13)-2*get_incidence(14) );

    double Exz =  Axz*Cz*(get_incidence(2)+get_incidence(9)  + Yz*get_incidence(14)-2*get_incidence(13) );

    double Eyx =  Ayx*Cx*(get_incidence(3)+get_incidence(11) + Yx*get_incidence(15)-2*get_incidence(16) );

    double Eyz =  Ayz*Cz*(get_incidence(4)+ get_incidence(8) + Yz*get_incidence(16)-2*get_incidence(15) );

    double Ezx =  Azx*Cx*(get_incidence(6)+ get_incidence(10) + Yx*get_incidence(17)-2*get_incidence(18) );

    double Ezy =  Azy*Cy*(get_incidence(5)+ get_incidence(7) + Yy*get_incidence(18)-2*get_incidence(17) );

/*
    double Hxy =  2*(get_incidence(5)- get_incidence(7) + (1+2/Zx)*get_incidence(19)-2*get_incidence(20)/Zx );
    Hxy = Hxy/ ( Zx+Rxy+4);

    double Hxz =  2*(get_incidence(8)- get_incidence(4) + (1+2/Zx)*get_incidence(20)-2*get_incidence(19)/Zx );
    Hxz = Hxz/ ( Zx+Rxz+4);

    Zx = -1.6e-16;
    cout<<" ?? "<< -2*get_incidence(19)/Zx +(1+2/Zx)*get_incidence(19) -2*get_incidence(20)/Zx +(1+2/Zx)*get_incidence(20) << "..."<<endl;
    cin>>Hxz;
*/
    double Hxy =  Amxy*Dy*(get_incidence(5)- get_incidence(7) + Zy*get_incidence(19)-2*get_incidence(20));

    double Hxz =  Amxz*Dz*(get_incidence(8)- get_incidence(4) + Zz*get_incidence(20)-2*get_incidence(19) );

    double Hyx =  Amyx*Dx*(get_incidence(10)- get_incidence(6) + Zx*get_incidence(21)-2*get_incidence(22) );

    double Hyz =  Amyz*Dz*(get_incidence(2)- get_incidence(9) + Zz*get_incidence(22)-2*get_incidence(21));

    double Hzx =  Amzx*Dx*(get_incidence(3)- get_incidence(11) + Zx*get_incidence(23)-2*get_incidence(24));

    double Hzy =  Amzy*Dy*(get_incidence(12)- get_incidence(1) + Zy*get_incidence(24)-2*get_incidence(23));

    double Vr1 = Exy + Exz + Hzx + Hzy - get_incidence(12);
    set_reflected(1,Vr1);

    double Vr2 = Exy + Exz - Hyx - Hyz - get_incidence(9);
    set_reflected(2,Vr2);

    double Vr3 = Eyx + Eyz - Hzx - Hzy - get_incidence(11);
    set_reflected(3,Vr3);

    double Vr4 = Eyx + Eyz + Hxy + Hxz - get_incidence(8);
    set_reflected(4,Vr4);

    double Vr5 = Ezx + Ezy - Hxy - Hxz - get_incidence(7);
    set_reflected(5,Vr5);

    double Vr6 = Ezx + Ezy + Hyx + Hyz - get_incidence(10);
    set_reflected(6,Vr6);

    double Vr7 = Ezx + Ezy + Hxy + Hxz - get_incidence(5);
    set_reflected(7,Vr7);

    double Vr8 = Eyx + Eyz - Hxy - Hxz - get_incidence(4);
    set_reflected(8,Vr8);

    double Vr9 = Exy + Exz + Hyx + Hyz - get_incidence(2);
    set_reflected(9,Vr9);

    double Vr10 = Ezx + Ezy - Hyx - Hyz - get_incidence(6);
    set_reflected(10,Vr10);

    double Vr11 = Eyx + Eyz + Hzx + Hzy - get_incidence(3);
    set_reflected(11,Vr11);

    double Vr12 = Exy + Exz - Hzx - Hzy - get_incidence(1);
    set_reflected(12,Vr12);

    double  Vr13 = Exy - get_incidence(13);
    set_reflected(13,Vr13);

    double Vr14 = Exz - get_incidence(14);
    set_reflected(14,Vr14);

    double Vr15 = Eyx - get_incidence(15);
    set_reflected(15,Vr15);

    double Vr16 = Eyz - get_incidence(16);
    set_reflected(16,Vr16);

    double Vr17 = Ezx - get_incidence(17);
    set_reflected(17,Vr17);

    double Vr18 = Ezy - get_incidence(18);
    set_reflected(18,Vr18);

    double Vr19 = Hxy - get_incidence(19);
    set_reflected(19,Vr19);

    double Vr20 = Hxz - get_incidence(20);
    set_reflected(20,Vr20);

    double Vr21 = Hyx - get_incidence(21);
    set_reflected(21,Vr21);

    double Vr22 = Hyz - get_incidence(22);
    set_reflected(22,Vr22);

    double Vr23 = Hzx - get_incidence(23);
    set_reflected(23,Vr23);

    double Vr24 = Hzy - get_incidence(24);
    set_reflected(24,Vr24);

    /*if ( SCN_id==1){

    cout<<" Y: "<< Yx<<" "<<Yy<<" "<<Yz<<" "<<Dx<<" "<<Dy<<" "<<Dz<<endl;
    cout<<" A: "<< Axy <<" "<< Axz << " "<< Ayz << " " << Ayx << " " << Azx << " "<< Azy<<endl;
    cout<<" Am: "<< Amxy <<" "<< Amxz << " "<< Amyz << " " << Amyx << " " << Amzx << " "<< Amzy<<endl;
    //cout<<Vr3<<"  "<<Vr11<<endl;

    //this->print_node_voltages();
    cout<<endl;
    cout<<"Ey..."<<Eyx+Eyz<< " "<<endl;
    cout<<"Ex..."<<Exy+Exz<< " "<<endl;
    cout<<"Ez..."<<Ezy+Ezx<< " "<<endl;
    cout<<"Hy..."<<Hyx+Hyz<< " "<<endl;
    cout<<"Hx..."<<Hxy+Hxz<< " "<<endl;
    cout<<"Hz..."<<Hzy+Hzx<< " "<<endl;
    //cin>>Exz;
    }*/

    //STUB CONNECTION
    //STUB CONNECTION IS PERFORMED HERE SINCE IT DOES NOT AFFECT ADJACENT NODES.
    set_incidence(13,get_reflected(13));            // capacitive stubs become  open circuits
    set_incidence(14,get_reflected(14));
    set_incidence(15,get_reflected(15));
    set_incidence(16,get_reflected(16));
    set_incidence(17,get_reflected(17));
    set_incidence(18,get_reflected(18));
    set_incidence(19,1*get_reflected(19));    //inductive stubs become  short circuits
    set_incidence(20,1*get_reflected(20));
    set_incidence(21,1*get_reflected(21));
    set_incidence(22,1*get_reflected(22));
    set_incidence(23,1*get_reflected(23));
    set_incidence(24,1*get_reflected(24));

   //cout<<".........."<<endl;
    }
}



void SCN_node::scattering_PML_df(bool check)
{
     if(this->is_PML_SCN()!=true)
    {
        //cout<<"here"<<endl;
        this->scattering();

    }

    else
    {
    //if( check == true) { cout<<" before scatter : "<<endl<<endl; this->print_node_voltages();}
    //if(this->SCN_id == 25263) this->print_node_voltages();

    vector<long double> vr;
    vr.assign(13,0);
    double L_x =  1/(Z0*( 4 + this->cndtvty_y + this->cndtvty_z) );
    double L_y =  1/(Z0*( 4 + this->cndtvty_x + this->cndtvty_z) );
    double L_z =  1/(Z0*( 4 + this->cndtvty_y + this->cndtvty_x) );
    double Ix = this->compute_Ik(1), Iy  = this->compute_Ik(2) , Iz = this->compute_Ik(3);
    double Vx = this->compute_Vj(1), Vy  = this->compute_Vj(2) , Vz = this->compute_Vj(3);
    double alpha_x = this->compute_alpha(1), alpha_y = this->compute_alpha(2) , alpha_z = this->compute_alpha(3);

    int test= 0;

    //REFLECTED PULSES!!
    vr[3] = alpha_y*Vy + Z0*alpha_z*Iz - get_incidence(11);  //3
    vr[11] = alpha_y*Vy - Z0*alpha_z*Iz - get_incidence(3);  //11
    vr[5] = alpha_z*Vz + Z0*alpha_x*Ix - get_incidence(7);  //5
    vr[7] = alpha_z*Vz - Z0*alpha_x*Ix - get_incidence(5);  //7

    vr[2] = alpha_x*Vx + Z0*alpha_y*Iy - get_incidence(9);  //2
    vr[9] = alpha_x*Vx - Z0*alpha_y*Iy - get_incidence(2);  //9
    vr[4] = alpha_y*Vy - Z0*alpha_x*Ix - get_incidence(8);  //4
    vr[8] = alpha_y*Vy + Z0*alpha_x*Ix - get_incidence(4);  //8

    vr[6]  = alpha_z*Vz - Z0*alpha_y*Iy - get_incidence(10);  //6
    vr[10] = alpha_z*Vz + Z0*alpha_y*Iy - get_incidence(6); //10
    vr[1] = alpha_x*Vx - Z0*alpha_z*Iz - get_incidence(12); //1
    vr[12] = alpha_x*Vx + Z0*alpha_z*Iz - get_incidence(1); //12


    //PML related components V_j_PML
    double v_j_pml_x = this->compute_Vj_pml(1);
    double v_j_pml_y = this->compute_Vj_pml(2);
    double v_j_pml_z = this->compute_Vj_pml(3);

    //storing for the next time step
    V_j_tot_n_1[0] = alpha_x*Vx + v_j_pml_x;
    V_j_tot_n_1[1] = alpha_y*Vy + v_j_pml_y;
    V_j_tot_n_1[2] = alpha_z*Vz + v_j_pml_z;

    V_j_n_1[0] = Vx;
    V_j_n_1[1] = Vy;
    V_j_n_1[2] = Vz;

    //PML related components V_dummy

    double v_series_x = this->compute_v_series_k(1);
    double v_series_y = this->compute_v_series_k(2);
    double v_series_z = this->compute_v_series_k(3);

    double i_k_pml_x = this->compute_i_k_pml(1, v_series_x);
    double i_k_pml_y = this->compute_i_k_pml(2, v_series_y);
    double i_k_pml_z = this->compute_i_k_pml(3, v_series_z);


//SCN_node::compute_v_dummy_kj(double cndtvty_i,double cndtvty_j,double cndtvty_k, double i_k_pml_k,double Ik, int axis, double v_dummy_k_n_1)
    double v_dummy_yx = this->compute_v_dummy_kj(this->cndtvty_z,  this->cndtvty_x, this->cndtvty_y, i_k_pml_y,Iy,2, V_kj_dummy_tot_n_1[0] );
    double v_dummy_zx = this->compute_v_dummy_kj(this->cndtvty_y,  this->cndtvty_x, this->cndtvty_z, i_k_pml_z,Iz,3, V_kj_dummy_tot_n_1[1] );
    double v_dummy_xy = this->compute_v_dummy_kj(this->cndtvty_z,  this->cndtvty_y, this->cndtvty_x, i_k_pml_x,Ix,1, V_kj_dummy_tot_n_1[2] );
    double v_dummy_zy = this->compute_v_dummy_kj(this->cndtvty_x,  this->cndtvty_y, this->cndtvty_z, i_k_pml_z,Iz,3, V_kj_dummy_tot_n_1[3] );
    double v_dummy_xz = this->compute_v_dummy_kj(this->cndtvty_y,  this->cndtvty_z, this->cndtvty_x, i_k_pml_x,Ix,1, V_kj_dummy_tot_n_1[4] );
    double v_dummy_yz = this->compute_v_dummy_kj(this->cndtvty_x,  this->cndtvty_z, this->cndtvty_y, i_k_pml_y,Iy,2, V_kj_dummy_tot_n_1[5] );

//scatter + connect combined here...
    this->set_reflected(3,  exp(-cndtvty_x) *(vr[3]  + v_j_pml_y + v_dummy_zy) );  //3  , zy
    this->set_reflected(11, exp(-cndtvty_x) *(vr[11] +  v_j_pml_y - v_dummy_zy) );  //11  , zy
    this->set_reflected(5,  exp(-cndtvty_y) *(vr[5]  + v_j_pml_z + v_dummy_xz) ) ;   //5 , xz
    this->set_reflected(7,  exp(-cndtvty_y) *(vr[7]   + v_j_pml_z - v_dummy_xz) );    //7 , xz

    this->set_reflected(2,  exp(-cndtvty_z) *(vr[2] + v_j_pml_x + v_dummy_yx) );    //2 , yx
    this->set_reflected(9,  exp(-cndtvty_z) *(vr[9] + v_j_pml_x - v_dummy_yx) );   //9 ,yx
    this->set_reflected(4,  exp(-cndtvty_z) *(vr[4] + v_j_pml_y - v_dummy_xy) );   //4 , xy
    this->set_reflected(8,  exp(-cndtvty_z) *(vr[8] + v_j_pml_y + v_dummy_xy) );   //8 ,xy

    this->set_reflected(6,  exp(-cndtvty_x) *( vr[6]  + v_j_pml_z - v_dummy_yz ) );    //6 , yz
    this->set_reflected(10, exp(-cndtvty_x) *( vr[10] + v_j_pml_z + v_dummy_yz) ); //10 , yz
    this->set_reflected(1,  exp(-cndtvty_y) *(vr[1]  + v_j_pml_x - v_dummy_zx) ); //1 , zx
    this->set_reflected(12, exp(-cndtvty_y) *(vr[12] + v_j_pml_x + v_dummy_zx) ); //12 , zx

    //store current
    I_k_pml_n_1[0] = i_k_pml_x;
    I_k_pml_n_1[1] = i_k_pml_y;
    I_k_pml_n_1[2] = i_k_pml_z;

    I_k_n_1[0] = Ix;
    I_k_n_1[1] = Iy;
    I_k_n_1[2] = Iz;

    I_k_tot_n_1[0] = 2*L_x*v_series_x + i_k_pml_x;
    I_k_tot_n_1[1] = 2*L_y*v_series_y + i_k_pml_y;
    I_k_tot_n_1[2] = 2*L_z*v_series_z + i_k_pml_z;

    V_k_series_n_1[0] = v_series_x;
    V_k_series_n_1[1] = v_series_y;
    V_k_series_n_1[2] = v_series_z;
    // store Vdummy
    V_kj_dummy_tot_n_1[0] = v_dummy_yx + alpha_y*Z0*Iy ;
    V_kj_dummy_tot_n_1[1] = v_dummy_zx + alpha_z*Z0*Iz ;
    V_kj_dummy_tot_n_1[2] = v_dummy_xy + alpha_x*Z0*Ix ;
    V_kj_dummy_tot_n_1[3] = v_dummy_zy + alpha_z*Z0*Iz ;
    V_kj_dummy_tot_n_1[4] = v_dummy_xz + alpha_x*Z0*Ix ;
    V_kj_dummy_tot_n_1[5] = v_dummy_yz + alpha_y*Z0*Iy ;

    V_kj_dummy_n_1[0] = v_dummy_yx  ;
    V_kj_dummy_n_1[1] = v_dummy_zx;
    V_kj_dummy_n_1[2] = v_dummy_xy ;
    V_kj_dummy_n_1[3] = v_dummy_zy  ;
    V_kj_dummy_n_1[4] = v_dummy_xz  ;
    V_kj_dummy_n_1[5] = v_dummy_yz  ;

   // if(this->SCN_id == 25263) { cout<< exp(-cndtvty_z)<< endl<<endl; this->print_node_voltages();
    //                            cout<<"vr4...."<<(vr[4] + v_j_pml_y - v_dummy_xy)<<endl; }

/*
if( SCN_id == 0 )
{
   bool checks = false;
    cout<< "PML ? " <<this->is_PML_node <<endl;

    cout<< "The scaling applied is " << exp(-cndtvty_z) <<endl<<endl;
    cout<< "The voltages are ......"<<endl;
    this->print_node_voltages();
    cout<<endl<<endl;
    cin>>checks;

}
*/
/*
     if ( check == true )
    {
        cout<<endl<<"After scatter " <<endl<<endl;

        cout<<"alpha_j : " <<alpha_x<<" "<<alpha_y<<" "<<alpha_z<<endl;

        cout<<"L_j*Z0 :   "<<Z0*L_x<<" "<<Z0*L_y<<" "<<Z0*L_z<<" "<<endl;

        cout<<"Vj :   "<<Vx<<" "<<Vy<<" "<<Vz<<endl;

        cout<<"Ik :   "<<Ix<<" "<<Iy<<" "<<Iz<<endl;

        cout<<"Vj_pml :   "<< v_j_pml_x <<" "<<v_j_pml_y <<" "<<v_j_pml_z<<endl;

        cout<<"Vseries_k :   "<< v_series_x <<" "<<v_series_y <<" "<<v_series_z<<endl;

        cout<<"ik_pml :   "<< i_k_pml_x <<" "<<i_k_pml_y <<" "<<i_k_pml_z<<endl;

        cout<<"V dummy kj _pml :   "<< v_dummy_yx <<" "<<v_dummy_zx <<" "<<v_dummy_xy <<" "<<v_dummy_zy<<" "<<v_dummy_xz <<" "<<v_dummy_yz<<endl<<endl;

        cout<<" After scatter :"<<endl; this->print_node_voltages();
*/
    //cin>>test;

//cout<<"here"<<endl;

    //store incident pulses
//this->print_node_voltages();}
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

    }
}

double SCN_node::compute_Vj(int j)
{
    if(j==1)  //x
    {
        return( 0.5* ( this->get_incidence(12)+this->get_incidence(2)+this->get_incidence(9)+this->get_incidence(1) ) );
    }
    else if( j==2) //y
    {
        return( 0.5* ( this->get_incidence(3)+this->get_incidence(8)+this->get_incidence(4)+this->get_incidence(11) ) );
    }
    else if(j==3)  //z
    {
        return( 0.5* ( this->get_incidence(7)+this->get_incidence(6)+this->get_incidence(5)+this->get_incidence(10) ) );
    }

    int k;
    cout<<"failing in vj " <<j <<endl;
    cin>>k;
    return 0;
}

double SCN_node::compute_Ik(int k)
{
    if(k==1)  //x
    {
        return( 0.5* ( this->get_incidence(7)-this->get_incidence(8)+this->get_incidence(4)-this->get_incidence(5) )/Z0 );
    }
    else if( k==2) //y
    {
        return( 0.5* ( this->get_incidence(6)-this->get_incidence(10)+this->get_incidence(9)-this->get_incidence(2) )/Z0 );
    }
    else if(k==3)  //z
    {
        return( 0.5* ( this->get_incidence(1)-this->get_incidence(12)+this->get_incidence(11)-this->get_incidence(3) )/Z0 );
    }
    int j;
    cout<<"failing in ik " <<k<<endl;
    cin>>k;
    return 0;
}

double SCN_node::compute_alpha(int j)
{
    if ( j == 1) return( 8/(    (4 + this->cndtvty_y + this->cndtvty_z )*(2 + this->cndtvty_x) ) );

    else if ( j==2) return( 8/(    (4 + this->cndtvty_x + this->cndtvty_z )*(2 + this->cndtvty_y) ));

    else if ( j==3) return( 8/(    (4 + this->cndtvty_y + this->cndtvty_x )*(2 + this->cndtvty_z) ));

return 0;
}


double SCN_node::compute_Vj_pml(int j)
{
    //this->set_all_to_one();

    double Vx = this->compute_Vj(1), Vy  = this->compute_Vj(2) , Vz = this->compute_Vj(3);

    if(j ==1) //x
    {
        double F_x = 1/(  8   +  2*this->cndtvty_z + 2*this->cndtvty_y   );
        double H_x = 1/(    Z0*(2 + cndtvty_x)  );
        double G_x = (  8   -  2*this->cndtvty_z - 2*this->cndtvty_y   );
        double i_shunt_pml_x = this->compute_i_shunt_pml(1);

        double v_j_pml_x =  F_x*(2*Z0*i_shunt_pml_x + Z0*cndtvty_x*(8*H_x*Vx + i_shunt_pml_x) - (2*Z0 - Z0*cndtvty_x)*(8*H_x*this->V_j_n_1[0] + this->I_j_sh_pml_n_1[0])
                                 + G_x*this->V_j_tot_n_1[0]);

      /*  if( v_j_pml_x != 0 )
        {
            cout<< "id : " <<SCN_id<<endl;
            cout<<" debugging " <<endl;
            cout<<" i shunt pml x : "<< i_shunt_pml_x <<endl;
            cout<<" V_x[n-1] : " << V_j_n_1[0] <<endl;
            cout<<" I shunt pml x[n-1] " << I_j_sh_pml_n_1[0]<<endl;
            cout<<" V tot x : " << V_j_tot_n_1[0]<<endl;
        }
*/
        //store for next time step
        this->I_j_sh_pml_n_1[0] = i_shunt_pml_x;



        return v_j_pml_x;

    }
    else if( j ==2)//y
    {
        double F_y = 1/(  8   +  2*this->cndtvty_x + 2*this->cndtvty_z   );
        double H_y = 1/(    Z0*(2 + cndtvty_y)  );
        double G_y = (  8   -  2*this->cndtvty_z - 2*this->cndtvty_x   );
        double i_shunt_pml_y = this->compute_i_shunt_pml(2);

        double v_j_pml_y =  F_y*(2*Z0*i_shunt_pml_y + Z0*cndtvty_y*(8*H_y*Vy + i_shunt_pml_y) - (2*Z0 - Z0*cndtvty_y)*(8*H_y*this->V_j_n_1[1] + this->I_j_sh_pml_n_1[1])
                                 + G_y*this->V_j_tot_n_1[1]);

       /* if( v_j_pml_y != 0 )
        {
            cout<< "id : " <<SCN_id<<endl;
            cout<<" debugging " <<endl;
            cout<<" i shunt pml y : "<< i_shunt_pml_y <<endl;
            cout<<" V_y[n-1] : " << V_j_n_1[1] <<endl;
            cout<<" I shunt pml y[n-1] " << I_j_sh_pml_n_1[1]<<endl;
            cout<<" V tot y : " << V_j_tot_n_1[1]<<endl;
        }
*/
        //store for next time step
        this->I_j_sh_pml_n_1[1] = i_shunt_pml_y;

        return v_j_pml_y;

    }
    else if(j ==3)//z
    {
        double F_z = 1/(  8   +  2*this->cndtvty_x + 2*this->cndtvty_y   );
        double H_z = 1/(    Z0*(2 + cndtvty_z)  );
        double G_z = (  8   -  2*this->cndtvty_y - 2*this->cndtvty_x   );
        double i_shunt_pml_z = this->compute_i_shunt_pml(3);

        double v_j_pml_z =  F_z*(2*Z0*i_shunt_pml_z + Z0*cndtvty_z*(8*H_z*Vz + i_shunt_pml_z) - (2*Z0 - Z0*cndtvty_z)*(8*H_z*this->V_j_n_1[2] + this->I_j_sh_pml_n_1[2])
                                 + G_z*this->V_j_tot_n_1[2]);

       /* if( v_j_pml_z != 0 )
        {
            cout<< "id : " <<SCN_id<<endl;
            cout<<" debugging " <<endl;
            cout<<" i shunt pml z : "<< i_shunt_pml_z <<endl;
            cout<<" V_z[n-1] : " << V_j_n_1[2] <<endl;
            cout<<" I shunt pml z [n-1] " << I_j_sh_pml_n_1[2]<<endl;
            cout<<" V tot z : " << V_j_tot_n_1[2]<<endl;
        }
*/
         //store for next time step
        this->I_j_sh_pml_n_1[2] = i_shunt_pml_z;

        return v_j_pml_z;
    }
    int k;
    cout<<"failing in v_pml " <<j <<endl;
    cin>>k;
    return 0;

}

double SCN_node::compute_i_shunt_pml(int j)
{
    //update i_shunt_j_pml_n_1_x,y,z here

    if( j == 1)//x
    {
        double H_x = 1/(    Z0*(2 + cndtvty_x)  );
        double i_shunt_pml_12 = this->compute_i_shunt_pml_a(12,H_x ,cndtvty_x) ,
        i_shunt_pml_1 = this->compute_i_shunt_pml_a(1, H_x , cndtvty_x) ,
        i_shunt_pml_9 = this->compute_i_shunt_pml_a(9, H_x, cndtvty_x),
        i_shunt_pml_2 = this->compute_i_shunt_pml_a(2, H_x, cndtvty_x);

        double i_shunt_12 = 4*H_x*this->get_incidence(12),
        i_shunt_1 = 4*H_x*this->get_incidence(1),
        i_shunt_9 = 4*H_x*this->get_incidence(9),
        i_shunt_2 = 4*H_x*this->get_incidence(2);

        //total shunt current in each line for the Vx
        double i_shunt_tot_12 = i_shunt_pml_12 + i_shunt_12;
        double i_shunt_tot_1 = i_shunt_pml_1 + i_shunt_1;
        double i_shunt_tot_9 = i_shunt_pml_9 + i_shunt_9;
        double i_shunt_tot_2 = i_shunt_pml_2 + i_shunt_2;

        //storing for next run
        this->I_j_sh_tot_a_n_1[1-1] = i_shunt_tot_1;
        this->I_j_sh_tot_a_n_1[12-1] = i_shunt_tot_12;
        this->I_j_sh_tot_a_n_1[9-1] = i_shunt_tot_9;
        this->I_j_sh_tot_a_n_1[2-1] = i_shunt_tot_2;

        //store for next run
        this->I_j_sh_pml_a_n_1[1-1] =  i_shunt_pml_1;
        this->I_j_sh_pml_a_n_1[12-1] = i_shunt_pml_12;
        this->I_j_sh_pml_a_n_1[9-1] =  i_shunt_pml_9;
        this->I_j_sh_pml_a_n_1[2-1] =  i_shunt_pml_2;

        //this->I_j_sh_pml_n_1[1-1] = i_shunt_pml_12 + i_shunt_pml_1 + i_shunt_pml_9 + i_shunt_pml_2;

        return(i_shunt_pml_12+i_shunt_pml_1+i_shunt_pml_9+i_shunt_pml_2);

    }
    else if( j == 2)//y
    {
        double H_y = 1/(    Z0*(2 + cndtvty_y)  );

        double i_shunt_pml_3 = this->compute_i_shunt_pml_a(3,H_y ,cndtvty_y) ,
        i_shunt_pml_11 = this->compute_i_shunt_pml_a(11, H_y , cndtvty_y) ,
        i_shunt_pml_8 = this->compute_i_shunt_pml_a(8,H_y, cndtvty_y),
        i_shunt_pml_4 = this->compute_i_shunt_pml_a(4, H_y, cndtvty_y);

        double i_shunt_4 = 4*H_y*this->get_incidence(4),
        i_shunt_3 = 4*H_y*this->get_incidence(3),
        i_shunt_11 = 4*H_y*this->get_incidence(11),
        i_shunt_8 = 4*H_y*this->get_incidence(8);

        //total shunt current in each line for the Vx
        double i_shunt_tot_4 = i_shunt_pml_4 + i_shunt_4;
        double i_shunt_tot_8 = i_shunt_pml_8 + i_shunt_8;
        double i_shunt_tot_11 = i_shunt_pml_11 + i_shunt_11;
        double i_shunt_tot_3 = i_shunt_pml_3 + i_shunt_3;

        //storing for next run
        this->I_j_sh_tot_a_n_1[4-1] = i_shunt_tot_4;
        this->I_j_sh_tot_a_n_1[11-1] = i_shunt_tot_11;
        this->I_j_sh_tot_a_n_1[3-1] = i_shunt_tot_3;
        this->I_j_sh_tot_a_n_1[8-1] = i_shunt_tot_8;

        //store for next run
        this->I_j_sh_pml_a_n_1[4-1] =  i_shunt_pml_4;
        this->I_j_sh_pml_a_n_1[11-1] = i_shunt_pml_11;
        this->I_j_sh_pml_a_n_1[3-1] =  i_shunt_pml_3;
        this->I_j_sh_pml_a_n_1[8-1] =  i_shunt_pml_8;

        //store for next run
        //this->I_j_sh_pml_n_1[2-1] = i_shunt_pml_11 + i_shunt_pml_3 + i_shunt_pml_4 + i_shunt_pml_8;

        return(i_shunt_pml_3+i_shunt_pml_11+i_shunt_pml_8+i_shunt_pml_4);


    }
    else if ( j ==3 )//z
    {
        double H_z = 1/(    Z0*(2 + cndtvty_z)  );

        double i_shunt_pml_6 = this->compute_i_shunt_pml_a(6,H_z ,cndtvty_z) ,
        i_shunt_pml_10 = this->compute_i_shunt_pml_a(10, H_z , cndtvty_z) ,
        i_shunt_pml_5 = this->compute_i_shunt_pml_a(5, H_z, cndtvty_z),
         i_shunt_pml_7 = this->compute_i_shunt_pml_a(7, H_z, cndtvty_z);

        double i_shunt_6 = 4*H_z*this->get_incidence(6),
        i_shunt_10 = 4*H_z*this->get_incidence(10),
         i_shunt_7 = 4*H_z*this->get_incidence(7),
        i_shunt_5 = 4*H_z*this->get_incidence(5);

        //total shunt current in each line for the Vx
        double i_shunt_tot_6 = i_shunt_pml_6 + i_shunt_6;
        double i_shunt_tot_10 = i_shunt_pml_10 + i_shunt_10;
        double i_shunt_tot_5 = i_shunt_pml_5 + i_shunt_5;
        double i_shunt_tot_7 = i_shunt_pml_7 + i_shunt_7;

        //storing for next run
        this->I_j_sh_tot_a_n_1[6-1] = i_shunt_tot_6;
        this->I_j_sh_tot_a_n_1[10-1] = i_shunt_tot_10;
        this->I_j_sh_tot_a_n_1[5-1] = i_shunt_tot_5;
        this->I_j_sh_tot_a_n_1[7-1] = i_shunt_tot_7;

        //store for next run
        this->I_j_sh_pml_a_n_1[6-1]  =  i_shunt_pml_6;
        this->I_j_sh_pml_a_n_1[10-1] = i_shunt_pml_10;
        this->I_j_sh_pml_a_n_1[5-1]  =  i_shunt_pml_5;
        this->I_j_sh_pml_a_n_1[7-1]  =  i_shunt_pml_7;

        //this->I_j_sh_pml_n_1[3-1] = i_shunt_pml_6 + i_shunt_pml_10 + i_shunt_pml_5 + i_shunt_pml_7;


        return(i_shunt_pml_6+i_shunt_pml_10+i_shunt_pml_5+i_shunt_pml_7);

    }
    int k;
    cout<<"failing in ishunt_pml " <<j <<endl;
    cin>>k;
    return 0;
}

double SCN_node ::compute_i_shunt_pml_a(int port_id, double H, double cndtvty)
{
    if( (port_id==5)||(port_id==7)||(port_id==4)||(port_id==8) )
    {
        double i_sh = 2*H*cndtvty_x*this->get_incidence(port_id) - 2*H*(2-cndtvty_x)*this->V_inc_n_1[port_id-1] +

        Z0*H*(2-cndtvty)*this->I_j_sh_tot_a_n_1[port_id-1] ;

        return i_sh;

    }
    else if( (port_id==2)||(port_id==9)||(port_id==10)||(port_id==6) )
    {
        double i_sh = 2*H*cndtvty_y*this->get_incidence(port_id) - 2*H*(2-cndtvty_y)*this->V_inc_n_1[port_id-1] +

        Z0*H*(2-cndtvty)*this->I_j_sh_tot_a_n_1[port_id-1] ;

        return i_sh;

    }
    else if( (port_id==1)||(port_id==12)||(port_id==11)||(port_id==3) )
    {

        double i_sh = 2*H*cndtvty_z*this->get_incidence(port_id) - 2*H*(2-cndtvty_z)*this->V_inc_n_1[port_id-1] +

        Z0*H*(2-cndtvty)*this->I_j_sh_tot_a_n_1[port_id-1] ;

        return i_sh;

    }

    return 0;
}

double  SCN_node::compute_i_k_pml(int k,double v_series_k)
{

    if(k ==1)//x
    {   double v_series_x = v_series_k;
        double J_x = Z0*(4 - this->cndtvty_y - this->cndtvty_z);
        double L_x = 1 / (Z0*(4 + this->cndtvty_y + this->cndtvty_z) );
        double ik_pml = this->cndtvty_x*v_series_x - this->V_k_series_n_1[0]*(2-this->cndtvty_x) + I_k_tot_n_1[0]*J_x;
        ik_pml = L_x*ik_pml;

        return ik_pml;
    }
    else if ( k==2 ) //y
    {
        double v_series_y = v_series_k;
        double J_y = Z0*(4 - this->cndtvty_x - this->cndtvty_z);
        double L_y = 1 / (Z0*(4 + this->cndtvty_x + this->cndtvty_z));
        double ik_pml = this->cndtvty_y*v_series_y - this->V_k_series_n_1[1]*(2-this->cndtvty_y) + I_k_tot_n_1[1]*J_y;
        ik_pml = L_y*ik_pml;

        return ik_pml;
    }
    else if ( k==3 ) //z
    {
        double v_series_z = v_series_k;
        double J_z = Z0*(4 - this->cndtvty_x - this->cndtvty_y);
        double L_z = 1 / (Z0*(4 + this->cndtvty_x + this->cndtvty_y) );
        double ik_pml = this->cndtvty_z*v_series_z - this->V_k_series_n_1[2]*(2-this->cndtvty_z) + I_k_tot_n_1[2]*J_z;
        ik_pml = L_z*ik_pml;

        return ik_pml;
    }

    int j;
    cout<<"failing in i_pml " <<k <<endl;
    cin>>j;
    return 0;

}

double SCN_node::compute_v_dummy_kj(double cndtvty_i,double cndtvty_j,double cndtvty_k, double i_k_pml,double Ik, int axis, double v_dummy_k_n_1)
{
    double L = 1 / ( Z0*(4 + cndtvty_i + cndtvty_j) );
    double K = 1/(2 + cndtvty_k);

    double v_dummy = 2*Z0*K*i_k_pml + cndtvty_j*Z0*K*(4*L*Z0*Ik) + cndtvty_j*Z0*K*i_k_pml  - Z0*K*(2-cndtvty_j)*(4*L*Z0*I_k_n_1[axis-1] + I_k_pml_n_1[axis-1])+
    K*(2-cndtvty_k)*v_dummy_k_n_1;
/*
    if( v_dummy != 0)
    {cout<<" troublesome point" <<endl<<endl;
    cout<<"axis" <<axis<<endl;
    cout<<" ik pml"<< i_k_pml<<endl;
    cout<<" I_k_n_1[axis-1] : " << I_k_n_1[axis-1]<<endl;
    cout<<" 2*Z0*K*i_k_pml :   "<<2*Z0*K*i_k_pml <<endl;
    cout<<" cndtvty_j*Z0*K*(4*L*Z0*Ik) :  " <<cndtvty_j*Z0*K*(4*L*Z0*Ik)<<endl;
    cout<<" cndtvty_j*Z0*K*i_k_pml :   "<<cndtvty_j*Z0*K*i_k_pml <<endl;
    cout<<" - Z0*K*(2-cndtvty_j)*(4*L*Z0*I_k_n_1[axis-1] + I_k_pml_n_1[axis-1]):    " <<-Z0*K*(2-cndtvty_j)*(4*L*Z0*I_k_n_1[axis-1] + I_k_pml_n_1[axis-1])<<endl;


    cout<<" K*(2-cndtvty_k)*v_dummy_k_n_1 :   "<< K*(2-cndtvty_k)*v_dummy_k_n_1<<endl<<endl;

    cin>>K;
    }*/


    return v_dummy;
}


double SCN_node::compute_v_series_k(int k)
{

   if( k==1 )  //x
    {
        return ( this->get_incidence(7)-this->get_incidence(8)+this->get_incidence(4)-this->get_incidence(5) );
    }
    else if( k==2 ) //y
    {
        return ( this->get_incidence(6)-this->get_incidence(10)+this->get_incidence(9)-this->get_incidence(2) ) ;
    }
    else if( k==3 )  //z
    {
        return ( this->get_incidence(1)-this->get_incidence(12)+this->get_incidence(11)-this->get_incidence(3) );
    }

    int j;
    cout<<"failing in v series " <<k <<endl;
    cin>>j;
    return 0;
}

void SCN_node::compute_source_voltage_y(const double v_inc,int source_res,double &Vy,double &Vy_prime)
{
    double Y0 = 1/Z0;
    int Rs = source_res;
    double vi3 = this->get_incidence(3);
    double vi11 = this->get_incidence(11);
    double vi4 = this->get_incidence(4);
    double vi8 = this->get_incidence(8);

    Vy =  ( vi3 + vi11 + vi4 + vi8  ) /2;
    Vy_prime = ( 2*Y0*Rs*( vi3 + vi11 + vi4 + vi8  )  + v_inc )/(4*Rs*Y0 + 1);

}



//.....................................................................................................
/*
mesh_handler_SCN::mesh_handler_SCN(double width, double height, double length, double dl, double tfactor, int n_PMl, int conduct_prof,double Refn_factor)
{

    if ( dl <= 0 || width <=0 || height <=0 || length <=0)
    {
        WG_width =  0;
        WG_height = 0;
        WG_length = 0;
        WG_dl = 0;
        Nx = 0;
        Ny = 0;
        Nz = 0;
        Ntotal = 0;
        check_meshed=0;
        cout<< " troublesome space discretization value chosen" <<endl;
    }
    else
    {
        WG_width = width;
        WG_height = height;
        WG_length = length;
        WG_dl = dl;

        Nx = WG_width / WG_dl;
        Ny = WG_height / WG_dl;
        Nz = WG_length / WG_dl;
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in waveguide

        WG_PML_length = n_PMl * dl;
        Nzz =  2*n_PMl + Nz;
        Ntotal_ = Nzz * Nx *Ny;            // Ntotal + nodes in PML

        int PML_boundary1 = n_PMl;
        int PML_boundary2 = n_PMl+Nz;
        WG_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nt=Nx*Ny;

        for(int z=0; z<Nzz; z++)
            for (int y=0 ; y<Ny ; y++)
            {
                for (int x=0 ; x<Nx ; x++)
                {
                    if((z >= PML_boundary1 )&& (z < PML_boundary2 ) )
                    {
                        WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,false));

                    }
                    else
                    {
                        if ( z < PML_boundary1) WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,true,n_PMl,Refn_factor,conduct_prof,PML_boundary1-z-1));

                        else WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,true,n_PMl,Refn_factor,conduct_prof,z-PML_boundary2));

                    }
                    id = id + 1;
                }
            }
    }


}

double mesh_handler_SCN::Ex_output_at_node(int node_id)
{

    double Gxy(this->WG_nodes[node_id].get_G(1)), Gxz(this->WG_nodes[node_id].get_G(2));

    double Yx(this->WG_nodes[node_id].get_Y(0));

    double Exy =  2*(this->WG_nodes[node_id].get_incidence(1)+this->WG_nodes[node_id].get_incidence(12) + (Yx+2)*this->WG_nodes[node_id].get_incidence(13)-2*this->WG_nodes[node_id].get_incidence(14) );
    Exy = Exy/ ( Yx+Gxy+4);

    double Exz =  2*(this->WG_nodes[node_id].get_incidence(2)+this->WG_nodes[node_id].get_incidence(9) + (Yx+2)*this->WG_nodes[node_id].get_incidence(14)-2*this->WG_nodes[node_id].get_incidence(13) );
    Exz = Exz/ ( Yx+Gxz+4);

    return ((Exy + Exz)/(-1e-3*this->WG_dl));

}

double mesh_handler_SCN::Ey_output_at_node(int node_id)
{

   double Gyx(this->WG_nodes[node_id].get_G(0)), Gyz(this->WG_nodes[node_id].get_G(2));
    //double Gyx(0), Gyz(0);
    double Yy(this->WG_nodes[node_id].get_Y(1));       // note: this is assuming the grid dimensions are all equal to dl

    double Eyx =  2*(this->WG_nodes[node_id].get_incidence(3)+this->WG_nodes[node_id].get_incidence(11) + (Yy+2)*this->WG_nodes[node_id].get_incidence(15)-2*this->WG_nodes[node_id].get_incidence(16) );
    Eyx = Eyx/ ( Yy+Gyx+4);

    double Eyz =  2*(this->WG_nodes[node_id].get_incidence(4)+ this->WG_nodes[node_id].get_incidence(8) + (Yy+2)*this->WG_nodes[node_id].get_incidence(16)-2*this->WG_nodes[node_id].get_incidence(15) );
    Eyz = Eyz/ ( Yy+Gyz+4);

    return ((Eyx + Eyz)/(-1e-3*this->WG_dl));
    //  return (-0.5*(this->WG_nodes[node_id].get_incidence(3) + this->WG_nodes[node_id].get_incidence(11) + this->WG_nodes[node_id].get_incidence(4) +this->WG_nodes[node_id].get_incidence(8))/(1e-3*this->WG_dl) );
}



double mesh_handler_SCN::Ez_output_at_node(int node_id)
{
    double Gzy(this->WG_nodes[node_id].get_G(1)), Gzx(this->WG_nodes[node_id].get_G(0));

    double Yz(this->WG_nodes[node_id].get_Y(2));

    double Ezx =  2*(this->WG_nodes[node_id].get_incidence(6)+ this->WG_nodes[node_id].get_incidence(10) + (Yz+2)*this->WG_nodes[node_id].get_incidence(17)-2*this->WG_nodes[node_id].get_incidence(18) );
    Ezx = Ezx/ ( Yz + Gzx+4);

    double Ezy =  2*(this->WG_nodes[node_id].get_incidence(5)+ this->WG_nodes[node_id].get_incidence(7) + (Yz+2)*this->WG_nodes[node_id].get_incidence(18)-2*this->WG_nodes[node_id].get_incidence(17) );
    Ezy = Ezy/ ( Yz + Gzy+4);

    return ((Ezx + Ezy)/(-1e-3*this->WG_dl));

}

double mesh_handler_SCN::Hx_output_at_node(int node_id)
{
    double Rxy(this->WG_nodes[node_id].get_R(1)), Rxz(this->WG_nodes[node_id].get_R(2));

    double Zx(this->WG_nodes[node_id].get_Z(0));

    double Hxy =  2*(this->WG_nodes[node_id].get_incidence(5)- this->WG_nodes[node_id].get_incidence(7) + (1+2/Zx)*this->WG_nodes[node_id].get_incidence(19)-2*this->WG_nodes[node_id].get_incidence(20)/Zx );
    Hxy = Hxy/ ( Zx+Rxy+4);

    double Hxz =  2*(this->WG_nodes[node_id].get_incidence(8)- this->WG_nodes[node_id].get_incidence(4) + (1+2/Zx)*this->WG_nodes[node_id].get_incidence(20)-2*this->WG_nodes[node_id].get_incidence(19)/Zx );
    Hxz = Hxz/ ( Zx+Rxz+4);

    return (-(Hxy + Hxz)/(Z0*1e-3*this->WG_dl));

}

double mesh_handler_SCN::Hy_output_at_node(int node_id)
{
    double Ryx(this->WG_nodes[node_id].get_R(0)), Ryz(this->WG_nodes[node_id].get_R(2));

    double Zy(this->WG_nodes[node_id].get_Z(1));

    double Hyx =  2*(this->WG_nodes[node_id].get_incidence(10)- this->WG_nodes[node_id].get_incidence(6) + (1+2/Zy)*this->WG_nodes[node_id].get_incidence(21)-2*this->WG_nodes[node_id].get_incidence(22)/Zy );
    Hyx = Hyx/ ( Zy+Ryx+4);

    double Hyz =  2*(this->WG_nodes[node_id].get_incidence(2)- this->WG_nodes[node_id].get_incidence(9) + (1+2/Zy)*this->WG_nodes[node_id].get_incidence(22)-2*this->WG_nodes[node_id].get_incidence(21)/Zy);
    Hyz = Hyz/ ( Zy+Ryz+4);

    return (-(Hyx + Hyz)/(Z0*1e-3*this->WG_dl));
}

double mesh_handler_SCN::Hz_output_at_node(int node_id){

    double Rzy(this->WG_nodes[node_id].get_R(1)), Rzx(this->WG_nodes[node_id].get_R(0));

    double Zz(this->WG_nodes[node_id].get_Z(2));

    double Hzx =  2*(this->WG_nodes[node_id].get_incidence(3)- this->WG_nodes[node_id].get_incidence(11) + (1+2/Zz)*this->WG_nodes[node_id].get_incidence(23)-2*this->WG_nodes[node_id].get_incidence(24)/Zz);
    Hzx = Hzx/ ( Zz+Rzx+4);

    double Hzy =  2*(this->WG_nodes[node_id].get_incidence(12)- this->WG_nodes[node_id].get_incidence(1) + (1+2/Zz)*this->WG_nodes[node_id].get_incidence(24)-2*this->WG_nodes[node_id].get_incidence(23)/Zz);
    Hzy = Hzy/ ( Zz+Rzy+4);

    return (-(Hzx + Hzy)/(Z0*1e-3*this->WG_dl));
    //return (0.5*(this->WG_nodes[node_id].get_incidence(1) - this->WG_nodes[node_id].get_incidence(3) + this->WG_nodes[node_id].get_incidence(11) -this->WG_nodes[node_id].get_incidence(12))/(Z0*1e-3*this->WG_dl) );

   }

void mesh_handler_SCN::set_WG_neighbour()
{
    int id = 0;
    int Nt = this->Nx*this->Ny;
    for(int z=0; z<this->Nzz; z++)
    {
        for (int y=0 ; y<this->Ny ; y++)
        {
            for (int x=0 ; x<this->Nx ; x++)
            {
                if (x!=0)
                {
                    WG_node_neighbours[id].push_back(&WG_nodes[id-1]);
                    WG_node_neighbours[id-1].push_back(&WG_nodes[id]);
                }
                if (y!=0)
                {
                    WG_node_neighbours[id].push_back(&WG_nodes[id-this->Nx]);
                    WG_node_neighbours[id-this->Nx].push_back(&WG_nodes[id]);
                }
                if(z!=0)
                {
                    WG_node_neighbours[id].push_back(&WG_nodes[id-Nt]);
                    WG_node_neighbours[id-Nt].push_back(&WG_nodes[id]);
                }
                id = id+1;
            }
        }
    }
}

void mesh_handler_SCN::print_SCN_neighbour(int node_id)
{
    int countt = WG_node_neighbours[node_id].size();
    for ( int i = 0; i<countt; i++)
        cout<<*(WG_node_neighbours[node_id][i]);
}

void mesh_handler_SCN::print_WG_nodes()
{
    for( int i=0; i< this->WG_nodes.size(); i++)
    {
        cout<<this->WG_nodes[i]<<endl;
    }
}
void mesh_handler_SCN::print_WG_node(int _node)

{
    cout<<" Node : " <<_node << endl;
    for(int i=0; i<12; i++)
    {
        cout<<this->WG_nodes[_node].get_incidence(i+1)<<endl;
    }

}

void mesh_handler_SCN::breakpoint_(int i)
{
    cout<<" breakpoint " << i <<endl;
    cin>>i;
}

void mesh_handler_SCN::print_this_neighbour(int excited_node)
{
    int ii = this->WG_node_neighbours[excited_node].size();

    cout<<" Node : " << excited_node << endl;
    for( int i=1; i<=12; i++)  cout<<this->WG_nodes[excited_node].get_reflected(i)<<endl;
    for (int j =0 ; j<ii ; j++)
    {
        cout<<" Neighbour "<< this->WG_node_neighbours[excited_node][j]->get_iD() <<" : "<<endl;
        cout<<"1: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(1)<<endl;
        cout<<"2: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(2)<<endl;
        cout<<"3: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(3)<<endl;
        cout<<"4: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(4)<<endl;
        cout<<"5: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(5)<<endl;
        cout<<"6: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(6)<<endl;
        cout<<"7: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(7)<<endl;
        cout<<"8: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(8)<<endl;
        cout<<"9: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(9)<<endl;
        cout<<"10: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(10)<<endl;
        cout<<"11: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(11)<<endl;
        cout<<"12: "<<this->WG_node_neighbours[excited_node][j]->get_incidence(12)<<" ? "<<endl;

        //cout<<this->WG_node_neighbours[excited_node][j].get_reflected()<<endl;
    }
}

void mesh_handler_SCN::TLM_simulation1(int dt_total,double tfactor, int output_node,int excited_node,double Zz1, double Zz2, double Zx1, double Zx2,double Zy1, double Zy2,double z,double f0, double sigma)
{

   // if ( ( excited_node == -1 )&&(z == -1) ) cin>> z;

    int nny = this->Ny;
    int nnx = this->Nx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = 0.5*(nnzz - nnz);
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;
    double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
    double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
    double reflct_x = -1;  // (Zx1 - Z0) /( Zx1 + Z0);
    double reflct_y = -1;  //(Zy1 - Z0) /( Zy1 + Z0);

    double test=0;

// Guassian Pulse Parameters
    double f = 0;
    double ff = 0;
    double te = 0;
    double t = tfactor*(1e-3)*dl/(c*2);
    double excitation_length = 0.2*dt_total;
    double variance= pow(sigma,2);

// Compute the start node in desired excitation plane
    int start_node= 0, start_nodee = 0;
    if( z != -1 ) start_node = return_coordinates (w,h,l,n_PML,dl,0,0,z);

    int xy = return_coordinates (w,h,l,n_PML,dl,0,0,dl);                //
    //else excited_node = nxy*n_PML  + excited_node;

    start_nodee = start_node;

// Time loop
    for( int dt=0; dt<dt_total; dt++)
    {

        int checks = 0;

        if((dt<300))   // EXCITE WHOLE PLANE LOCATED AT Z
        {
            //EXCITATION

            if ( z==-1)
            {
                f =(dt*t-(t*20))*sigma* exp( -pow( sigma*(dt*t-(t*100) ),2))*sin(2*pi*f0*dt*t);;

                f = f*(-1e-3*dl);

                te = WG_nodes[excited_node].get_incidence(3);
                ff=f+te;
                WG_nodes[excited_node].set_incidence(3,ff);

                te = WG_nodes[excited_node].get_incidence(11);
                ff=f+te;
                WG_nodes[excited_node].set_incidence(11,ff);

                te = WG_nodes[excited_node].get_incidence(4);
                ff=f+te;
                WG_nodes[excited_node].set_incidence(4,ff);

                te = WG_nodes[excited_node].get_incidence(8);
                ff=f+te;
                WG_nodes[excited_node].set_incidence(8,ff);

                te = WG_nodes[excited_node].get_incidence(15);
                ff=f+te;
                WG_nodes[excited_node].set_incidence(15,ff);

                te = WG_nodes[excited_node].get_incidence(16);
                ff=f+te;
                WG_nodes[excited_node].set_incidence(16,ff);


            }
            else if(z==-2)
            {
                for (int i = xy; i<(xy+nxy);i++ ){

                    f = sin((i)*pi/(Nx-1))*(dt*t-(t*20))*sigma* exp( -pow( sigma*(dt*t-(t*100) ),2))*sin(2*pi*f0*dt*t);

                    f = f*(1e-3*dl);

                    for(int j = 0; j<Ny; j++){

                        te = WG_nodes[i+j*nnx].get_incidence(3);
                        ff=f+te;
                        WG_nodes[i+j*nnx].set_incidence(3,ff);

                        te = WG_nodes[i+j*nnx].get_incidence(11);
                        ff=f+te;
                        WG_nodes[i+j*nnx].set_incidence(11,ff);

                        te = WG_nodes[i+j*nnx].get_incidence(4);
                        ff=f+te;
                        WG_nodes[i+j*nnx].set_incidence(4,ff);

                        te = WG_nodes[i+j*nnx].get_incidence(8);
                        ff=f+te;
                        WG_nodes[i+j*nnx].set_incidence(8,ff);

                        te = WG_nodes[i+j*nnx].get_incidence(15);
                        ff=f+te;
                        WG_nodes[i+j*nnx].set_incidence(15,ff);

                        te = WG_nodes[i+j*nnx].get_incidence(16);
                        ff=f+te;
                        WG_nodes[i+j*nnx].set_incidence(16,ff);

                    }
                }
            }
            else if( z >=0 )
            {
               f =(dt*t-(t*20))*sigma* exp( -pow( sigma*(dt*t-(t*100) ),2))*sin(2*pi*f0*dt*t);;
                //f = 1;
                f = f*(double(-1e-3)*dl);

                for ( int j = 0; j < nxy ; j++)
                {
                    //EXCITE EY field component
                    te = WG_nodes[start_node].get_incidence(3);
                    ff=f+te;
                    WG_nodes[start_node].set_incidence(3,ff);

                    te = WG_nodes[start_node].get_incidence(11);
                    ff=f+te;
                    WG_nodes[start_node].set_incidence(11,ff);

                    te = WG_nodes[start_node].get_incidence(4);
                    ff=f+te;
                    WG_nodes[start_node].set_incidence(4,ff);

                    te = WG_nodes[start_node].get_incidence(8);
                    ff=f+te;
                    WG_nodes[start_node].set_incidence(8,ff);

                    te = WG_nodes[start_node].get_incidence(15);
                    ff=f+te;
                    WG_nodes[start_node].set_incidence(15,ff);

                    te = WG_nodes[start_node].get_incidence(16);
                    ff=f+te;
                    WG_nodes[start_node].set_incidence(16,ff);

                    start_node = start_node + 1;

                }
                start_node = start_nodee;
            }

        }


        //OUTPUT FIELDS
        WG_Ex.push_back(Ex_output_at_node(output_node));

        WG_Ey.push_back(Ey_output_at_node(output_node));

        WG_Ez.push_back(Ez_output_at_node(output_node));

        WG_Hx.push_back(Hx_output_at_node(output_node));

        WG_Hy.push_back(Hy_output_at_node(output_node));

        WG_Hz.push_back(Hz_output_at_node(output_node));

        //SCATTERING PROCESS
        for (int i=0; i<this->Ntotal_; i++)
        {
            WG_nodes[i].scattering_PML() ;
        }



        //CONNECTION PROCESS
        int ii=0;
        int xx = 0;
        int yy = 0;
        int zz = 0;

        for(int z=0; z<nnzz; z++)
        {
            for (int y=0 ; y<this->Ny ; y++)
            {
                for (int x=0 ; x<this->Nx ; x++)
                {
                    if ( (x==0)||(x==nnx-1) )
                    {
                       if ((WG_nodes[ii].is_PML_SCN()==0)){

                            if (x==0)//left wall rule.
                            {
                                WG_nodes[ii].set_incidence(6, reflct_x*(WG_nodes[ii].get_reflected(6)));
                                WG_nodes[ii].set_incidence(3, reflct_x*(WG_nodes[ii].get_reflected(3)));
                                WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                                WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));

                            }
                            else    //right wall rule
                            {
                                WG_nodes[ii].set_incidence(11,reflct_x*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,reflct_x*(WG_nodes[ii].get_reflected(10)));
                                WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                                WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));

                            }
                        }

                          if ((WG_nodes[ii].is_PML_SCN()==1)){

                            if (x==0)//left wall rule.
                            {
                                WG_nodes[ii].set_incidence(6, reflct_x*(WG_nodes[ii].get_reflected(6)));
                                WG_nodes[ii].set_incidence(3, reflct_x*(WG_nodes[ii].get_reflected(3)));
                                WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                                WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));

                            }
                            else    //right wall rule
                            {
                                WG_nodes[ii].set_incidence(11,reflct_x*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,reflct_x*(WG_nodes[ii].get_reflected(10)));
                                WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                                WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                            }
                        }
                    }
                    else
                    {
                        // reflect into adjacent nodes
                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                    }

                    if (( y==0) || (y==nny-1))
                    {
                        if ((WG_nodes[ii].is_PML_SCN()==0)){
                            if (y==0)   //bottom wall rule
                            {
                                WG_nodes[ii].set_incidence(1, reflct_y*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, reflct_y*WG_nodes[ii].get_reflected(5));
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));

                            }
                            else        //top wall rule
                            {
                                WG_nodes[ii].set_incidence(12, reflct_y*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7,  reflct_y*WG_nodes[ii].get_reflected(7));
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                            }
                        }

                        if ((WG_nodes[ii].is_PML_SCN()==1)){
                         if (y==0)   //bottom wall rule
                            {
                                WG_nodes[ii].set_incidence(1, reflct_y*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, reflct_y*WG_nodes[ii].get_reflected(5));
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                            }
                            else        //top wall rule
                            {
                                WG_nodes[ii].set_incidence(12,reflct_y*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, reflct_y*WG_nodes[ii].get_reflected(7));
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                            }
                        }
                    }

                    else       //reflect into adjacent nodes
                    {
                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                    }


                    if ( (z==0) || (z==nnzz-1))
                    {
                        if (z==0)   // front face rule
                        {
                            WG_nodes[ii].set_incidence(8, reflct_z1*WG_nodes[ii].get_reflected(8));
                            WG_nodes[ii].set_incidence(9, reflct_z1*WG_nodes[ii].get_reflected(9));
                            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));

                        }
                        else        //back face rule
                        {
                            WG_nodes[ii].set_incidence(4,reflct_z2*WG_nodes[ii].get_reflected(4));
                            WG_nodes[ii].set_incidence(2,reflct_z2*WG_nodes[ii].get_reflected(2));
                            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                        }
                    }
                    else        //normal reflection rules
                    {
                        // reflect into adjacent nodes
                        WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                        WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                        WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                        WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                    }
                    ii=ii+1;
                }
            }

        }
        ii=0;
    }

}
void mesh_handler_SCN:: write_output_file (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;
    ofstream out_fileEX;
    ofstream out_fileEZ;

    ofstream out_fileHY;
    ofstream out_fileHX;
    ofstream out_fileHZ;


    string nodefileEX = "PEC_EX.txt";
    string nodefileEZ = "PEC_EZ.txt";

    string nodefileHX = "PEC_HX9.txt";
    string nodefileHY = "PEC_HY.txt";
    string nodefileHZ = "PEC_HZ.txt";

    out_fileEX.open((nodefileEX+"").c_str());
    out_fileEZ.open((nodefileEZ+"").c_str());
    out_fileEY.open((nodefileEY+"").c_str());

    out_fileHX.open((nodefileHX+"").c_str());
    out_fileHY.open((nodefileHY+"").c_str());
    out_fileHZ.open((nodefileHZ+"").c_str());


    int end_ = this->WG_Ex.size();
    int end__= this->WG_Ey.size();

    cout<< " Writing to file...." << endl;

    for(int ii = 0 ; ii< end_; ii++)
    {
     out_fileEY<<setprecision(numeric_limits<double>::digits10+1) <<WG_Ey[ii]<<endl;
     out_fileHX<<setprecision(numeric_limits<double>::digits10+1)<<WG_Hx[ii]<<endl;
     out_fileEX<<setprecision(numeric_limits<double>::digits10+1)<<WG_Ex[ii]<<endl;
     out_fileHY<<setprecision(numeric_limits<double>::digits10+1)<<WG_Hy[ii]<<endl;
     out_fileEZ<<setprecision(numeric_limits<double>::digits10+1)<<WG_Ez[ii]<<endl;
     out_fileHZ<<setprecision(numeric_limits<double>::digits10+1)<<WG_Hz[ii]<<endl;
    }
}


void mesh_handler_SCN::print_Ex_output()
{
    cout<< " PRINTING Ey......." <<endl;
    for(int i= 0; i<WG_Ex.size(); i++)
        cout<<i << ": "<< setprecision(15)<<WG_Ey[i] <<endl;
}

void mesh_handler_SCN::analytical_fc()
{
  cout<< " The analytical resonant frequencies:  " <<endl;  // a = width ; b = height ; c = length ;

    int W = WG_width/WG_dl;
    W = W*WG_dl;
    int L = WG_length/WG_dl;
    L= L*WG_dl;

    cout<< W<<endl;
    cout << L<<endl;

            cout<<"TE101 : ";
            cout<< sqrt((pi/WG_width)*(pi/WG_width) + (pi/WG_length)*(pi/WG_length))*300000000000/(2*pi)<<endl;

            cout<<"TE102 : ";
            cout<< sqrt((pi/WG_width)*(pi/WG_width) + (2*pi/WG_length)*(2*pi/WG_length))*300000000000/(2*pi)<<endl;

            cout<<"TE103 : ";
            cout<< sqrt((pi/WG_width)*(pi/WG_width) + (3*pi/WG_length)*(3*pi/WG_length))*300000000000/(2*pi)<<endl;

            cout<<"TE104 : ";
            cout<< sqrt((pi/WG_width)*(pi/WG_width) + (4*pi/WG_length)*(4*pi/WG_length))*300000000000/(2*pi)<<endl;

            cout<<"TE105 : ";
            cout<< sqrt((pi/WG_width)*(pi/WG_width) + (5*pi/WG_length)*(5*pi/WG_length))*300000000000/(2*pi)<<endl;


}

int return_coordinates(double width, double height, double length,int npml, double dl, double x, double y, double z)
{
    int h = height/dl;
    int l = length/dl;
    int w = width/ dl;
    int s = h*w*npml;

    int c = x/dl;
    int d = y/dl;
    int e = z/dl;

    if(l == e) e = e-1;
    if (h == d) d = d-1;
    if (w == c) c= c - 1;

    int node_id =  c + d*w + e*(h*w);
    //cout<<node_id;

    return node_id+s;

}
*/
