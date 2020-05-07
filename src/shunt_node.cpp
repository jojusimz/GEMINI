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


const double pi = 3.14159265;
const double  u0 = 12566370614e-16;
const double e0 = 88541878176e-22;
const double err= 1;
const double c = 299792458.;
const double Z0 = sqrt(u0/e0);//
const double Z0_shunt = Z0*sqrt(2);
double er(1),ur(1);

//CONSTRUCTOR
shunt_node::shunt_node (bool freespace_, float err,float sigma_e, int id_, double dt_,int xx, int yy, double dl_,
                        bool PML, int Lx, int Ly, double sigmax,double sigmay,int cndct_prof)
{
    shunt_port.assign( 5,make_pair(0,0));
    shunt_coords.push_back(xx);
    shunt_coords.push_back(yy);

    id = id_;
    dt = dt_;
    dl = dl_;
    v_total = 0;
    i_pml_sum = 0;
    isPML = PML;
    cnductvty_x = 0 ;//*PML_x;
    cnductvty_y = 0;//*PML_y;
    PEC =false;
    Ge = 1e-3*dl*sigma_e;
    //cout<<err<<endl;
    er_ = err;
    Yst = (2*sqrt(2)*(er_-1))/Z0;
    freespace = freespace_;
    cnductvty_x = 0;
    cnductvty_y = 0;
    //er= er_*e0;

    if(PML)
    {
        double sigma_ex(0),sigma_ey(0);
        if (Lx != -1)
        {
            sigma_ex = sigmax*(pow((Lx+1),cndct_prof+1)- pow(Lx,cndct_prof+1));
            cnductvty_x = sigma_ex/e0;          //*PML_x;
            //cnductvty_x = 0;
        }
        if (Ly != -1)
        {
            sigma_ey = sigmay*(pow((Ly+1),cndct_prof+1)- pow(Ly,cndct_prof+1));
            cnductvty_y = sigma_ey/e0;  //*PML_x;
            //cnductvty_y =0;
        }
    }

    // v and i per timestep
    v_i_ndt__1.assign(5, make_pair(0,0) );             // stored voltage and current from previous time step
    v_ndt__2.assign(5,0 );                               // stored voltage from previous previous time step
    v_i_total_ndt__1 = make_pair(0,0);                    // stored total voltage and current from previous time step.
    v_total_ndt__2 = 0;

}

shunt_node::shunt_node(const shunt_node& copy_shunt_node)                      // copy constructor
{
    shunt_port = copy_shunt_node.shunt_port;
    shunt_coords = copy_shunt_node.shunt_coords;
    id = copy_shunt_node.id;
    dt = copy_shunt_node.dt;
    dl = copy_shunt_node.dl;
    v_total =copy_shunt_node.v_total;
    PEC=copy_shunt_node.PEC;
    //STUBS
    freespace = copy_shunt_node.freespace;
    Yst = copy_shunt_node.Yst;
    er_ = copy_shunt_node.er_;
    Ge  = copy_shunt_node.Ge;
        //PML
    isPML = copy_shunt_node.isPML;                                             //
    PML_Lx = copy_shunt_node.PML_Lx;
    PML_Ly = copy_shunt_node.PML_Ly;
    cnductvty_x = copy_shunt_node.cnductvty_x;                                        //conductivity
    cnductvty_y = copy_shunt_node.cnductvty_y;
    i_pml_sum  = copy_shunt_node.i_pml_sum;
    v_i_ndt__1 = copy_shunt_node.v_i_ndt__1;             // stored voltage and current from previous time step
    v_ndt__2 = copy_shunt_node.v_ndt__2;
    v_i_total_ndt__1 = copy_shunt_node.v_i_total_ndt__1;              // stored total voltage and current from previous time step.
    v_total_ndt__2 = copy_shunt_node.v_total_ndt__2;
}

shunt_node & shunt_node::operator=(const shunt_node& c_shunt_node)                      // assignment operator
{
   if (this!= &c_shunt_node)
    {
    shunt_port = c_shunt_node.shunt_port;
    shunt_coords = c_shunt_node.shunt_coords;
    id = c_shunt_node.id;
    dt = c_shunt_node.dt;
    dl = c_shunt_node.dl;
    v_total =c_shunt_node.v_total;
    PEC=c_shunt_node.PEC;
     //STUBS
    freespace = c_shunt_node.freespace;
    Yst = c_shunt_node.Yst;
    er_ = c_shunt_node.er_;
    Ge  = c_shunt_node.Ge;
        //PML
    isPML = c_shunt_node.isPML;                                             //
    PML_Lx = c_shunt_node.PML_Lx;
    PML_Ly = c_shunt_node.PML_Ly;
    cnductvty_x = c_shunt_node.cnductvty_x;                                        //conductivity
    cnductvty_y = c_shunt_node.cnductvty_y;
    i_pml_sum  = c_shunt_node.i_pml_sum;
    v_i_ndt__1 = c_shunt_node.v_i_ndt__1;             // stored voltage and current from previous time step
    v_ndt__2 = c_shunt_node.v_ndt__2;
    v_i_total_ndt__1 = c_shunt_node.v_i_total_ndt__1 ;             // stored total voltage and current from previous time step.
    v_total_ndt__2 = c_shunt_node.v_total_ndt__2 ;             // stored total voltage and current from previous time step.

    }
    return (*this);

}


void shunt_node::shunt_excitation (double V_i)
{
    this->shunt_port[0].first = this->shunt_port[0].first + V_i;
    this->shunt_port[1].first = this->shunt_port[1].first + V_i;
    this->shunt_port[2].first = this->shunt_port[2].first + V_i;
    this->shunt_port[3].first = this->shunt_port[3].first + V_i;
    if( !freespace ) this->shunt_port[4].first = this->shunt_port[4].first + V_i;

}

//OPERATOR OVERLOAD <<
ostream& operator <<(std::ostream& out,shunt_node& node)
{
    out<<" [ "<< node.id<<" | "<< node.shunt_coords[0] <<" "<< node.shunt_coords[1] <<" "<<node.shunt_coords[2]<< " ] "<<endl;
    return out;
}

//GETS
double shunt_node::get_reflected(int port_id)
{
    return this->shunt_port[port_id].second;
}

double shunt_node::get_incidence( int port_id)
{
    return this->shunt_port[port_id].first;
}

int shunt_node::get_iD()
{
    return this->id;
}

double shunt_node::get_dt()
{
    return this->dt;
}

float shunt_node::get_cnductvty_x()
{
    return this->cnductvty_x;
}

float shunt_node::get_cnductvty_y()
{
    return this->cnductvty_y;
}

int shunt_node::get_coord(int i)
{
    return this->shunt_coords[i];
}

double shunt_node::get_Ez()
{  // cout<<"v: "<<v_total<<endl;
    return( ( -this->v_total)/(1e-3*this->dl)   );
}

double shunt_node::get_Hx()
{
    if (this->isPML) return 0;
    else return(( shunt_port[3].first - shunt_port[1].first )/(Z0_shunt*1e-3*this->dl)) ;
}

double shunt_node::get_Hy()
{
    if (this->isPML) return 0;
    else return(( shunt_port[0].first - shunt_port[2].first )/(Z0_shunt*1e-3*this->dl)) ;
}

//SETS
void  shunt_node::set_reflected( int port_id,double Vr)
{
    this->shunt_port[port_id].second = Vr;
}

void  shunt_node::set_incidence(int port_id,double Vi)
{
    this->shunt_port[port_id].first = Vi;
}

//MISC
 bool shunt_node:: is_PML_shunt()
 {
     return this->isPML;
 }
 //TLM PROCESS
void shunt_node::shunt_scatter()
{
    double vi0(0),vi1(0),vi2(0),vi3(0);
    double vr0(0),vr1(0),vr2(0),vr3(0);
    double vr0_(0),vr1_(0),vr2_(0),vr3_(0);

    vi0 = this->shunt_port[0].first;
    vi1 = this->shunt_port[1].first;
    vi2 = this->shunt_port[2].first;
    vi3 = this->shunt_port[3].first;


    if(!this->isPML)
    {
       if(this->freespace )
       {
            vr0 = 0.5*(-vi0+vi1+vi2+vi3);
            vr1 = 0.5*(vi0-vi1+vi2+vi3);
            vr2 = 0.5*(vi0+vi1-vi2+vi3);
            vr3 = 0.5*(vi0+vi1+vi2-vi3);

            this->v_total = 0.5*(vi0 + vi1 + vi2 + vi3);                //useful in observation member function.
        }

       else
        {
                //cout<<"dielectric"<<endl;
            double vr4;
            double vi4 = this->shunt_port[4].first;
            //cout<<Yst<<endl;
            double v_ttl = ( 2*(vi0 + vi1 + vi2 + vi3)/Z0_shunt + 2*vi4*Yst ) / (4/Z0_shunt + Yst + Ge ) ;
            vr0 = v_ttl - vi0;
            vr1 = v_ttl - vi1;
            vr2 = v_ttl - vi2;
            vr3 = v_ttl - vi3;
            vr4 = v_ttl - vi4;


            this->v_total = v_ttl;

            this->shunt_port[4].second = vr4;
        }

        this->shunt_port[0].second = vr0;
        this->shunt_port[1].second = vr1;
        this->shunt_port[2].second = vr2;
        this->shunt_port[3].second = vr3;
    }

    else
    {
        //freespace  and stubs
        if ( freespace )
        {
            double vz_n_1(this->v_i_total_ndt__1.first);
            double iz_n_1(this->v_i_total_ndt__1.second);
            double i_pml_total(0);
            double sigmax = this->cnductvty_x;
            double sigmay = this->cnductvty_y;
            double dt_scale_x = exp(-sigmax*dt);
            double dt_scale_y = exp(-sigmay*dt);
            double v_pml_terms(0);
            int i(0);

            double v_sum = 2*(vi0 + vi1 + vi2 + vi3);

            double F = 1/( 4 + sigmax*dt + sigmay*dt);

            double v_sum_scaled = F*v_sum;

            //usual tlm scattering
            vr0_ = v_sum_scaled - vi0 ;
            vr1_ = v_sum_scaled - vi1 ;
            vr2_ = v_sum_scaled - vi2 ;
            vr3_ = v_sum_scaled - vi3 ;

            this->compute_total_i();                            //updates current and computes ipml
            i_pml_total = this->i_pml_sum;

            v_pml_terms = F*( Z0_shunt*i_pml_total - Z0_shunt*iz_n_1 + (4 - sigmax*dt - sigmay*dt)*vz_n_1 );

            //stores TOTAL VOLTAGE for next time step.
            this->v_i_total_ndt__1.first = v_pml_terms + v_sum_scaled;
            this->v_total = v_pml_terms + v_sum_scaled;                //useful in observation member function.

          //reflected voltage
        this->shunt_port[0].second = dt_scale_y*(vr0_ + v_pml_terms);
        this->shunt_port[1].second = dt_scale_x*(vr1_ + v_pml_terms);
        this->shunt_port[2].second = dt_scale_y*(vr2_ + v_pml_terms);
        this->shunt_port[3].second = dt_scale_x*(vr3_ + v_pml_terms);

        }
        else this->shunt_scatter_with_stubs();

    }

}

void shunt_node::shunt_scatter_with_stubs ()
{
    //constants
    double sigmax = this->cnductvty_x;
    double sigmay = this->cnductvty_y;
    double vr0(0),vr1(0),vr2(0),vr3(0), vr4(0);

    double vz_n_1 = this->v_i_total_ndt__1.first;
    double vz_n_2= this->v_total_ndt__2;
    double vi0 = this->shunt_port[0].first;
    double vi1 = this->shunt_port[1].first;
    double vi2 = this->shunt_port[2].first;
    double vi3 = this->shunt_port[3].first;
    double voz = this->shunt_port[4].first;

    double v_n_1_0 = this->v_i_ndt__1[0].first;
    double v_n_1_1 = this->v_i_ndt__1[1].first;
    double v_n_1_2 = this->v_i_ndt__1[2].first;
    double v_n_1_3 = this->v_i_ndt__1[3].first;
    double v_n_1_4 = this->v_i_ndt__1[4].first;

    double v_n_2_0 = this->v_ndt__2[0];
    double v_n_2_1 = this->v_ndt__2[1];
    double v_n_2_2 = this->v_ndt__2[2];
    double v_n_2_3 = this->v_ndt__2[3];
    double v_n_2_4 = this->v_ndt__2[4];

    double dt_ = (dl*1e-3/(c*sqrt(2)));

    double alpha = 4*(er_-1);
    double G = alpha + Ge*Z0_shunt;
    double Zst = 4 + G;

    double A = alpha*(( sigmax +sigmay )*dt_ + 0.5*sigmax*sigmay*dt_*dt_);
    double A2 = alpha*(-( sigmax +sigmay )*dt_ + 0.5*sigmax*sigmay*dt_*dt_);
    double B = alpha*sigmax*sigmay*dt_*dt_;
    double C = 8 + G*(2-0.5*sigmax*sigmay*dt_*dt_);
    double D = 4 + G -( 1 + G*0.5)*(sigmax + sigmay)*dt_ + G*0.25*sigmax*sigmay*dt_*dt_;
    double F = 1 / ( 4 + G +( 1 + 0.5*G)*(sigmax+sigmay)*dt_ + G*0.25*sigmax*sigmay*dt_*dt_);

        //cmputing the denominator
    double V_part1 = 2*(vi1 + vi3 + vi0 + vi2 + alpha*voz) - 4*(v_n_1_0 + v_n_1_1 + v_n_1_2 + v_n_1_3 + alpha*v_n_1_4)+ 2*(v_n_2_0 + v_n_2_1 + v_n_2_2 + v_n_2_3 + alpha*v_n_2_4);
    double V_part2 = (vi1 + vi3)*sigmay*dt_ + (vi0 + vi2)*sigmax*dt_ - (v_n_2_1 + v_n_2_3)*sigmay*dt_ -(v_n_2_0 + v_n_2_2)*sigmax*dt_ ;
    double V_part3 = A*voz + B*v_n_1_4 + A2*v_n_2_4;
    double V_part4 = C*vz_n_1 - D*vz_n_2;
    double Vz = V_part1 + V_part2 + V_part3 + V_part4 ;
    Vz = Vz*F;

    this->v_total= Vz;

    // compute reflected voltages
    vr0 = Vz - vi0;
    vr1 = Vz - vi1;
    vr2 = Vz - vi2;
    vr3 = Vz - vi3;
    vr4 = Vz - voz;

    //reflected voltage

    double dt_scale_x = exp(-sigmax*dt_);
    double dt_scale_y = exp(-sigmay*dt_);
    //cout<<dt_scale_x<<endl;
    //cout<<dt_scale_y<<endl;
   // cin>>dt_scale_x;

    this->shunt_port[0].second = dt_scale_y*vr0;
    this->shunt_port[1].second = dt_scale_x*vr1;
    this->shunt_port[2].second = dt_scale_y*vr2;
    this->shunt_port[3].second = dt_scale_x*vr3;
    this->shunt_port[4].second = 1*vr4;

    //updated voltages
    this->update_voltages_shunt_node();

}

void shunt_node::shunt_scatter_with_stubs2 ()
{
    //constants
    double sigmax = this->cnductvty_x;
    double sigmay = this->cnductvty_y;
    double vr0(0),vr1(0),vr2(0),vr3(0), vr4(0);

    double vz_n_1 = this->v_i_total_ndt__1.first;
    double vz_n_2= this->v_total_ndt__2;
    double vi0 = this->shunt_port[0].first;
    double vi1 = this->shunt_port[1].first;
    double vi2 = this->shunt_port[2].first;
    double vi3 = this->shunt_port[3].first;
    double voz = this->shunt_port[4].first;

    double v_n_1_0 = this->v_i_ndt__1[0].first;
    double v_n_1_1 = this->v_i_ndt__1[1].first;
    double v_n_1_2 = this->v_i_ndt__1[2].first;
    double v_n_1_3 = this->v_i_ndt__1[3].first;
    double v_n_1_4 = this->v_i_ndt__1[4].first;

    double v_n_2_0 = this->v_ndt__2[0];
    double v_n_2_1 = this->v_ndt__2[1];
    double v_n_2_2 = this->v_ndt__2[2];
    double v_n_2_3 = this->v_ndt__2[3];
    double v_n_2_4 = this->v_ndt__2[4];

    double dt_ = (dl*1e-3/(c*sqrt(2)));
    double Z = 4*er + Ge*Z0_shunt;
    double A = sigmax*sigmay*dt_*dt_;
    double F = 1 + 0.5*( sigmax + sigmay )*dt_ + 0.25*A;

    //if( id == 0){
        /*
      vi0 = 1;
      vi1 = 1;
      vi2 = 1;
      vi3 = 1;
      voz = 1;
      this->shunt_port[0].first = 1;
      this->shunt_port[1].first = 1;
      this->shunt_port[2].first = 1;
      this->shunt_port[3].first = 1;
      this->shunt_port[4].first = 1;
*/
    double Vx = vi1 + vi3 - 2*voz;
    double Vy = vi2 + vi0 - 2*voz;
    double Vx_n_1 = v_n_1_1 + v_n_1_3 - 2*v_n_1_4;
    double Vy_n_1 = v_n_1_2 + v_n_1_0 - 2*v_n_1_4;
    double Vx_n_2 = v_n_2_1 + v_n_2_3 - 2*v_n_2_4;
    double Vy_n_2 = v_n_2_2 + v_n_2_0 - 2*v_n_2_4;


    double V_part1 = Vx*( 2+sigmay*dt_ ) + Vy*( 2+sigmax*dt_ ) + 8*voz*er*F ;
    double V_part2 = -4*(Vx_n_1 + Vy_n_1) - (8*v_n_1_4*er - vz_n_1*Z)*(2 - 0.5*A) ;
    double V_part3 = Vx_n_2*(2-sigmay*dt_) + Vy_n_2*(2-sigmax*dt_) + ( 8*v_n_2_4*er -  Z*vz_n_2 )*(1 - 0.5*(sigmax+sigmay)*dt_ + 0.25*A);
    double Vz = V_part1 + V_part2 + V_part3 ;
    Vz = Vz/(Z*F);

/*
    cout<< "F " << F<<endl;
    cout<< "Z " << Z<<endl;
    cout<< "A " << A<<endl;
    cout<< " Vx "<<Vx<<endl;
    cout<< " Vy "<<Vy<<endl;
    cout<< " V_part1 " << V_part1<<endl;
    cout<< " V_part2 " << V_part2<<endl;
    cout<< " V_part3 " << V_part3<<endl;
    cout<< " Vz " << Vz<<endl;
    cout<< "--->>> " << Vx_n_2*(2-sigmay*dt_) + Vy_n_2*(2-sigmax*dt_) <<endl;
    cout<< "---->>> " << 8*v_n_2_4*er -  Z*vz_n_2 <<endl;
    cout<< " ---->> " << v_n_2_4<<endl;
    cout<< " Vz_n_1 "<< vz_n_1<<endl;
    cout<< " Vz_n_2 "<< vz_n_2<<endl;
    cout<< Vx_n_2  << " " << Vy_n_2 <<endl;
    cout<<endl;
    cout<<endl;
    cin>>Z;
*/
    this->v_total= Vz;

    // compute reflected voltages
    vr0 = Vz - vi0;
    vr1 = Vz - vi1;
    vr2 = Vz - vi2;
    vr3 = Vz - vi3;
    vr4 = Vz - voz;

    //reflected voltage

    double dt_scale_x = exp(-sigmax*dt_);
    double dt_scale_y = exp(-sigmay*dt_);
    //cout<<dt_scale_x<<endl;
    //cout<<dt_scale_y<<endl;
   // cin>>dt_scale_x;

    this->shunt_port[0].second = dt_scale_y*vr0;
    this->shunt_port[1].second = dt_scale_x*vr1;
    this->shunt_port[2].second = dt_scale_y*vr2;
    this->shunt_port[3].second = dt_scale_x*vr3;
    this->shunt_port[4].second = 1*vr4;

    //updated voltages
    this->update_voltages_shunt_node();



}

void shunt_node::shunt_scatter_with_stubs3 ()
{
    //constants
    double sigmax = this->cnductvty_x;
    double sigmay = this->cnductvty_y;
    double vr0(0),vr1(0),vr2(0),vr3(0), vr4(0);

    double vz_n_1 = this->v_i_total_ndt__1.first;
    double vz_n_2= this->v_total_ndt__2;
    double vi0 = this->shunt_port[0].first;
    double vi1 = this->shunt_port[1].first;
    double vi2 = this->shunt_port[2].first;
    double vi3 = this->shunt_port[3].first;
    double voz = this->shunt_port[4].first;

    double v_n_1_0 = this->v_i_ndt__1[0].first;
    double v_n_1_1 = this->v_i_ndt__1[1].first;
    double v_n_1_2 = this->v_i_ndt__1[2].first;
    double v_n_1_3 = this->v_i_ndt__1[3].first;
    double v_n_1_4 = this->v_i_ndt__1[4].first;

    double v_n_2_0 = this->v_ndt__2[0];
    double v_n_2_1 = this->v_ndt__2[1];
    double v_n_2_2 = this->v_ndt__2[2];
    double v_n_2_3 = this->v_ndt__2[3];
    double v_n_2_4 = this->v_ndt__2[4];

    double dt_ = (dl*1e-3/(c*sqrt(2)));
    //double Z = 4*er + Ge*Z0_shunt;
    //double A = sigmax*sigmay*dt_*dt_;
    double F = 16 + 16*sigmay*dt_ + 2*sigmay*sigmay*dt_*dt_;
/*
    if( id == 0){

      vi0 = 1;
      vi1 = 1;
      vi2 = 1;
      vi3 = 1;
      voz = 1;
      this->shunt_port[0].first = 1;
      this->shunt_port[1].first = 1;
      this->shunt_port[2].first = 1;
      this->shunt_port[3].first = 1;
      this->shunt_port[4].first = 1;
*/
    double Vx = vi1 + vi3;
    double Vy = vi2 + vi0 + 2*voz;
    double Vx_n_1 = v_n_1_1 + v_n_1_3 ;
    double Vy_n_1 = v_n_1_2 + v_n_1_0 + 2*v_n_1_4;
    double Vx_n_2 = v_n_2_1 + v_n_2_3 ;
    double Vy_n_2 = v_n_2_2 + v_n_2_0 + 2*v_n_2_4;


    double V_part1 = Vx*( 8 + 8*sigmay*dt_ + 2*sigmay*sigmay*dt_*dt_ ) + Vy*(8 + 4*sigmay*dt_) -16*voz;
    double V_part2 = Vx_n_1*(-16 + 4*sigmay*sigmay*dt_*dt_) - 16*Vy_n_1 + 32*v_n_1_4 - vz_n_1*(-32 + 4*sigmay*sigmay*dt_*dt_);
    double V_part3 = Vx_n_2*( 8 - 8*sigmay*dt_ + 2*sigmay*sigmay*dt_*dt_ ) + Vy_n_2*(8 - 4*sigmay*dt_) - 16*v_n_2_4 -  vz_n_2*(16 - 16*sigmay*dt_ + 2*sigmay*sigmay*dt_*dt_);
    double Vz = V_part1 + V_part2 + V_part3 ;
    Vz = Vz/F;

/*
    cout<< "F " << F<<endl;
    cout<< " Vx "<<Vx<<endl;
    cout<< " Vy "<<Vy<<endl;
    cout<< " V_part1 " << V_part1<<endl;
    cout<< " V_part2 " << V_part2<<endl;
    cout<< " V_part3 " << V_part3<<endl;
    cout<< " Vz " << Vz<<endl;
    cout<< "--->>> " <<endl;
    cout<< "---->>> " <<endl;
    cout<< " ---->> " << v_n_2_4<<endl;
    cout<< " Vz_n_1 "<< vz_n_1<<endl;
    cout<< " Vz_n_2 "<< vz_n_2<<endl;
    cout<< Vx_n_2  << " " << Vy_n_2 <<endl;
    cout<<endl;
    cout<<endl;
    cin>>F;
*/
    this->v_total= Vz;

    // compute reflected voltages
    vr0 = Vz - vi0;
    vr1 = Vz - vi1;
    vr2 = Vz - vi2;
    vr3 = Vz - vi3;
    vr4 = Vz - voz;

    //reflected voltage

    double dt_scale_x = exp(-sigmax*dt_);
    double dt_scale_y = exp(-sigmay*dt_);
    //cout<<dt_scale_x<<endl;
    //cout<<dt_scale_y<<endl;
   // cin>>dt_scale_x;

    this->shunt_port[0].second = dt_scale_y*vr0;
    this->shunt_port[1].second = dt_scale_x*vr1;
    this->shunt_port[2].second = dt_scale_y*vr2;
    this->shunt_port[3].second = dt_scale_x*vr3;
    this->shunt_port[4].second = 1*vr4;

    //updated voltages
    this->update_voltages_shunt_node();

}


 void shunt_node::update_voltages_shunt_node()
 {

   this->v_ndt__2[0] = this->v_i_ndt__1[0].first;
   this->v_ndt__2[1] = this->v_i_ndt__1[1].first;
   this->v_ndt__2[2] = this->v_i_ndt__1[2].first;
   this->v_ndt__2[3] = this->v_i_ndt__1[3].first;
   this->v_ndt__2[4] = this->v_i_ndt__1[4].first;

    this->v_i_ndt__1[0].first = this->shunt_port[0].first;
    this->v_i_ndt__1[1].first = this->shunt_port[1].first;
    this->v_i_ndt__1[2].first = this->shunt_port[2].first;
    this->v_i_ndt__1[3].first = this->shunt_port[3].first;
    this->v_i_ndt__1[4].first = this->shunt_port[4].first;

    this->v_total_ndt__2 = this->v_i_total_ndt__1.first;
    this->v_i_total_ndt__1.first = this->v_total;

 }

 void shunt_node::compute_total_i()
 {

    double i_pml0(0), i_pml1(0), i_pml2(0), i_pml3(0);              // pml current from current time step
    double v_n_1_0(0), v_n_1_1(0), v_n_1_2(0), v_n_1_3(0);          // incident voltage in node from previous time step
    double i_n_1_0(0), i_n_1_1(0), i_n_1_2(0), i_n_1_3(0);          // total current from previous time step
    double vi0(0),vi1(0),vi2(0),vi3(0);
    double sigmax( this->cnductvty_x), sigmay (this->cnductvty_y);

    vi0 = this->shunt_port[0].first;
    vi1 = this->shunt_port[1].first;
    vi2 = this->shunt_port[2].first;
    vi3 = this->shunt_port[3].first;


// the voltage and current stored in the previous time step
    v_n_1_0 = this->v_i_ndt__1[0].first;
    v_n_1_1 = this->v_i_ndt__1[1].first;
    v_n_1_2 = this->v_i_ndt__1[2].first;
    v_n_1_3 = this->v_i_ndt__1[3].first;

    i_n_1_0 = this->v_i_ndt__1[0].second;
    i_n_1_1 = this->v_i_ndt__1[1].second;
    i_n_1_2 = this->v_i_ndt__1[2].second;
    i_n_1_3 = this->v_i_ndt__1[3].second;


    i_pml0 = ( -2*v_n_1_0 + sigmax*dt*vi0 + sigmax*dt*v_n_1_0 ) /Z0_shunt + i_n_1_0;
    i_pml1 = ( -2*v_n_1_1 + sigmay*dt*vi1 + sigmay*dt*v_n_1_1 ) /Z0_shunt + i_n_1_1;
    i_pml2 = ( -2*v_n_1_2 + sigmax*dt*vi2 + sigmax*dt*v_n_1_2 ) /Z0_shunt + i_n_1_2;
    i_pml3 = ( -2*v_n_1_3 + sigmay*dt*vi3 + sigmay*dt*v_n_1_3 ) /Z0_shunt + i_n_1_3;

    // total ipml in shunt node
    this->i_pml_sum = i_pml0 + i_pml1 + i_pml2 + i_pml3;
    double v_sum = 2*(vi0 + vi1 + vi2 + vi3);

    //total current in each line in present time step is stored for next time step
    this->v_i_ndt__1[0].second = i_pml0 + 2*vi0/Z0_shunt;
    this->v_i_ndt__1[1].second = i_pml1 + 2*vi1/Z0_shunt;
    this->v_i_ndt__1[2].second = i_pml2 + 2*vi2/Z0_shunt;
    this->v_i_ndt__1[3].second = i_pml3 + 2*vi3/Z0_shunt;

    // voltage in present time step stored for next time step.
    this->v_i_ndt__1[0].first = vi0;
    this->v_i_ndt__1[1].first = vi1;
    this->v_i_ndt__1[2].first = vi2;
    this->v_i_ndt__1[3].first = vi3;

    //total current in present node and updated for next time step
    this->v_i_total_ndt__1.second = this->i_pml_sum + v_sum/Z0_shunt;

 }

void shunt_node::compute_source_voltages(float ez_inc,int source_res,float &Vz,float &Vz_prime)
{

    double Y0 = 1/Z0_shunt;
    int Rs = source_res;
    double vi0 = this->shunt_port[0].first;
    double vi1 = this->shunt_port[1].first;
    double vi2 = this->shunt_port[2].first;
    double vi3 = this->shunt_port[3].first;
    double voz = this->shunt_port[4].first;

    Vz = (  2*Y0*( vi0 + vi1 + vi2 + vi3  ) + 2*Yst*voz)/(4*Y0 + Yst);

    Vz_prime = (  2*Y0*Rs*( vi0 + vi1 + vi2 + vi3  ) + 2*Yst*Rs*voz + ez_inc   )/(4*Rs*Y0 + Rs*Yst + 1);
}

shunt_node::~shunt_node()
{
    //dtor
}
