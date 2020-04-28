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
#include <omp.h>
#include <time.h>

#include "SCN_node.h"
#include "WG_16.h"
#include "Geometry_2D.h"
#include "Automate_TLM.h"

using namespace std;

const double pi=acos(-1.0L);
const double  u0 = 12566370614e-16;
const double e0 = 88541878176e-22;
const double c = 299792458;
const double Z0 = sqrt(u0/e0);//
//const double er(1),ur(1);
//.....................................................................................................
//constructor  for waveguide
WG_16::WG_16(int simtype, double width, double height, double length, double dl, double tfactor, float er, float sigma_ey,int scn_typ, int n_PMl, int conduct_prof,double Refn_factor)
{

    cout<<" Creating     Waveguide ......"<<endl;
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
        neighbours_set = false;
        scn_type = scn_typ;
        sim_type = simtype;
        centre_node = return_coordinates_WG( width,height,length,npml,dl,width/2,height/2,length/2 );
        WG_width = width;
        WG_height = height;
        WG_length = length;
        WG_dl = dl;
        npml = n_PMl;

        Nx = int ( (WG_width / WG_dl) + 0.5 );
        Nxx = Nx;
        Ny = int ( ( WG_height / WG_dl) + 0.5);
        if( Ny ==1 )cout<<" 2D SCN waveguide mode!! " <<endl;
        Nyy = Ny;
        Nz = int ( (WG_length / WG_dl) + 0.5);
        Nzz = Nz;
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in waveguide

        cout<<"Ntotal:"<<Ntotal<<endl;
        cout<<"Nx * Ny "<< Nx*Ny<<endl;
        //cin>>Nx;

        WG_PML_length = n_PMl * dl;
        Nzz =  2*n_PMl + Nz;
        Ntotal_ = Nzz * Nx *Ny;            // Ntotal + nodes in PML
        cout<<" Ntotal with PML:  "<<Ntotal_<<endl;
        cout<<" Centre node : " <<centre_node<<endl;

        int PML_boundary1 = n_PMl;
        int PML_boundary2 = n_PMl+Nz;
        WG_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nt=Nx*Ny;

        bool special_nodex = false, special_nodey = false, special_nodez=false , spec = false;

        for(int z=0; z<Nzz; z++)
            for (int y=0 ; y<Nyy ; y++)
            {
                for (int x=0 ; x<Nxx ; x++)
                {
                    //Setting flag for nodes at boundaries/ edge of domain
                    special_nodex = (( x==0)||( x==Nxx-1));
                    special_nodey = (( y==0)||( y==Nyy-1));
                    special_nodez = (( z==0)||( z==Nzz-1));
                    spec = special_nodex || special_nodey || special_nodez ;
                    // Creating SCN_nodes;
                    if((z >= PML_boundary1 )&& (z < PML_boundary2 ) )
                    {

                        WG_nodes.push_back(SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,false));
                       //if( id == 38879 ) {cout<<dl*x<<" "<<dl*y<<" "<<dl*z;cin>>z;}

                    }
                    else
                    {
                        if ( z < PML_boundary1) WG_nodes.push_back(SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,-1,-1,PML_boundary1-z-1));

                        else WG_nodes.push_back(SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,-1,-1,z-PML_boundary2)); //cout<<"total_id : "<<id<<endl;cin>>id; }

                    }
                    id = id + 1;
                }
            }
          // cin>>id;// cout<<"id"<<id<<endl;cin>>id;
    }


}

//constructor for 3D problem
WG_16::WG_16(double width, double height, double length, int simtype, double dl, double tfactor, float er, float sigma_ey,int scn_typ, int n_PMl, int conduct_prof,double Refn_factor)
{
    cout<<" Creating Cubic domain ......"<<endl;
    cout<<" dl =  "<< dl<<" length = "<<length<<endl;

    if ( dl <= 0 || int(length) <=0 )
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
        cout<< " Troublesome space discretization value chosen " <<endl;
        cin>>WG_width;
    }
    else
    {   neighbours_set = false;
        sim_type = simtype;
        scn_type = scn_typ;
        WG_width = width;
        WG_height = height;
        WG_length = length;
        WG_dl = dl;
        npml = n_PMl;


        //Number of nodes in the problem domain along each coordinate axis
        Nx = int ( (WG_width / WG_dl)  + 0.5) +1;                     // note the added node which indicates position where the SD is placed.
        Ny = int ( (WG_height / WG_dl) + 0.5) +1;
        Nz = int ( (WG_length / WG_dl) + 0.5) +1;
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in waveguide

        int mid_x = int(0.5*width/dl  + 0.5)+ n_PMl;
        int mid_y = int(0.5*height/dl + 0.5)+ n_PMl;
        int mid_z = int(0.5*length/dl + 0.5)+ n_PMl;

        //Number of nodes in entire computation along each coordinate axis
        Nxx = Nx + 2*n_PMl;
        Nyy = Ny + 2*n_PMl;
        Nzz = Nz + 2*n_PMl;
        Ntotal_ = Nxx*Nyy*Nzz;                  // total number of nodes in waveguide

       //Reserve memory for WG_16 nodes
        //WG_nodes.reserve(Ntotal_);

        cout<<Nxx<<": "<<Nyy<<": "<<Nzz<<endl;
        WG_PML_length = n_PMl * dl;

        int PML_boundaryx1 = n_PMl;
        int PML_boundaryx2 = Nx + n_PMl;
        int PML_boundaryy1 = n_PMl;
        int PML_boundaryy2 = Ny + n_PMl;
        int PML_boundaryz1 = n_PMl;
        int PML_boundaryz2 = Nz + n_PMl;

        //WG_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nxy = Nxx*Nyy;

        cout<<"total node in problem domain:"<<Ntotal<<endl;
        cout<<"total node in whole domain "<< Ntotal_<<endl;
        WG_nodes.resize(Ntotal_);
        //int il =0;
        //cin>>il;

         bool special_nodex = false, special_nodey = false, special_nodez=false , spec = false;

        for(int z=0; z<Nzz; z++){
            for (int y=0 ; y<Nyy ; y++)
            {
                for (int x=0 ; x<Nxx ; x++)
                {
                    //WG_nodes.resize(id+1);
                    //Setting flag for nodes at boundaries/ edge of domain
                    special_nodex = (( x==0)||( x==Nxx-1));
                    special_nodey = (( y==0)||( y==Nyy-1));
                    special_nodez = (( z==0)||( z==Nzz-1));
                    spec = special_nodex || special_nodey || special_nodez ;
                    // Creating SCN_nodes;

                    if(((x >= PML_boundaryx1 )&& (x < PML_boundaryx2 ))&& ((y >= PML_boundaryy1 )&& (y < PML_boundaryy2 )) &&((z >= PML_boundaryz1 )&& (z < PML_boundaryz2 )) )
                    {
                       if( (x == mid_x )&&(z == mid_z)&& (y== mid_y)) { centre_node = id ; cout<<" ID of centre node : "<<centre_node<<" placed in z plane : "<<mid_z<<endl; }
                       //if( (x == mid_x )&&(z == mid_z)&& (y>=n_PMl)&&(y!= mid_y) ) WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,spec,true,false));
                       //else
                       WG_nodes[id] = SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,false);
                      //WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,spec,false,false));
                    }
                    else
                    {
                        int Lx = -1;
                        int Ly = -1;
                        int Lz = -1;

                        if (( x < PML_boundaryx1)||( x >= PML_boundaryx2))
                        {
                            if ( x < PML_boundaryx1 ) Lx = PML_boundaryx1-x-1;
                            else Lx = x-PML_boundaryx2;
                        }
                        if (( y < PML_boundaryy1)||( y >= PML_boundaryy2))
                        {
                            if ( y < PML_boundaryy1 ) Ly = PML_boundaryy1-y-1;
                            else Ly = y-PML_boundaryy2;
                        }
                        if (( z < PML_boundaryz1)||( z >= PML_boundaryz2))
                        {
                            if ( z < PML_boundaryz1 ) Lz = PML_boundaryz1-z-1;
                            else Lz = z-PML_boundaryz2;
                        }

                       // WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz));
                        WG_nodes[id] = (SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz));

                    }
                    id = id + 1;
                }
            }

    }
    }

}

//2D_SCN
WG_16::WG_16(double width, double height, int simtype, double dl, double tfactor, float er, float sigma_ey, int scn_typ, int n_PMl, int conduct_prof,double Refn_factor)
{

    cout<<" Creating Planar domain for SCN ......"<<endl;
    cout<<" dl "<< dl<<" length "<< dl<<endl;

    if ( dl <= 0  )
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
        cout<< " Troublesome space discretization value chosen " <<endl;
        cin>>WG_width;
    }
    else
    {
        neighbours_set = false;
        sim_type = simtype;         //sim type - waveguide, dipole etc.
        scn_type = scn_typ;         // dubard or new
        WG_width = width;
        WG_height = height;
        WG_length = dl;
        WG_dl = dl;
        npml = n_PMl;
        centre_node = return_coordinates_2D(width,height,npml,npml,dl,width/2,height/2,dl);

//Number of nodes in the problem domain along each coordinate axis
        Nx = int ( (WG_width / WG_dl)  + 0.5);
        Ny = int ( (WG_height / WG_dl) + 0.5);
        Nz = int ( (WG_length / WG_dl) + 0.5);
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in computational domain

//Number of nodes in entire computation along each coordinate axis
        Nxx = Nx + 2*n_PMl;
        Nyy = Ny + 2*n_PMl;
        Nzz = Nz;
        Ntotal_ = Nxx*Nyy*Nzz;                  // total number of nodes

//Reserve memory for WG_16 nodes
        //WG_nodes.reserve(Ntotal_);

        cout<<Nxx<<" : "<<Nyy<<" : "<<Nz<<endl;
        WG_PML_length = n_PMl * dl;

        int PML_boundaryx1 = n_PMl;
        int PML_boundaryx2 = Nx + n_PMl;
        int PML_boundaryy1 = n_PMl;
        int PML_boundaryy2 = Ny + n_PMl;

        //WG_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nxy = Nxx*Nyy;

        cout<<"Total node in problem domain:"<<Ntotal<<endl;
        cout<<"Nx * Ny *Nz "<< Ntotal_<<endl;
        WG_nodes.resize(Ntotal_);

        bool special_nodex = false, special_nodey = false, spec = false;

        for(int z=0; z<Nz; z++){
            for (int y=0 ; y<Nyy ; y++)
            {
                for (int x=0 ; x<Nxx ; x++)
                {
                    //WG_nodes.resize(id+1);
                    //Setting flag for nodes at boundaries/ edge of domain
                    special_nodex = ( ( x==0)||( x==Nxx-1) );
                    special_nodey = ( ( y==0)||( y==Nyy-1) );
                    spec = special_nodex || special_nodey ;
                    // Creating SCN_nodes;

                    if(((x >= PML_boundaryx1 )&& (x < PML_boundaryx2 ))&& ((y >= PML_boundaryy1 )&& (y < PML_boundaryy2 )))
                    {
                        //cout<<id<<endl;
                        WG_nodes[id] = SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,false);
                    }
                    else
                    {
                        int Lx = -1;
                        int Ly = -1;
                        int Lz = -1;

                        if (( x < PML_boundaryx1)||( x >= PML_boundaryx2))
                        {
                            if ( x < PML_boundaryx1 ) Lx = PML_boundaryx1-x-1;
                            else Lx = x-PML_boundaryx2;
                        }
                        if (( y < PML_boundaryy1)||( y >= PML_boundaryy2))
                        {
                            if ( y < PML_boundaryy1 ) Ly = PML_boundaryy1-y-1;
                            else Ly = y-PML_boundaryy2;
                        }

                       // WG_nodes.push_back(SCN_node(dl,tfactor,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz));
                        WG_nodes[id] = (SCN_node(scn_typ,dl,tfactor,er,sigma_ey,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz));

                    }
                    id = id + 1;
                }
            }

    }
    }
    cout<<" Centre node : " <<centre_node<<endl;
}


WG_16::WG_16(const WG_16& copy_WG_16)
{
    WG_width = copy_WG_16.WG_width;
    WG_height = copy_WG_16.WG_height;
    WG_length = copy_WG_16.WG_length;
    WG_dl = copy_WG_16.WG_dl;            // discretization length
    WG_PML_length = copy_WG_16.WG_PML_length;
    sim_type = copy_WG_16.sim_type;
    Nx = copy_WG_16.Nx;                 //number of nodes on the x , y and z axis
    Ny = copy_WG_16.Ny;
    Nz = copy_WG_16.Nz;
    Nzz = copy_WG_16.Nzz;
    Nxx = copy_WG_16.Nxx;
    Nyy = copy_WG_16.Nyy;
    npml = copy_WG_16.npml;
    centre_node = copy_WG_16.centre_node;         //node of excitation for short dipole - default set to -1 if waveguide
    WG_Ex = copy_WG_16.WG_Ex;
    WG_Ey = copy_WG_16.WG_Ey;
    Ey_outputs = copy_WG_16.Ey_outputs;
    WG_Ez = copy_WG_16.WG_Ez;
    WG_Hx = copy_WG_16.WG_Hx;
    WG_Hy = copy_WG_16.WG_Hy;
    WG_Hz = copy_WG_16.WG_Hz;
    average_Ey = copy_WG_16.average_Ey;
    WG_Ey_plane_xy = copy_WG_16.WG_Ey_plane_xy;
    WG_Ey_plane_yz = copy_WG_16.WG_Ey_plane_yz;
    WG_Ey_plane_zx = copy_WG_16.WG_Ey_plane_zx;
    Ey_field_along_line = copy_WG_16.Ey_field_along_line;
    Ntotal= copy_WG_16.Ntotal;
    Ntotal_ = copy_WG_16.Ntotal_;
    WG_nodes = copy_WG_16.WG_nodes;   // represents the waveguide nodes; this is set in the constructor depending on the simulation parameters.
    WG_node_neighbours= copy_WG_16.WG_node_neighbours;
    check_meshed = copy_WG_16.check_meshed;         //checks if the waveguide has been meshed / filled with cubes
    neighbours_set = copy_WG_16.neighbours_set;
    scn_type = copy_WG_16.scn_type;
    V_incident_xy = copy_WG_16.V_incident_xy;
}

WG_16 &WG_16::operator=(const WG_16& copy_WG_16)
{
    if (this!= &copy_WG_16)
    {
    WG_width = copy_WG_16.WG_width;
    WG_height = copy_WG_16.WG_height;
    WG_length = copy_WG_16.WG_length;
    WG_dl = copy_WG_16.WG_dl;            // discretization length
    WG_PML_length = copy_WG_16.WG_PML_length;
    sim_type = copy_WG_16.sim_type;
    Nx = copy_WG_16.Nx;                 //number of nodes on the x , y and z axis
    Ny = copy_WG_16.Ny;
    Nz = copy_WG_16.Nz;
    Nzz = copy_WG_16.Nzz;
    Nxx = copy_WG_16.Nxx;
    Nyy = copy_WG_16.Nyy;
    npml = copy_WG_16.npml;
    centre_node = copy_WG_16.centre_node;         //node of excitation for short dipole - default set to -1 if waveguide
    WG_Ex = copy_WG_16.WG_Ex;
    WG_Ey = copy_WG_16.WG_Ey;
    Ey_outputs = copy_WG_16.Ey_outputs;
    WG_Ez = copy_WG_16.WG_Ez;
    WG_Hx = copy_WG_16.WG_Hx;
    WG_Hy = copy_WG_16.WG_Hy;
    WG_Hz = copy_WG_16.WG_Hz;
    average_Ey = copy_WG_16.average_Ey;
    WG_Ey_plane_xy = copy_WG_16.WG_Ey_plane_xy;
    WG_Ey_plane_yz = copy_WG_16.WG_Ey_plane_yz;
    WG_Ey_plane_zx = copy_WG_16.WG_Ey_plane_zx;
    Ey_field_along_line = copy_WG_16.Ey_field_along_line;
    Ntotal= copy_WG_16.Ntotal;
    Ntotal_ = copy_WG_16.Ntotal_;
    WG_nodes = copy_WG_16.WG_nodes;   // represents the waveguide nodes; this is set in the constructor depending on the simulation parameters.
    WG_node_neighbours= copy_WG_16.WG_node_neighbours;
    check_meshed = copy_WG_16.check_meshed;         //checks if the waveguide has been meshed / filled with cubes
    neighbours_set = copy_WG_16.neighbours_set;
    scn_type = copy_WG_16.scn_type;
    V_incident_xy = copy_WG_16.V_incident_xy;
    }
    return (*this);
}

 WG_16::~WG_16()
{
  cout<<"destroyed"<<endl;//destructor
}

double WG_16::Ex_output_at_node(int node_id)
{
    if ( scn_type==1)
    {
    double Gxy(this->WG_nodes[node_id].get_G(1)), Gxz(this->WG_nodes[node_id].get_G(2));

    double Yx(this->WG_nodes[node_id].get_Y(0));

    double Exy =  2*(this->WG_nodes[node_id].get_incidence(1)+this->WG_nodes[node_id].get_incidence(12) + (Yx+2)*this->WG_nodes[node_id].get_incidence(13)-2*this->WG_nodes[node_id].get_incidence(14) );
    Exy = Exy/ ( Yx+Gxy+4);

    double Exz =  2*(this->WG_nodes[node_id].get_incidence(2)+this->WG_nodes[node_id].get_incidence(9) + (Yx+2)*this->WG_nodes[node_id].get_incidence(14)-2*this->WG_nodes[node_id].get_incidence(13) );
    Exz = Exz/ ( Yx+Gxz+4);

    return ((Exy + Exz)/(-1e-3*this->WG_dl));
    }
     if(scn_type ==3)
    {
        double Gxy(this->WG_nodes[node_id].get_G(1)), Gxz(this->WG_nodes[node_id].get_G(2)) ;
        double Axy = 4/(4+ Gxy) ,Axz = 4/(4+ Gxz) ;
        double Yy(0);//this->Y[0]);//,Zy(this->Z[0]);
        double Yz(0);//this->Y[0]);//,Zz(this->Z[0]);

        Yy = 4*( this->WG_nodes[node_id].get_alpha_y() - 0.5);
        Yz = 4*( this->WG_nodes[node_id].get_alpha_z() - 0.5);

        float Cy = 0.5/this->WG_nodes[node_id].get_alpha_y();
        float Cz = 0.5/this->WG_nodes[node_id].get_alpha_z();

        double Exy =  Axy*Cy*(this->WG_nodes[node_id].get_incidence(1)+this->WG_nodes[node_id].get_incidence(12) + Yy*this->WG_nodes[node_id].get_incidence(13)-2*this->WG_nodes[node_id].get_incidence(14) );

        double Exz =  Axz*Cz*(this->WG_nodes[node_id].get_incidence(2)+this->WG_nodes[node_id].get_incidence(9)  + Yz*this->WG_nodes[node_id].get_incidence(14)-2*this->WG_nodes[node_id].get_incidence(13) );
        return ((Exy + Exz)/(-1e-3*this->WG_dl));
    }

    return 0;
}

double WG_16::Ey_output_at_node(int node_id)
{
    if ( scn_type==1)
    {
    double Gyx(this->WG_nodes[node_id].get_G(0));

    double Gyz(this->WG_nodes[node_id].get_G(2));

    //double Gyx(0), Gyz(0);
    double Yy(this->WG_nodes[node_id].get_Y(1));       // note: this is assuming the grid dimensions are all equal to dl

    double Eyx =  2*(this->WG_nodes[node_id].get_incidence(3)+this->WG_nodes[node_id].get_incidence(11) + (Yy+2)*this->WG_nodes[node_id].get_incidence(15)-2*this->WG_nodes[node_id].get_incidence(16) );
    Eyx = Eyx/ ( Yy+Gyx+4);

    double Eyz =  2*(this->WG_nodes[node_id].get_incidence(4)+ this->WG_nodes[node_id].get_incidence(8) + (Yy+2)*this->WG_nodes[node_id].get_incidence(16)-2*this->WG_nodes[node_id].get_incidence(15) );
    Eyz = Eyz/ ( Yy+Gyz+4);

    //cout<<Gyx<<endl<<Gyz<<endl;cin>>Eyz;

    return ((Eyx + Eyz)/(-1e-3*this->WG_dl));
      }

    else if (scn_type==2)
    {
        if(!this->WG_nodes[node_id].is_PML_SCN())
        {
        double Ey =  -0.5*(this->WG_nodes[node_id].get_incidence(3)+this->WG_nodes[node_id].get_incidence(11) + this->WG_nodes[node_id].get_incidence(4)+this->WG_nodes[node_id].get_incidence(8) );

        return (Ey /(1e-3*this->WG_dl) );
        }
        else return (0);
    }

    if(scn_type ==3)
    {
        double Gyx(this->WG_nodes[node_id].get_G(0)), Gyz(this->WG_nodes[node_id].get_G(2)) ;
        double Ayx = 4/(4+ Gyx) ,Ayz = 4/(4+ Gyz) ;
        double Yx(0);//this->Y[0]);//,Zy(this->Z[0]);
        double Yz(0);//this->Y[0]);//,Zz(this->Z[0]);

        Yx = 4*( this->WG_nodes[node_id].get_alpha_x() - 0.5);
        Yz = 4*( this->WG_nodes[node_id].get_alpha_z() - 0.5);

        float Cx = 0.5/this->WG_nodes[node_id].get_alpha_x();
        float Cz = 0.5/this->WG_nodes[node_id].get_alpha_z();

        double Eyx =  Ayx*Cx*(this->WG_nodes[node_id].get_incidence(3)+this->WG_nodes[node_id].get_incidence(11) + Yx*this->WG_nodes[node_id].get_incidence(15)-2*this->WG_nodes[node_id].get_incidence(16) );

        double Eyz =  Ayz*Cz*(this->WG_nodes[node_id].get_incidence(4)+ this->WG_nodes[node_id].get_incidence(8) + Yz*this->WG_nodes[node_id].get_incidence(16)-2*this->WG_nodes[node_id].get_incidence(15) );

        return (-(Eyx + Eyz))/(1e-3*this->WG_dl) ;
    }

    return 0;
}



double WG_16::Ez_output_at_node(int node_id)
{

    if ( scn_type==1)
    {
    double Gzy(this->WG_nodes[node_id].get_G(1)), Gzx(this->WG_nodes[node_id].get_G(0));

    double Yz(this->WG_nodes[node_id].get_Y(2));

    double Ezx =  2*(this->WG_nodes[node_id].get_incidence(6)+ this->WG_nodes[node_id].get_incidence(10) + (Yz+2)*this->WG_nodes[node_id].get_incidence(17)-2*this->WG_nodes[node_id].get_incidence(18) );
    Ezx = Ezx/ ( Yz + Gzx+4);

    double Ezy =  2*(this->WG_nodes[node_id].get_incidence(5)+ this->WG_nodes[node_id].get_incidence(7) + (Yz+2)*this->WG_nodes[node_id].get_incidence(18)-2*this->WG_nodes[node_id].get_incidence(17) );
    Ezy = Ezy/ ( Yz + Gzy+4);

    return ((Ezx + Ezy)/(-1e-3*this->WG_dl));
    }

     else if( scn_type==2 )
    {

    double Ez =  -0.5*(this->WG_nodes[node_id].get_incidence(5)+this->WG_nodes[node_id].get_incidence(7) + this->WG_nodes[node_id].get_incidence(6)+this->WG_nodes[node_id].get_incidence(10) );
    //cout<< Ez<<endl;
    return (Ez /(1e-3*this->WG_dl) );

    }
    if(scn_type ==3)
    {
        double Gzx(this->WG_nodes[node_id].get_G(0)), Gzy(this->WG_nodes[node_id].get_G(1)) ;
        double Azx = 4/(4+ Gzx) ,Azy = 4/(4+ Gzy) ;
        double Yx(0);//this->Y[0]);//,Zy(this->Z[0]);
        double Yy(0);//this->Y[0]);//,Zz(this->Z[0]);

        Yx = 4*( this->WG_nodes[node_id].get_alpha_x() - 0.5);
        Yy = 4*( this->WG_nodes[node_id].get_alpha_y() - 0.5);

        float Cx = 0.5/this->WG_nodes[node_id].get_alpha_x();
        float Cy = 0.5/this->WG_nodes[node_id].get_alpha_y();

        double Ezx =  Azx*Cx*(this->WG_nodes[node_id].get_incidence(6)+ this->WG_nodes[node_id].get_incidence(10) + Yx*this->WG_nodes[node_id].get_incidence(17)-2*this->WG_nodes[node_id].get_incidence(18) );

        double Ezy =  Azy*Cy*(this->WG_nodes[node_id].get_incidence(5)+ this->WG_nodes[node_id].get_incidence(7) + Yy*this->WG_nodes[node_id].get_incidence(18)-2*this->WG_nodes[node_id].get_incidence(17) );

        return ((Ezx + Ezy)/(-1e-3*this->WG_dl));
    }
    return 0;
}

double WG_16::Hx_output_at_node(int node_id)
{
    if ( scn_type==1)
    {
    double Rxy(this->WG_nodes[node_id].get_R(1)), Rxz(this->WG_nodes[node_id].get_R(2));

    double Zx(this->WG_nodes[node_id].get_Z(0));

    double Hxy =  2*(this->WG_nodes[node_id].get_incidence(5)- this->WG_nodes[node_id].get_incidence(7) + (1+2/Zx)*this->WG_nodes[node_id].get_incidence(19)-2*this->WG_nodes[node_id].get_incidence(20)/Zx );
    Hxy = Hxy/ ( Zx+Rxy+4);

    double Hxz =  2*(this->WG_nodes[node_id].get_incidence(8)- this->WG_nodes[node_id].get_incidence(4) + (1+2/Zx)*this->WG_nodes[node_id].get_incidence(20)-2*this->WG_nodes[node_id].get_incidence(19)/Zx );
    Hxz = Hxz/ ( Zx+Rxz+4);

    return (-(Hxy + Hxz)/(Z0*1e-3*this->WG_dl));
    }

    if( scn_type ==2)
    {
     double Hx =  0.5*( this->WG_nodes[node_id].get_incidence(5)- this->WG_nodes[node_id].get_incidence(7) + this->WG_nodes[node_id].get_incidence(8)- this->WG_nodes[node_id].get_incidence(4) );

     return  ( -Hx/(Z0*1e-3*this->WG_dl));
    }
    if(scn_type ==3)
    {
        double Rxy(this->WG_nodes[node_id].get_R(1)), Rxz(this->WG_nodes[node_id].get_R(2)) ;
        double Amxy = 4/(4+ Rxy) ,Amxz = 4/(4+ Rxz) ;
        double Zy(0);//this->Y[0]);//,Zy(this->Z[0]);
        double Zz(0);//this->Y[0]);//,Zz(this->Z[0]);

        Zz = 4*( this->WG_nodes[node_id].get_alpha_z() - 0.5);
        Zy = 4*( this->WG_nodes[node_id].get_alpha_y() - 0.5);

        float Dy = 0.5/this->WG_nodes[node_id].get_alpha_y();
        float Dz = 0.5/this->WG_nodes[node_id].get_alpha_z();

        double Hxy =  Amxy*Dy*(this->WG_nodes[node_id].get_incidence(5)- this->WG_nodes[node_id].get_incidence(7) + Zy*this->WG_nodes[node_id].get_incidence(19)-2*this->WG_nodes[node_id].get_incidence(20));

        double Hxz =  Amxz*Dz*(this->WG_nodes[node_id].get_incidence(8)- this->WG_nodes[node_id].get_incidence(4) + Zz*this->WG_nodes[node_id].get_incidence(20)-2*this->WG_nodes[node_id].get_incidence(19) );

        return ((Hxy + Hxz)/(Z0*1e-3*this->WG_dl));
    }
    return 0;
}

double WG_16::Hy_output_at_node(int node_id)
{
    if ( scn_type==1)
    {
    double Ryx(this->WG_nodes[node_id].get_R(0)), Ryz(this->WG_nodes[node_id].get_R(2));

    double Zy(this->WG_nodes[node_id].get_Z(0));

    double Hyx =  2*(this->WG_nodes[node_id].get_incidence(10)- this->WG_nodes[node_id].get_incidence(6) + (1+2/Zy)*this->WG_nodes[node_id].get_incidence(21)-2*this->WG_nodes[node_id].get_incidence(22)/Zy );
    Hyx = Hyx/ ( Zy+Ryx+4);

    double Hyz =  2*(this->WG_nodes[node_id].get_incidence(2)- this->WG_nodes[node_id].get_incidence(9) + (1+2/Zy)*this->WG_nodes[node_id].get_incidence(22)-2*this->WG_nodes[node_id].get_incidence(21)/Zy);
    Hyz = Hyz/ ( Zy+Ryz+4);

    return (-(Hyx + Hyz)/(Z0*1e-3*this->WG_dl));
    }
    return 0;
}

double WG_16::Hz_output_at_node(int node_id)
{
    if ( scn_type==1)
    {
    double Rzy(this->WG_nodes[node_id].get_R(1)), Rzx(this->WG_nodes[node_id].get_R(0));

    double Zz(this->WG_nodes[node_id].get_Z(0));

    double Hzx =  2*(this->WG_nodes[node_id].get_incidence(3)- this->WG_nodes[node_id].get_incidence(11) + (1+2/Zz)*this->WG_nodes[node_id].get_incidence(23)-2*this->WG_nodes[node_id].get_incidence(24)/Zz);
    Hzx = Hzx/ ( Zz+Rzx+4);

    double Hzy =  2*(this->WG_nodes[node_id].get_incidence(12)- this->WG_nodes[node_id].get_incidence(1) + (1+2/Zz)*this->WG_nodes[node_id].get_incidence(24)-2*this->WG_nodes[node_id].get_incidence(23)/Zz);
    Hzy = Hzy/ ( Zz+Rzy+4);

    return (-(Hzx + Hzy)/(Z0*1e-3*this->WG_dl));
    }
    if(scn_type ==3)
    {
        double Rzx(this->WG_nodes[node_id].get_R(0)), Rzy(this->WG_nodes[node_id].get_R(1)) ;
        double Amzx = 4/(4+ Rzx) ,Amzy = 4/(4+ Rzy) ;
        double Zx(0);//this->Y[0]);//,Zy(this->Z[0]);
        double Zy(0);//this->Y[0]);//,Zz(this->Z[0]);

        Zx = 4*( this->WG_nodes[node_id].get_alpha_x() - 0.5);
        Zy = 4*( this->WG_nodes[node_id].get_alpha_y() - 0.5);

        float Dx = 0.5/this->WG_nodes[node_id].get_alpha_x();
        float Dy = 0.5/this->WG_nodes[node_id].get_alpha_y();

        double Hzx =  Amzx*Dx*(this->WG_nodes[node_id].get_incidence(3)- this->WG_nodes[node_id].get_incidence(11) + Zx*this->WG_nodes[node_id].get_incidence(23)-2*this->WG_nodes[node_id].get_incidence(24));
        double Hzy =  Amzy*Dy*(this->WG_nodes[node_id].get_incidence(12)- this->WG_nodes[node_id].get_incidence(1) + Zy*this->WG_nodes[node_id].get_incidence(24)-2*this->WG_nodes[node_id].get_incidence(23));

        return (-(Hzx + Hzy)/(Z0*1e-3*this->WG_dl));
    }
    return 0;
    //return (0.5*(this->WG_nodes[node_id].get_incidence(1) - this->WG_nodes[node_id].get_incidence(3) + this->WG_nodes[node_id].get_incidence(11) -this->WG_nodes[node_id].get_incidence(12))/(Z0*1e-3*this->WG_dl) );
}

double WG_16::compute_field_average_in_z_plane(int z_plane)
{
     //cout<<" Plane : " << z_plane << endl;
     int node_id = z_plane* this->Nxx*Nyy;
     double field_value = 0;
     int test =0;

     for( int y = 0; y< Nyy ; y++)
        for ( int x = 0; x< Nxx; x++)
        {
            field_value = field_value + Ey_output_at_node(node_id);//can test using print par?
            node_id = node_id+1;
        }

       // cout<<"...test..."<<test<<endl;
       // cout<<"...total... "<<Nxx*Nyy<<endl;

    return (field_value/(Nxx*Nyy));
}

vector<double> WG_16::compute_field_along_line( int x_plane_id, int y_plane, int z_plane_id)
{
    int start_node = x_plane_id + (y_plane)*Nxx + (npml+0)*Nxx*Nyy;
    double field = 0;
    vector <double> fields_along_line;

    for(int i = 0 ; i < Nz; i++)
    {
      field = Ey_output_at_node(start_node);
      //cout<< " node iD: " <<start_node<<" "<< WG_nodes[start_node].is_PML_SCN() << " ";
      fields_along_line.push_back(field);
      start_node = start_node +(Nxx*Nyy);
    }

    //cout<<endl;
    return( fields_along_line);
}

int WG_16::get_true_centre_node()
{
    return(centre_node);
}

void WG_16::set_special_nodes()
{
    if(neighbours_set == false)
    {   int temp;
        cout<<" neighbours have not been set!!" <<endl;
        cout<<" press any number to proceed...neighbours will now be set"<<endl;
        cin>>temp;
        this->set_WG_neighbour();               // potential problems could arise here when parallelization of the code
    }

    for (int wg_i = 0; wg_i<this->WG_nodes.size(); wg_i++)
    {
    int numbr_neigh = WG_node_neighbours[wg_i].size();

        for (int i =0; i<numbr_neigh;i++)
        {
            if ( (  this->WG_nodes[ WG_node_neighbours[wg_i][i] ].is_PEC() ) &&   (   this->WG_nodes[wg_i].check_special_node())  )
            {
                int temp;
                //cout<<" |[][x] A PEC node is one node away from boundary!!"<<endl;
                //  cout<<"This will create problems in the connection process!!"<<endl;
                // cout<<"Suggestion is to extend the domain in the " << i <<"direction to create |[][][x]"<<endl;

            }
            if ( (  this->WG_nodes[ WG_node_neighbours[wg_i][i] ].is_PEC() ) &&   (  !this->WG_nodes[wg_i].check_special_node())  )
            {
                this->WG_nodes[wg_i].set_special_SCN_node(true);
                break;
            }
        }
    }
}


void WG_16::place_dipole_antenna(double dipole_length, int x, int y, int z )
{
        int sd_node = this->centre_node;
        int antenna_node_no = 0.5*dipole_length/this->WG_dl + 0.5;
        int max_l = ( Ny-1 );

        cout<<endl<<" Inserting Dipole Antenna of length "<<0.5*dipole_length<<endl<<endl;

        if (antenna_node_no>max_l)
        {
            cout<<"CAUTION: Antenna length bigger than geometry!!"<<endl;
            cout<<"Setting to full length mode"<<endl;
            antenna_node_no = max_l;
        }

        for(int i = 1; i<=int(antenna_node_no*0.5+0.5) ; i++)
        {
           cout<<"Location ID: "<< sd_node -i*Nxx<<" ; "<<sd_node +i*Nxx <<endl;
           this->WG_nodes[sd_node +i*Nxx].set_PEC();
           this->WG_nodes[sd_node -i*Nxx].set_PEC();

        }

         int zz = this->WG_nodes[sd_node].get_coord(2);
         this->print_WG_plane_par(zz,3);
}

void WG_16::insert_iris(int thickness , float z_plane)
{
    if (2*thickness < Nyy)
    {
        cout<<" Inserting Iris......in ";
        cout<< " Plane : "<< z_plane <<endl;
        double dl = WG_dl;
        int start_node = return_coordinates_WG (WG_width,WG_height,WG_length,npml,dl,0,0,z_plane); // on the left
        cout<<"start node bottom"<<start_node<<endl;
        int width_node = start_node + Nxx*Nyy-1; // on the right
        cout<<"start node top "<<width_node<<endl;

        for(int i = 0; i<thickness;i++)
        {
            int ii=start_node + i*Nxx;
            int jj=width_node - i*Nxx;
            for(int counter =0 ; counter<Nx; counter++)
            {
                WG_nodes[ii].set_PEC();
                WG_nodes[jj].set_PEC();
                ii = ii + 1;
                jj = jj - 1;
            }
        }
    }

    else cout<<"Cannot insert Iris "<<endl;
}


void WG_16::insert_perfect_cube(int cube_l,int nx0,int ny0, int nz0) // cube_l == full dimension of cube
{
    int Nxy = Nxx*Nyy;
    int ref_loc = nx0 + npml + (npml+ny0)*Nxx + (npml+nz0)*Nxy;
    cout<<" reference location of PEC cube "<<ref_loc<<endl;

    float ref_x =nx0 + cube_l, ref_y =ny0 + cube_l ,ref_z =nz0 + cube_l;
    int loc1 = ref_loc;
    int loc2 = 0;

// simple test is cout<< loc1<<endl; cout<<
    if( (ref_x < Nx)&&( ref_y < Ny)&&( ref_z < Nz) )
    {
        for(int k = 0; k<cube_l;k++)// z
        {
            loc1 = ref_loc + k*Nxy;
            for(int j =0;j<cube_l;j++)   // y
            {
                loc2 = loc1 + j*Nxx;
                for(int i =0;i<cube_l;i++)  //x
                {
                    WG_nodes[loc2 + i].set_PEC();   //
                }
            }
        }
    }
    else
    {
        int a=0;
        cout<< "...cannot insert cube ... change dimensions..."<<endl;
        cout<<"x length "<<ref_x << " y length "<<ref_y << "z length "<<ref_z<<endl;
        cin>>a;
    }
}


void WG_16::insert_FSS_square(double dims, float pos_z)
{
    cout<< "....INSERTING FREQUENCY SELECTIVE SURFACE ( SQUARE ) IN PLANE...." <<pos_z<<endl<<endl;
    int start_node = return_coordinates_WG (WG_width,WG_height,WG_length,npml,WG_dl,0,0,pos_z);
    int start_node2= start_node;

    int N_dims_B = int(dims/WG_dl + 0.5); // bottom limit
    int N_dims_T = Nyy-N_dims_B;
    int N_dims_R = Nxx-N_dims_B;

    for(int i=0; i<Nyy; i++)
        for( int j =0; j<Nxx; j++)
        {
            start_node2 = start_node + Nxx*i;

            if(( i<N_dims_B) ||( i>=N_dims_T)) this->WG_nodes[start_node2 + j].set_PEC();
            else
                if(( j<N_dims_B) || (j>=N_dims_R)) this->WG_nodes[start_node2 + j].set_PEC();
        }

    int zz = this->WG_nodes[start_node].get_coord(2);
    this->print_WG_plane_par(zz,3);

}

void WG_16::insert_FSS_jerusalem_cross(double dims, float pos_z)
{
    cout<< "....INSERTING FREQUENCY SELECTIVE SURFACE (JC) IN PLANE...." <<pos_z<<endl<<endl;
    int start_node = return_coordinates_WG (WG_width,WG_height,WG_length,npml,WG_dl,0,0,pos_z);

    //cout<< " dims "<< Nxx << " .... " << Nyy <<endl;

    if ( (Nxx*Nyy % 2)==0 )         // change dimensions
    {

            cout<< " CHANGE DIMENSIONS TO SUPPORT INSERTION OF JC!! "<<endl;
            cout<< Nxx << " x " << Nyy <<endl;
            cout<< " CHANGE L or W to odd number of cells" <<endl;
            cout<<"  TERMINATING.........."<<endl;
            cin>> dims;
            return;
    }

    int mid_x = int(Nxx*.5 + 0.5); //
    int mid_y = int(Nyy*.5 + 0.5);
    int N_dims = int(dims/WG_dl +0.5);
    int p = int(3/WG_dl +0.5);

    cout<< " p number of nodes : " << p <<endl;

    start_node = start_node+ mid_x-1 + (mid_y-1)*Nxx; // finds mid point

    int i=1,j=1;
    this->WG_nodes[ start_node ].set_PEC(); //

    while(i<=N_dims)
    {
        this->WG_nodes[start_node + i].set_PEC(); // along x
        this->WG_nodes[start_node - i].set_PEC(); // along x

        this->WG_nodes[start_node + i*Nxx].set_PEC(); // along y
        this->WG_nodes[start_node - i*Nxx].set_PEC(); // along y

       if ( i == N_dims)
       {
            while(j<=p)
            {

                this->WG_nodes[start_node + i + j*Nxx].set_PEC(); // moving along x
                this->WG_nodes[start_node + i - j*Nxx].set_PEC();

                this->WG_nodes[start_node - i + j*Nxx].set_PEC();
                this->WG_nodes[start_node - i - j*Nxx].set_PEC();

                this->WG_nodes[start_node + i*Nxx + j].set_PEC(); // moving along y
                this->WG_nodes[start_node - i*Nxx + j].set_PEC(); //

                this->WG_nodes[start_node + i*Nxx - j].set_PEC();
                this->WG_nodes[start_node - i*Nxx - j].set_PEC();

                j = j+1;
            }
       }
       i = i+1;
    }
    int zz = this->WG_nodes[start_node].get_coord(2);
    this->print_WG_plane_par(zz,3);
}

void WG_16::insert_square_loop(double sq_length, int x, int y, int z )
{
        int sd_node = this->centre_node;// - 0.125*sq_length/this->WG_dl;

        int zz = this->WG_nodes[sd_node].get_coord(2);
        this->print_WG_plane_par(zz,3);

        int total_no = 0.5*sq_length/this->WG_dl + 0.5;
        int max_l = ( Ny-1 );

        cout<<endl<<" Inserting Square loop of length = "<<total_no<<endl<<endl;

        if ((total_no>max_l)||(total_no > Nx-1))
        {
            cout<<"CAUTION: Square loop length bigger than geometry!!"<<endl;
            cout<<"Setting to full length mode"<<endl;
             total_no = max_l;
        }

        for(int i = 1; i<=int(total_no*0.5+0.5) ; i++)
        {
           cout<<"Location ID: "<<endl<< sd_node -i*Nxx<<" <-- -ve y--> "<<sd_node -i*Nxx + total_no <<" | "<<endl
           <<sd_node +i*Nxx <<"<-- +ve y -->"<<sd_node -i*Nxx + total_no<<endl<<endl;

           this->WG_nodes[sd_node +i*Nxx].set_PEC();
           this->WG_nodes[sd_node +i*Nxx + total_no].set_PEC();

           this->WG_nodes[sd_node -i*Nxx].set_PEC();
           this->WG_nodes[sd_node -i*Nxx + total_no].set_PEC();

           if(i==1) //inserting PEC
           {
               this->WG_nodes[sd_node + total_no].set_PEC();
           }

           if(i==int(total_no*0.5+0.5)) // inserting PEC along width
           {
               for(int j=1; j<total_no; j++)
               {
                   this->WG_nodes[sd_node +i*Nxx + j].set_PEC();
                   this->WG_nodes[sd_node -i*Nxx + j].set_PEC();
               }
           }
        }

        this->print_WG_plane_par(zz,3);
}

void WG_16::insert_PML(int num_PML, double loc, double Reflctn_f, int c_profile, int scn_typ)
{
        if( loc > WG_length)
        {
            cout<< " error !" <<endl;
            cout<< " attempting to insert PML in location beyond waveguide length" <<endl;
            cout<< " default position -> centre ... will be set" <<endl;

            loc = WG_length/2;
        }

        int start_node = return_coordinates_WG (WG_width,WG_height,WG_length,npml,WG_dl,0,0,loc);
        int start_node2 = start_node;

        int z_loc = this->WG_nodes[start_node].get_coord(2);

        this->print_WG_plane_par(z_loc,2);  // print the PML properties of the node

        //int max_l = ( Ny-1 );

        scn_type = scn_typ;

        cout<<endl<<" Inserting " << scn_type<< " PML in centre of domain with NUMBER OF LAYERS = "<<num_PML<<endl<<endl;

        int total_no = Nx*Ny;

        if( scn_type < 1)
        {
            cout<< " setting SCN TYPE to default... mapped SCN-PML " <<endl;
           // scn_type = 2;
        }

        for(int i=0; i<num_PML; i++)
        {
            start_node2 = start_node + Nx*Ny*i;

            for( int j =0; j<total_no; j++)
            {
                this->WG_nodes[start_node2+j].set_PML(i,Reflctn_f,c_profile, scn_type, num_PML);
            }
            cout<<endl;
        }

        this->print_WG_plane_par(z_loc,1);
        this->print_WG_plane_par(z_loc+1,1);
        this->print_WG_plane_par(z_loc+2,1);


        //cin>> num_PML;

}


void WG_16::set_WG_neighbour()
{
    int id = 0;
    int Nt = this->Nxx*this->Nyy;

    WG_node_neighbours.resize(this->Ntotal_);

    for(int z=0; z<this->Nzz; z++)
    {
        //cout<<id<<endl;
        for (int y=0 ; y<this->Nyy ; y++)
        {
            for (int x=0 ; x<this->Nxx ; x++)
            {
                if (x!=0)
                {   //cout<<WG_node_neighbours.size()<<endl;
                //cin>>id;
                    WG_node_neighbours[id].push_back(WG_nodes[id-1].get_iD());
                    WG_node_neighbours[id-1].push_back(WG_nodes[id].get_iD());
                }
                if (y!=0)
                {
                    WG_node_neighbours[id].push_back(WG_nodes[id-this->Nxx].get_iD());
                    WG_node_neighbours[id-this->Nxx].push_back(WG_nodes[id].get_iD());
                }
                if(z!=0)
                {
                    WG_node_neighbours[id].push_back(WG_nodes[id-Nt].get_iD());
                    WG_node_neighbours[id-Nt].push_back(WG_nodes[id].get_iD());
                }
                id = id+1;
            }
        }
    }
    neighbours_set = true;
    cout<<"neighbours set"<<endl;

}

void WG_16::print_SCN_neighbour(int node_id)
{
    int countt = WG_node_neighbours[node_id].size();
    for ( int i = 0; i<countt; i++)
        cout<<WG_node_neighbours[node_id][i];
}

void WG_16::print_WG_nodes()
{
    for( int i=0; i< this->WG_nodes.size(); i++)
    {
        cout<<this->WG_nodes[i]<<endl;
    }
}


void WG_16::print_WG_node(int _node)
{
    cout<<" Node : " <<_node << endl;
    for(int i=0; i<12; i++)
    {
        cout<<this->WG_nodes[_node].get_incidence(i+1)<<endl;
    }

}

void WG_16::print_WG_plane_par(int z_plane, int par )
{
    cout<<" Plane : " << z_plane << endl;

 int node_id = (z_plane) * this->Nxx*Nyy;

 for( int y = 0; y< Nyy ; y++){
    for ( int x = 0; x< Nxx; x++){
        if( par == 1) cout<<this->WG_nodes[node_id].is_PML_SCN()<<" :  ";
        else if( par == 2) cout<<this->WG_nodes[node_id].get_iD()<<" :  ";
        else if( par == 3) cout<<this->WG_nodes[node_id].is_PEC()<<" :  ";
        else if( par == 4) cout<<this->WG_nodes[node_id].get_G(0)<<"," <<this->WG_nodes[node_id].get_G(1)<<","<<this->WG_nodes[node_id].get_G(2) <<": ";
        node_id = node_id+1;
    }
    cout<<endl<<endl;
 }
}


void WG_16::breakpoint_(int i)
{
    cout<<" breakpoint " << i <<endl;
    cin>>i;
}

void WG_16::print_this_neighbour(int excited_node)
{
    int ii = this->WG_node_neighbours[excited_node].size();

    cout<<" Node : " << excited_node << endl;
    for( int i=1; i<=12; i++)  cout<<this->WG_nodes[excited_node].get_reflected(i)<<endl;
    for (int j =0 ; j<ii ; j++)
    {
        cout<<" Neighbour "<< this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_iD() <<" : "<<endl;
        cout<<"1: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(1)<<endl;
        cout<<"2: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(2)<<endl;
        cout<<"3: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(3)<<endl;
        cout<<"4: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(4)<<endl;
        cout<<"5: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(5)<<endl;
        cout<<"6: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(6)<<endl;
        cout<<"7: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(7)<<endl;
        cout<<"8: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(8)<<endl;
        cout<<"9: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(9)<<endl;
        cout<<"10: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(10)<<endl;
        cout<<"11: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(11)<<endl;
        cout<<"12: "<<this->WG_nodes[ WG_node_neighbours[excited_node][j] ].get_incidence(12)<<" ? "<<endl;

        //cout<<this->WG_node_neighbours[excited_node][j].get_reflected()<<endl;
    }
}

void WG_16::TLM_excitation_te10(int node , double v_inc)
{
        double start_node = 0;
        double v = 0, v_ = 0;

        for(int j = 0; j<Nyy; j++)
        {
            start_node = node + j*Nxx;
            if((!WG_nodes[start_node].is_PML_SCN())||(!WG_nodes[start_node].is_PEC()))
            {

                v = WG_nodes[start_node].get_incidence(3);
                v_ = v_inc + v;
                WG_nodes[start_node].set_incidence(3,v_);

                v = WG_nodes[start_node].get_incidence(11);
                v_ = v_inc + v;
                WG_nodes[start_node].set_incidence(11,v_);

                v = WG_nodes[start_node].get_incidence(4);
                v_= v_inc + v;
                WG_nodes[start_node].set_incidence(4,v_);

                v = WG_nodes[start_node].get_incidence(8);
                v_= v_inc + v;
                WG_nodes[start_node].set_incidence(8,v_);

                if( scn_type == 1){

                    v = WG_nodes[start_node].get_incidence(15);
                    v_ = v_inc + v;
                    WG_nodes[start_node].set_incidence(15,v_);

                    v = WG_nodes[start_node].get_incidence(16);
                    v_ = v_inc + v;
                    WG_nodes[start_node].set_incidence(16,v_);

                    }
            }
            else cout<< " TROUBLESOME EXCITATION CHECK PEC / PML " <<endl;
        }

}

void WG_16::TLM_excitation_single_node(int node , double v_inc, int axis,int source_res)
{

    double v = 0, v_ = 0;
    int start_node = node;

    //EXCITE EY field component
    if( axis ==2){
/*
        if( source_res > 0)
            {
                double Vy=0,Vy_prime=0, delta_Vy;
                this->WG_nodes[node].compute_source_voltage_y(v_inc,source_res,Vy,Vy_prime);
                delta_Vy = Vy_prime - Vy;
                v_inc = delta_Vy;
            }
*/
    v = WG_nodes[start_node].get_incidence(3);
    v_ = v_inc + v;
    WG_nodes[start_node].set_incidence(3,v_);

    v = WG_nodes[start_node].get_incidence(11);
    v_ = v_inc + v;
    WG_nodes[start_node].set_incidence(11,v_);

    v = WG_nodes[start_node].get_incidence(4);
    v_= v_inc + v;
    WG_nodes[start_node].set_incidence(4,v_);

    v = WG_nodes[start_node].get_incidence(8);
    v_= v_inc + v;
    WG_nodes[start_node].set_incidence(8,v_);

        if( scn_type == 1){

            v = WG_nodes[start_node].get_incidence(15);
            v_ = v_inc + v;
            WG_nodes[start_node].set_incidence(15,v_);

            v = WG_nodes[start_node].get_incidence(16);
            v_ = v_inc + v;
            WG_nodes[start_node].set_incidence(16,v_);

        }
    }
    else if(axis == 3)
    {
        v = WG_nodes[start_node].get_incidence(7);
        v_ = v_inc + v;
        WG_nodes[start_node].set_incidence(7,v_);

        v = WG_nodes[start_node].get_incidence(5);
        v_ = v_inc + v;
        WG_nodes[start_node].set_incidence(5,v_);

        v = WG_nodes[start_node].get_incidence(6);
        v_= v_inc + v;
        WG_nodes[start_node].set_incidence(6,v_);

        v = WG_nodes[start_node].get_incidence(10);
        v_= v_inc + v;
        WG_nodes[start_node].set_incidence(10,v_);

            if( scn_type == 1){

                v = WG_nodes[start_node].get_incidence(17);
                v_ = v_inc + v;
                WG_nodes[start_node].set_incidence(17,v_);

                v = WG_nodes[start_node].get_incidence(18);
                v_ = v_inc + v;
                WG_nodes[start_node].set_incidence(18,v_);

            }
        }
}

void WG_16::TLM_excitation_plane_node(int node, double v_inc)
{
   int nxy = Nyy*Nxx;
   double v = 0, v_ = 0;
   int start_node = node;

    for ( int j = 0; j < nxy ; j++){
            //EXCITE EY field component
            if( !WG_nodes[start_node].is_PEC())
            {
            v = WG_nodes[start_node].get_incidence(3);
            v_ = v_inc + v;
            WG_nodes[start_node].set_incidence(3,v_);

            v = WG_nodes[start_node].get_incidence(11);
            v_ = v_inc + v;
            WG_nodes[start_node].set_incidence(11,v_);

            v = WG_nodes[start_node].get_incidence(4);
            v_= v_inc + v;
            WG_nodes[start_node].set_incidence(4,v_);

            v = WG_nodes[start_node].get_incidence(8);
            v_= v_inc + v;
            WG_nodes[start_node].set_incidence(8,v_);

            if( scn_type == 1){

                v = WG_nodes[start_node].get_incidence(15);
                v_ = v_inc + v;
                WG_nodes[start_node].set_incidence(15,v_);

                v = WG_nodes[start_node].get_incidence(16);
                v_ = v_inc + v;
                WG_nodes[start_node].set_incidence(16,v_);

            }
            //cout<<"plane excitation "<<endl;
            start_node = start_node + 1;

        }
        else cout<< " CAUTION: PROBLEM!! EXCITING PEC !!" <<endl;
    }
}

void WG_16::TLM_excitation_by_imposing_inc(int start_node_z, vector <double> scn_inc )
{
   int nxy = Nyy*Nxx;
   int temp;
   double v = 0, v_ = 0, v_inc = 0, element = 0;
    //int start_node = node;
    //cout<<" size "<< scn_inc.size() <<endl;
    //cout<< scn_inc[0] << " "<< scn_inc[1];

    int start_node = 0 + start_node_z*nxy;
    start_node = return_coordinates_WG(WG_width,WG_height,WG_length,npml,WG_dl,0,0,2);

    for (int j = start_node; j < start_node+nxy ;j++)
    {

        for (int i = 1; i <= 12; i++)
        {
            //if(( i==10)||(i==11))
           //{
            //cout<<scn_inc.size()<<endl;cin>>j;
            v_inc = scn_inc[element];             // port voltage
            v = WG_nodes[j].get_incidence(i);
            v_ =  v_inc + v;

            WG_nodes[j].set_incidence(i,v_);
            //}
            element = element + 1;

          // cout<< " id " << element-1 << " incident at " << i<<" "<< v_ << endl;
        }
            //cout<< " " << WG_nodes[j].get_coord(1)<< "  " ;
            //if( WG_nodes[j].get_coord(0)== 19 ) cout<<endl<<endl;
        //}
       // cin>>temp;
        //if(j == start_node+nxy-1) cout<< WG_nodes[j].get_incidence(12)<<endl;

    }

 //cin>>temp;
}

void WG_16::TLM_simulation1(int dt_total,double tfactor, int output_node,int excited_node,
                                         double Zz1, double Zz2, double Zy1, double Zy2,double Zx1, double Zx2,
                                                    f_gaussian func, double exctn_type, int srce_res, bool scn_xy)
{

    int temp = 0;

    int nxy = Nyy*Nxx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;
    double Isrc = 0,  Vy=0,Vy_prime=0;;

// Guassian Pulse Parameters
    double f = 0;           //frequency
    double ff = 0;          //dummy variable
    double te = 0;          //dummy variable
    double t = tfactor*(1e-3)*dl/(c*2);     //time factor

// Compute the start node in desired excitation plane
    int plane_start_node = excited_node, start_nodee = 0;
 // Used in exciting TE10 mode
    int te10_nd = excited_node;



// Exciting SCN ports in xy.
     vector< vector <double> > v_inc ;//= parseIncident_2File("V_incident_xy.csv");
     //v_inc =  parseIncident_2File("V_incident_wg_xy.csv");
     int max_time = 1;

     max_time = v_inc.size()-1;

//plane excitation
  /*  if( exctn_type != -1 )
    {
           // start_node = return_coordinates ( w,h,l,npml,dl,0,0,exctn_type ); waveguide
           start_node = Nxx*Nyy*(npml+exctn_type) + Nxx*npml + npml;
    }

//return_coordinates ( w,h,l,npml,dl,0,0,dl);        // exciting from one node from the interface
//output_node = return_coordinates( w,h,l,npml,dl,3,2,dl);


//3D cubic domain....setting the excited node to be the centre node
    if( this->sim_type == 1)
    {
        excited_node = this->centre_node;
        output_node = excited_node;

    }
 */
 cout<<" here "<<endl;
    cout<<" source_Resistance "<< srce_res <<endl;
    cout<<" excitation type "<< exctn_type <<endl;
    cout<<" excitation at "<<excited_node<<" output at "<< output_node<<endl;

// setting the soft source node
    if( srce_res>0) this->WG_nodes[excited_node].set_soft_source_node();

// TIME LOOP
    for( int dt=0; dt<dt_total; dt++ )
    {
        //EXCITATION
        if( ( dt<900) )
        {
// single node excitation
            if ( exctn_type==-1 )
            {
                //cout<<" Point source Excitation  at " << excited_node << endl;
                if( !WG_nodes[excited_node].is_PEC())
                {
                   f = func(dt)*(-1e-3*dl);
                   //srand (time(NULL));

                    /* generate secret number between 1 and 10: */
                   //f = rand() % 10 + 1;

                   //f = (f/10)*-1e-3*dl;

                    output_node = return_coordinates_3D (w,h,l,npml,dl,6,4,8);
                    //output_node = return_coordinates_3D (w,h,l,npml,dl,1,1,1);

                    if (srce_res==0) this->TLM_excitation_single_node(excited_node, f,2,srce_res);
                }
                else { cout<<"ERROR: Attempted to excite a PEC!!"<<endl; return ;}

            }
// TE10 excitation
            else if( exctn_type==-2 )
            {
            //cout<<" TE10 Excitation  at " << te10_nd << endl;
                int xy_end = te10_nd + Nxx;
                for (int i = te10_nd; i< xy_end ;i++ )
                {
                    f = sin((i-te10_nd)*pi/(Nxx-1))*func(dt)*1e-3*dl;
                    //f =  -1e-3*dl;
                    this->TLM_excitation_te10(i, f);
                }
            }
//plane excitation
            else if ( exctn_type ==0 )
            {
            //cout<<" Plane Excitation  at " << start_node << endl;
               f = func(dt)*(-1e-3)*dl;
               //f = 1;
               //cout<<dt<<endl;
                //f = -1e-3*dl;
               //cout<<f<<endl;
               this->TLM_excitation_plane_node(plane_start_node, f );
            }
//BOOT STRAP scn xy plane
            else if ( exctn_type == -3 )
            {
              //cout<<" in boot strap "<<endl;
               vector < double > tes;

               if (dt < max_time )
               {
                   tes = v_inc[dt+1];
                   this->TLM_excitation_by_imposing_inc(WG_length/2,tes);
               }
            }
                //int temp; cin>>temp;
        }


//cout<< WG_nodes[output_node].get_coord(0) << " " <<WG_nodes[output_node].get_coord(1)<< " " << WG_nodes[output_node].get_coord(2)<<endl;
//cout<< WG_nodes[excited_node].get_coord(0) << " " <<WG_nodes[excited_node].get_coord(1)<< " " << WG_nodes[excited_node].get_coord(2)<<endl;

//cin>>dt;
// soft source node
//When a soft source is applied, the excitation node is treated separately from the previous nodes
//scattering is done separately and output
// connection is unchanged!!
//future CODE modifications would need to address ways of improving code efficiency

//cout<<endl<<" time step " << dt<<endl<<endl;

        if( srce_res>0)  // soft source
        {
            f = func(dt)*(-1e-3*dl);
            //f=-1;
            this->WG_nodes[excited_node].scattering_source_node(f,srce_res, Vy_prime); // scatter is carried out here- for ease of code

            Isrc = (Vy_prime - f)/srce_res;

            if( excited_node == output_node ) WG_Ey.push_back(Vy_prime/(-1e-3*dl)); // assumes output node is excitation node
            else WG_Ey.push_back(Ey_output_at_node(output_node));
        }

//OUTPUT FIELDS
        //if( WG_nodes[output_node].is_PEC()) cout<< "Output node is a PEC !! " <<endl;


        //cout<< dt <<endl;

        Is.push_back(Isrc);

        Vs.push_back(f);

        if( WG_nodes[output_node].is_PEC() ) { cout<<" error "<<endl;cout<<output_node<<endl; cin>>output_node;}

        WG_Ex.push_back(Ex_output_at_node(output_node));

        if ( srce_res == 0 ) WG_Ey.push_back(Ey_output_at_node(output_node));  // i.e.  hard source

        if ( excited_node != output_node) Ey_outputs.push_back(Ey_output_at_node(excited_node));

                   //cout<< dt << " .... " <<Ey_output_at_node(output_node+5*Nxx)<<endl;

        WG_Ez.push_back(Ez_output_at_node(output_node));

        WG_Hx.push_back(Hx_output_at_node(output_node));

        WG_Hy.push_back(Hy_output_at_node(output_node));

        WG_Hz.push_back(Hz_output_at_node(output_node));

 //gets the index of output.

 //cout<< "output node iD: " << WG_nodes[output_node].get_iD() <<endl;

        int x_plane_id = WG_nodes[output_node].get_coord(0);

        int z_plane_id = WG_nodes[output_node].get_coord(2);

        int y_plane_id = WG_nodes[output_node].get_coord(1);

//cout<< x_plane_id << " " << y_plane_id << " " << z_plane_id <<endl<<endl;

        average_Ey.push_back( this->compute_field_average_in_z_plane(z_plane_id) );

        Ey_field_along_line.push_back( this->compute_field_along_line( x_plane_id, y_plane_id,z_plane_id ));
 //cin>>dt;
/* CODE FOR STORING THE INCIDENT PULSES ACROSS THE XY PLANE - COULD ALSO BE APPLIED TO STORING THE FIELD VALUES

        // vector<double> Ey_v_plane_xy, Ex_v_plane_xy, Ez_v_plane_xy,Hy_v_plane_xy, Hx_v_plane_xy, Hz_v_plane_xy,
        vector<double> v_inc_ports;
        int fss_node = return_coordinates_WG (WG_width,WG_height,WG_length,npml,WG_dl,0,0,WG_length/2);
        for( int ii = 0; ii<Nxx*Nyy ; ii++)
            {
                // Ey_v_plane_xy.push_back(Ey_output_at_node(ii+fss_node));
                vector<double> temp;
                temp = WG_nodes[ii+fss_node].get_SCN_inc_ports();
                temp.resize(12);
                v_inc_ports.insert(v_inc_ports.end(),temp.begin(),temp.end());
            }
        V_incident_xy.push_back(v_inc_ports);
        //  WG_Ey_plane_xy.push_back(v_plane_xy);
*/

//OUTPUTING FIELD VALUES ACROSS A PLANE XZ or XY or YZ
/*
     if (( dt % 5) == 0)
     {
        vector<double> v_plane_xz;
        int start_node_xz = 0 + Nyy/2 * Nxx;
        for( int i = start_node_xz; i<Ntotal_ ; i = i + Nxx*Nyy )
        {
            for( int j = 0; j<Nxx; j++)
            {
                    v_plane_xz.push_back(Ey_output_at_node(i+j));
                    //cout<< i+j << " ";
            }
           // cout<<endl;
        }
        //cin>>dt;
        WG_Ey_plane_zx.push_back(v_plane_xz);
    }
*/
    //Far Field Computation variables
    //vector<SCN_node> v_plane_xy1,v_plane_xy2,v_plane_yz1,v_plane_yz2,v_plane_zx1,v_plane_zx2;

    bool check = false, ch = false;

    //cout<< " time step " << dt <<endl;
    //cin>>ch;


    //SCATTERING PROCESS
    if( scn_type ==1 )
    {
    #pragma omp parallel
    #pragma omp for
        for(int i=0; i< Ntotal_; i++)
                    if(!WG_nodes[i].is_PEC()) WG_nodes[i].scattering_PML() ;
    }

    else if ( scn_type == 2)
    {
    #pragma omp parallel
    #pragma omp for
        for(int i=0; i< Ntotal_; i++)
                   {
                      //cout<<" iii "<<i <<endl;
                     if(!WG_nodes[i].is_PEC()) WG_nodes[i].scattering_PML_df(check) ;
                   }

    }

    else if ( scn_type == 3 )
    {
    #pragma omp parallel
    #pragma omp for
        for(int i=0; i< Ntotal_; i++)
                    if(!WG_nodes[i].is_PEC()) WG_nodes[i].scattering_EPML() ;

    }

   if( sim_type == 0 )
    {
        if (Ny > 1 )this->TLM_connection_optimized(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);
        else if( Ny==1) this->TLM_connection_2D_SCN_wg(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);  // strictly for the 2D SCN waveguide
    }
    else if( sim_type ==1 ) this->TLM_connection_optimized(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);// this->TLM_connection_optimized(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);
    }

}

void WG_16::TLM_simulation_2D_SCN(int dt_total,double tfactor, int output_node,int excited_node,double Zz1, double Zz2, double Zy1, double Zy2,double Zx1, double Zx2, f_gaussian func, double z)
{

    int nxy = Nyy*Nxx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    output_node = centre_node;
    excited_node  = centre_node;

// Guassian Pulse Parameters
    double f = 0;           //frequency
    double t = tfactor*(1e-3)*dl/(c*2);     //time factor

// Compute the start node in desired excitation plan

// TIME LOOP
    for( int dt=0; dt<dt_total; dt++ )
    {
        //EXCITATION
        if( (dt<1000) )
        {
            // single node excitation
                if( !WG_nodes[excited_node].is_PEC())
                {

                    f = func(dt)*(-1e-3*dl);

                    this->TLM_excitation_single_node(centre_node, f,3);

                }
                else {
                cout<<"ERROR: Attempted to excite a PEC!!"<<endl;
                return ;
                }

        }

        //OUTPUT FIELDS

        WG_Ex.push_back(Ex_output_at_node(output_node));

        WG_Ey.push_back(Ey_output_at_node(output_node));

        WG_Ez.push_back(Ez_output_at_node(output_node));

        WG_Hx.push_back(Hx_output_at_node(output_node));

        WG_Hy.push_back(Hy_output_at_node(output_node));

        WG_Hz.push_back(Hz_output_at_node(output_node));

//cout<<Ez_output_at_node(output_node)<<endl;
        //Far Field Computation variables
        //vector<SCN_node> v_plane_xy1,v_plane_xy2;

        bool check = false;

    //SCATTERING PROCESS
    if( scn_type ==1)
    {
    #pragma omp parallel
    #pragma omp for

        for(int i=0; i< Ntotal_; i++)
                    if(!WG_nodes[i].is_PEC()) WG_nodes[i].scattering_PML() ;
    }


    else if (scn_type == 2)
    {
    #pragma omp parallel
    #pragma omp for
        for(int i=0; i< Ntotal_; i++)
                    if(!WG_nodes[i].is_PEC()) WG_nodes[i].scattering_PML_df(check) ;


    }

    this->TLM_connection_2D_SCN(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);

    }
}

void WG_16::TLM_connection(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 )
{
   // Variable Definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml;
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    // Boundary Conditions of Computational Domain
    double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
    double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
    double reflct_x = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
    double reflct_y = (Zy1 - Z0) /( Zy1 + Z0);

    if( Zz1 < 0) reflct_z1 = 1;
    if( Zz2 < 0) reflct_z2 = 1;
    if( Zy1 < 0) reflct_y  = 1;
    if( Zx1 < 0) reflct_x  = 1;

    //CONNECTION PROCESS
    int ii = 0;

    for(int z=0 ; z<this->Nzz; z++)
        {
            for (int y=0 ; y<this->Nyy ; y++)
            {
                for (int x=0 ; x<this->Nxx ; x++)
                {
                    if(!WG_nodes[ii].check_special_node()){

                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));

                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));

                        WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                        WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                        WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                        WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));

                    }

                    else {
                    if ( (x==0)||(x==nnx-1) )
                    {

                            if ( (x==0) )//left wall rule.
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
                    else
                    {
                        // reflect into adjacent nodes
                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                        //cout<<" Position : "<<" 5 "<<endl;
                    }

                    if (( y==0) || (y==nny-1))
                    {
                            if (y==0)   //bottom wall rule
                            {
                                WG_nodes[ii].set_incidence(1, reflct_y*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, reflct_y*WG_nodes[ii].get_reflected(5));
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                                //cout<<" Position : "<<" 9 "<<endl;
                            }
                            else       //top wall rule
                            {
                                WG_nodes[ii].set_incidence(12,reflct_y*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, reflct_y*WG_nodes[ii].get_reflected(7));
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                            }
                    }
                    else       //reflect into adjacent nodes
                    {
                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                        //cout<<" Position : "<<"12"<<endl;
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
                        //cout<<" Position : "<<"17"<<endl;
                    }
                }
                    //STUB CONNECTION
                        WG_nodes[ii].set_incidence(13,WG_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        WG_nodes[ii].set_incidence(14,WG_nodes[ii].get_reflected(14));
                        WG_nodes[ii].set_incidence(15,WG_nodes[ii].get_reflected(15));
                        WG_nodes[ii].set_incidence(16,WG_nodes[ii].get_reflected(16));
                        WG_nodes[ii].set_incidence(17,WG_nodes[ii].get_reflected(17));
                        WG_nodes[ii].set_incidence(18,WG_nodes[ii].get_reflected(18));
                        WG_nodes[ii].set_incidence(19,-1*WG_nodes[ii].get_reflected(19));    //inductive stubs become  short circuits
                        WG_nodes[ii].set_incidence(20,-1*WG_nodes[ii].get_reflected(20));
                        WG_nodes[ii].set_incidence(21,-1*WG_nodes[ii].get_reflected(21));
                        WG_nodes[ii].set_incidence(22,-1*WG_nodes[ii].get_reflected(22));
                        WG_nodes[ii].set_incidence(23,-1*WG_nodes[ii].get_reflected(23));
                        WG_nodes[ii].set_incidence(24,-1*WG_nodes[ii].get_reflected(24));
                     ii=ii+1;
                    }
                }
            }
}


void WG_16::TLM_connection_2D_SCN_wg(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 )
{
     if( Ny >1)
     {
        int a;
        cout<<" Error!! calling 2D SCN connection for 3D SCN node" <<endl;
        cin>>a;
     }

     // Variable Definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml;
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    // Boundary Conditions of Computational Domain
    double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
    double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
    double reflct_x = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
    double reflct_y = (Zy1 - Z0) /( Zy1 + Z0);

    if( Zz1 < 0) reflct_z1 = 1;
    if( Zz2 < 0) reflct_z2 = 1;
    if( Zy1 < 0) reflct_y  = 1;
    if( Zx1 < 0) reflct_x  = 1;

    //CONNECTION PROCESS

    int ii =0;

    for(int z=0 ; z<this->Nzz; z++)
        {
            for (int x=0 ; x<this->Nxx ; x++)
                {

                    if( !WG_nodes[ii].is_PEC() )
                    {

                        if(!WG_nodes[ii].check_special_node())
                        {

                            WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                            WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                            WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                            WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));

                            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));

                        }

                    else
                    {

                        if ( (x==0)||(x==nnx-1) )
                        {
                            if ( (x==0) )//left wall rule.
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

                        else if ((WG_nodes[ii+1].is_PEC())||(WG_nodes[ii-1].is_PEC()) )
                        {
                            if (WG_nodes[ii-1].is_PEC()) //left wall rule.
                            {
                                WG_nodes[ii].set_incidence(6, -1*(WG_nodes[ii].get_reflected(6)));
                                WG_nodes[ii].set_incidence(3, -1*(WG_nodes[ii].get_reflected(3)));
                                if (!WG_nodes[ii+1].is_PEC())
                                {
                                WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                                WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                               // cout<<" Position : "<<" 6a"<<endl;
                                }
                                else
                                {
                                WG_nodes[ii].set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
                                //cout<<" Position : "<<" 6 "<<endl;
                                }
                            }
                            else    //right wall rule
                            {
                                WG_nodes[ii].set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
                                WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                                WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                                //cout<<" Position : "<<" 6b"<<endl;
                            }
                        }
                        else
                        {
                            // reflect into adjacent nodes
                            WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                            WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                            WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                            WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                            //cout<<" Position : "<<" 5 "<<endl;
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

                        else if ((WG_nodes[ii+nxy].is_PEC())||(WG_nodes[ii-nxy].is_PEC()))
                        {

                            if ((WG_nodes[ii-nxy].is_PEC()))   // front face rule
                            {
                                WG_nodes[ii].set_incidence(8, -1*WG_nodes[ii].get_reflected(8));
                                WG_nodes[ii].set_incidence(9, -1*WG_nodes[ii].get_reflected(9));
                                if (!WG_nodes[ii+nxy].is_PEC())
                                {
                                WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                                WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                                //cout<<" Position : "<<"18a "<<endl;
                                }
                                else
                                {
                                WG_nodes[ii].set_incidence(4,-1*WG_nodes[ii].get_reflected(4));
                                WG_nodes[ii].set_incidence(2,-1*WG_nodes[ii].get_reflected(2));
                                //cout<<" Position : "<<"18"<<endl;
                                }
                            }
                            else       //back face rule
                            {
                                WG_nodes[ii].set_incidence(4,-1*WG_nodes[ii].get_reflected(4));
                                WG_nodes[ii].set_incidence(2,-1*WG_nodes[ii].get_reflected(2));
                                WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                                WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                                //cout<<" Position : "<<"18b"<<endl;
                            }
                        }
                        else        //normal reflection rules
                        {
                            // reflect into adjacent nodes
                            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                            //cout<<" Position : "<<"17"<<endl;
                        }
                        }
                     //Connection along y axis
                        WG_nodes[ii].set_incidence(12,-WG_nodes[ii].get_reflected(12));
                        WG_nodes[ii].set_incidence(1, -WG_nodes[ii].get_reflected(1));
                        WG_nodes[ii].set_incidence(7, -WG_nodes[ii].get_reflected(7));
                        WG_nodes[ii].set_incidence(5, -WG_nodes[ii].get_reflected(5));

                    //STUB CONNECTION
                        WG_nodes[ii].set_incidence(13,WG_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        WG_nodes[ii].set_incidence(14,WG_nodes[ii].get_reflected(14));
                        WG_nodes[ii].set_incidence(15,WG_nodes[ii].get_reflected(15));
                        WG_nodes[ii].set_incidence(16,WG_nodes[ii].get_reflected(16));
                        WG_nodes[ii].set_incidence(17,WG_nodes[ii].get_reflected(17));
                        WG_nodes[ii].set_incidence(18,WG_nodes[ii].get_reflected(18));
                        WG_nodes[ii].set_incidence(19,-1*WG_nodes[ii].get_reflected(19));    //inductive stubs become  short circuits
                        WG_nodes[ii].set_incidence(20,-1*WG_nodes[ii].get_reflected(20));
                        WG_nodes[ii].set_incidence(21,-1*WG_nodes[ii].get_reflected(21));
                        WG_nodes[ii].set_incidence(22,-1*WG_nodes[ii].get_reflected(22));
                        WG_nodes[ii].set_incidence(23,-1*WG_nodes[ii].get_reflected(23));
                        WG_nodes[ii].set_incidence(24,-1*WG_nodes[ii].get_reflected(24));
                    }
                    ii=ii+1;
                }
        }
        //cout<<"here"<<endl;
}

void WG_16::TLM_connection_2D_SCN(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 )
{
    // Variable Definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml;
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    // Boundary Conditions of Computational Domain
    double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
    double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
    double reflct_x  = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
    double reflct_y  = (Zy1 - Z0) /( Zy1 + Z0);

    if( Zz1 < 0) reflct_z1 = 1;
    if( Zz2 < 0) reflct_z2 = 1;
    if( Zy1 < 0) reflct_y  = 1;
    if( Zx1 < 0) reflct_x  = 1;

    //CONNECTION PROCESS
    int ii = 0;


            for (int y=0 ; y<this->Nyy ; y++)
            {
                for (int x=0 ; x<this->Nxx ; x++)
                {

                    if( !WG_nodes[ii].is_PEC() ) {

                    if(!WG_nodes[ii].check_special_node()){

                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));

                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));

                    }
                    else{

                    if ( (x==0)||(x==nnx-1) )
                    {

                            if ( (x==0) )//left wall rule.
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

                    else if ((WG_nodes[ii+1].is_PEC())||(WG_nodes[ii-1].is_PEC()) )
                    {
                            if (WG_nodes[ii-1].is_PEC()) //left wall rule.
                            {
                                WG_nodes[ii].set_incidence(6, -1*(WG_nodes[ii].get_reflected(6)));
                                WG_nodes[ii].set_incidence(3, -1*(WG_nodes[ii].get_reflected(3)));
                                if (!WG_nodes[ii+1].is_PEC())
                                {
                                WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                                WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                               // cout<<" Position : "<<" 6a"<<endl;
                                }
                                else
                                {
                                WG_nodes[ii].set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
                                //cout<<" Position : "<<" 6 "<<endl;
                                }
                            }
                            else    //right wall rule
                            {
                                WG_nodes[ii].set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
                                WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                                WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                                //cout<<" Position : "<<" 6b"<<endl;
                            }
                    }
                    else
                    {
                        // reflect into adjacent nodes
                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                        //cout<<" Position : "<<" 5 "<<endl;
                    }

                    if (( y==0) || (y==nny-1))
                    {
                            if (y==0)   //bottom wall rule
                            {
                                WG_nodes[ii].set_incidence(1, reflct_y*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, reflct_y*WG_nodes[ii].get_reflected(5));
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                                //cout<<" Position : "<<" 9 "<<endl;
                            }
                            else       //top wall rule
                            {
                                WG_nodes[ii].set_incidence(12,reflct_y*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, reflct_y*WG_nodes[ii].get_reflected(7));
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                            }
                    }

                    else if ((WG_nodes[ii+nnx].is_PEC())||(WG_nodes[ii-nnx].is_PEC()) )
                    {
                       if (WG_nodes[ii-nnx].is_PEC())   //bottom wall rule
                            {
                                WG_nodes[ii].set_incidence(1, -1*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, -1*WG_nodes[ii].get_reflected(5));
                                if (!WG_nodes[ii+nnx].is_PEC())
                                {
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                                }
                                else
                                {
                                WG_nodes[ii].set_incidence(12,-1*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, -1*WG_nodes[ii].get_reflected(7));
                                }
                            }
                        else        //top wall rule
                            {
                                WG_nodes[ii].set_incidence(12,-1*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, -1*WG_nodes[ii].get_reflected(7));
                                if (!WG_nodes[ii-nnx].is_PEC())
                                {
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                                //cout<<" Position : "<<" 11b"<<endl;
                                }
                                else
                                {
                                WG_nodes[ii].set_incidence(1, -1*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, -1*WG_nodes[ii].get_reflected(5));
                                //cout<<" Position : "<<"11"<<endl;

                                }
                            }
                    }

                    else       //reflect into adjacent nodes
                    {
                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                        //cout<<" Position : "<<"12"<<endl;
                    }

                    }
                    //connection in z axis
                        WG_nodes[ii].set_incidence(8, -1*WG_nodes[ii].get_reflected(8));
                        WG_nodes[ii].set_incidence(9, -1*WG_nodes[ii].get_reflected(9));
                        WG_nodes[ii].set_incidence(4, -1*WG_nodes[ii].get_reflected(4));
                        WG_nodes[ii].set_incidence(2, -1*WG_nodes[ii].get_reflected(2));

                    //STUB CONNECTION
                        WG_nodes[ii].set_incidence(13,WG_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        WG_nodes[ii].set_incidence(14,WG_nodes[ii].get_reflected(14));
                        WG_nodes[ii].set_incidence(15,WG_nodes[ii].get_reflected(15));
                        WG_nodes[ii].set_incidence(16,WG_nodes[ii].get_reflected(16));
                        WG_nodes[ii].set_incidence(17,WG_nodes[ii].get_reflected(17));
                        WG_nodes[ii].set_incidence(18,WG_nodes[ii].get_reflected(18));
                        WG_nodes[ii].set_incidence(19,-1*WG_nodes[ii].get_reflected(19));    //inductive stubs become  short circuits
                        WG_nodes[ii].set_incidence(20,-1*WG_nodes[ii].get_reflected(20));
                        WG_nodes[ii].set_incidence(21,-1*WG_nodes[ii].get_reflected(21));
                        WG_nodes[ii].set_incidence(22,-1*WG_nodes[ii].get_reflected(22));
                        WG_nodes[ii].set_incidence(23,-1*WG_nodes[ii].get_reflected(23));
                        WG_nodes[ii].set_incidence(24,-1*WG_nodes[ii].get_reflected(24));

                    }
                     ii=ii+1;

                }
            }

}

void WG_16::TLM_connection_sd(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 )
{

   // variable definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = 0.5*(nnzz - nnz);
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    // Boundary conditions of computational domain
    double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
    double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
    double reflct_x = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
    double reflct_y = (Zy1 - Z0) /( Zy1 + Z0);

    if( Zz1 < 0) reflct_z1 = 1;
    if( Zz2 < 0) reflct_z2 = 1;
    if( Zy1 < 0) reflct_y  = 1;
    if( Zx1 < 0) reflct_x  = 1;

    //double output_node = return_coordinates_3D( WG_width,WG_height,WG_length,npml,WG_dl,10*dl,10*dl,28*dl);


    //CONNECTION PROCESS
    int ii = 0;
   // cout<<" CONNECTION : "<<endl;

    for(int z=0 ; z<this->Nzz; z++)
        {
            for (int y=0 ; y<this->Nyy ; y++)
            {
                for (int x=0 ; x<this->Nxx ; x++)
                {
                     // cout<<"here : "<<endl;
                    if( !WG_nodes[ii].is_PEC() ) {

                    if(!WG_nodes[ii].check_special_node())
                    {

                       // if( ii == output_node ) cout<<"node id "<< ii<<" .... is PEC ? .. "<<WG_nodes[output_node].is_PEC()<<endl;
                        //Connect
                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));

                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));

                        WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                        WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                        WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                        WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));

                    }
                    else{

                    if ( (x==0)||(x==nnx-1) )
                    {

                            if ( (x==0) )//left wall rule.
                            {
                                //Connect B
                                WG_nodes[ii].set_incidence(6, reflct_x*(WG_nodes[ii].get_reflected(6)));
                                WG_nodes[ii].set_incidence(3, reflct_x*(WG_nodes[ii].get_reflected(3)));
                                WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                                WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                            }
                            else    //right wall rule
                            {
                                //Connect C
                                WG_nodes[ii].set_incidence(11,reflct_x*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,reflct_x*(WG_nodes[ii].get_reflected(10)));
                                WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                                WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                            }
                    }
                    else if ((WG_nodes[ii+1].is_PEC())||(WG_nodes[ii-1].is_PEC()) )
                    {
                            if (WG_nodes[ii-1].is_PEC()) //left wall rule.
                            {
                                //Connect
                                WG_nodes[ii].set_incidence(6, -1*(WG_nodes[ii].get_reflected(6)));
                                WG_nodes[ii].set_incidence(3, -1*(WG_nodes[ii].get_reflected(3)));
                                if (!WG_nodes[ii+1].is_PEC())
                                {//Connect D
                                WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                                WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                               // cout<<" Position : "<<" 6a"<<endl;
                                }
                                else
                                {
                                //Connect E
                                WG_nodes[ii].set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
                                //cout<<" Position : "<<" 6 "<<endl;
                                }
                            }
                            else    //right wall rule
                            {
                                //Connect F
                                WG_nodes[ii].set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
                                WG_nodes[ii].set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
                                WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                                WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                                //cout<<" Position : "<<" 6b"<<endl;
                            }
                    }
                    else
                    {
                        // reflect into adjacent nodes
                        //Connect A
                        WG_nodes[ii].set_incidence(11,WG_nodes[ii+1].get_reflected(3));
                        WG_nodes[ii].set_incidence(3,WG_nodes[ii-1].get_reflected(11));
                        WG_nodes[ii].set_incidence(10,WG_nodes[ii+1].get_reflected(6));
                        WG_nodes[ii].set_incidence(6,WG_nodes[ii-1].get_reflected(10));
                        //cout<<" Position : "<<" 5 "<<endl;
                    }

                    if (( y==0) || (y==nny-1))
                    {
                            if (y==0)   //bottom wall rule
                            {//Connect G
                                WG_nodes[ii].set_incidence(1, reflct_y*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, reflct_y*WG_nodes[ii].get_reflected(5));
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                                //cout<<" Position : "<<" 9 "<<endl;
                            }
                            else       //top wall rule
                            {//Connect H
                                WG_nodes[ii].set_incidence(12,reflct_y*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, reflct_y*WG_nodes[ii].get_reflected(7));
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                            }
                    }

                    else if ((WG_nodes[ii+nnx].is_PEC())||(WG_nodes[ii-nnx].is_PEC()) )
                    {
                       if (WG_nodes[ii-nnx].is_PEC())   //bottom wall rule
                            {//
                                WG_nodes[ii].set_incidence(1, -1*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, -1*WG_nodes[ii].get_reflected(5));
                                if (!WG_nodes[ii+nnx].is_PEC())
                                {//Connect I
                                WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                                WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                                }
                                else
                                {//Connect J
                                WG_nodes[ii].set_incidence(12,-1*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, -1*WG_nodes[ii].get_reflected(7));
                                }
                            }
                        else        //top wall rule
                            {

                                WG_nodes[ii].set_incidence(12,-1*WG_nodes[ii].get_reflected(12));
                                WG_nodes[ii].set_incidence(7, -1*WG_nodes[ii].get_reflected(7));
                                if (!WG_nodes[ii-nnx].is_PEC())
                                {//Connect K
                                WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                                WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                                //cout<<" Position : "<<" 11b"<<endl;
                                }
                                else
                                {
                                WG_nodes[ii].set_incidence(1, -1*WG_nodes[ii].get_reflected(1));
                                WG_nodes[ii].set_incidence(5, -1*WG_nodes[ii].get_reflected(5));
                                //cout<<" Position : "<<"11"<<endl;

                                }
                            }
                    }

                    else       //reflect into adjacent nodes
                    {//Connect L
                        WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
                        WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
                        WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
                        WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
                        //cout<<" Position : "<<"12"<<endl;
                    }

                    if ( (z==0) || (z==nnzz-1))
                    {
                        if (z==0)   // front face rule
                        {//Connect M
                            WG_nodes[ii].set_incidence(8, reflct_z1*WG_nodes[ii].get_reflected(8));
                            WG_nodes[ii].set_incidence(9, reflct_z1*WG_nodes[ii].get_reflected(9));
                            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                        }
                        else        //back face rule
                        {//Connect N
                            WG_nodes[ii].set_incidence(4,reflct_z2*WG_nodes[ii].get_reflected(4));
                            WG_nodes[ii].set_incidence(2,reflct_z2*WG_nodes[ii].get_reflected(2));
                            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                        }
                    }
                    else if ((WG_nodes[ii+nxy].is_PEC())||(WG_nodes[ii-nxy].is_PEC()))
                    {
                        if ((WG_nodes[ii-nxy].is_PEC()))   // front face rule
                        {
                            WG_nodes[ii].set_incidence(8, -1*WG_nodes[ii].get_reflected(8));
                            WG_nodes[ii].set_incidence(9, -1*WG_nodes[ii].get_reflected(9));
                            if (!WG_nodes[ii+nxy].is_PEC())
                            {//Connect O
                            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                            //cout<<" Position : "<<"18a "<<endl;
                            }
                            else
                            {//Connect P
                            WG_nodes[ii].set_incidence(4,-1*WG_nodes[ii].get_reflected(4));
                            WG_nodes[ii].set_incidence(2,-1*WG_nodes[ii].get_reflected(2));
                            //cout<<" Position : "<<"18"<<endl;
                            }
                        }
                        else       //back face rule
                        {//Connect Q
                            WG_nodes[ii].set_incidence(4,-1*WG_nodes[ii].get_reflected(4));
                            WG_nodes[ii].set_incidence(2,-1*WG_nodes[ii].get_reflected(2));
                            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                            //cout<<" Position : "<<"18b"<<endl;
                        }
                    }
                    else        //normal reflection rules
                    {
                        // reflect into adjacent nodes
                        //Connect R
                        WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
                        WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
                        WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
                        WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
                        //cout<<" Position : "<<"17"<<endl;

                    }
                }
                    //STUB CONNECTION
                    //set_incidence(12,get_reflected(12));
                        WG_nodes[ii].set_incidence(13,WG_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        WG_nodes[ii].set_incidence(14,WG_nodes[ii].get_reflected(14));
                        WG_nodes[ii].set_incidence(15,WG_nodes[ii].get_reflected(15));
                        WG_nodes[ii].set_incidence(16,WG_nodes[ii].get_reflected(16));
                        WG_nodes[ii].set_incidence(17,WG_nodes[ii].get_reflected(17));
                        WG_nodes[ii].set_incidence(18,WG_nodes[ii].get_reflected(18));
                        WG_nodes[ii].set_incidence(19,-1*WG_nodes[ii].get_reflected(19));    //inductive stubs become  short circuits
                        WG_nodes[ii].set_incidence(20,-1*WG_nodes[ii].get_reflected(20));
                        WG_nodes[ii].set_incidence(21,-1*WG_nodes[ii].get_reflected(21));
                        WG_nodes[ii].set_incidence(22,-1*WG_nodes[ii].get_reflected(22));
                        WG_nodes[ii].set_incidence(23,-1*WG_nodes[ii].get_reflected(23));
                        WG_nodes[ii].set_incidence(24,-1*WG_nodes[ii].get_reflected(24));

                   // if(ii == centre_node) cout<<"...11.."<< WG_nodes[ii-1].get_reflected(11)<<"...3.. "<<WG_nodes[ii+1].get_reflected(3)<<
                      //  " ... "<<endl<<endl;
                    }

                    ii=ii+1;
                }
            }
        }
}

void WG_16::create_connect_bins()
{

     cout<<"Creating Bins for TLM Connect process ...."<<endl;
    // variable definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml; //0.5*(nnzz - nnz);
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    connect_bins.resize(19);  //creates 19 different connect categories
    int ii=0;


  for(int z=0 ; z<this->Nzz; z++)
        {
            for (int y=0 ; y<this->Nyy ; y++)
            {
                for (int x=0 ; x<this->Nxx ; x++)
                {
                     // cout<<"here : "<<endl;
                    if( !WG_nodes[ii].is_PEC() )
                    {
                        if(!WG_nodes[ii].check_special_node()) connect_bins[0].push_back(&WG_nodes[ii]);// STORE POINTER TO THE REGULAR NODES

                        else
                        {
                            if ( (x==0)||(x==nnx-1) )
                            {
                                if ( (x==0) ) connect_bins[2].push_back(&WG_nodes[ii]);//left wall rule.   //Connect B

                                else  connect_bins[3].push_back(&WG_nodes[ii]);  //right wall rule  //Connect C
                            }
                            else if ((WG_nodes[ii+1].is_PEC())||(WG_nodes[ii-1].is_PEC()) )
                            {
                                if (WG_nodes[ii-1].is_PEC()) //left wall rule.
                                   {
                                    if (!WG_nodes[ii+1].is_PEC()) connect_bins[4].push_back(&WG_nodes[ii]); //Connect D

                                    else connect_bins[5].push_back(&WG_nodes[ii]);//Connect E
                                   }
                                else connect_bins[6].push_back(&WG_nodes[ii]);// Connect F
                            }
                            else connect_bins[1].push_back(&WG_nodes[ii]);////Connect A



                            if (( y==0) || (y==nny-1))
                            {
                                if (y==0)connect_bins[7].push_back(&WG_nodes[ii]);//Connect G   //bottom wall rule

                                else  connect_bins[8].push_back(&WG_nodes[ii]);//Connect H      //top wall rule
                            }

                            else if ((WG_nodes[ii+nnx].is_PEC())||(WG_nodes[ii-nnx].is_PEC()) )
                            {
                               if (WG_nodes[ii-nnx].is_PEC())   //bottom wall rule
                                {
                                    if (!WG_nodes[ii+nnx].is_PEC()) connect_bins[9].push_back(&WG_nodes[ii]);//Connect I

                                    else connect_bins[10].push_back(&WG_nodes[ii]);//Connect J
                                }
                                else        //top wall rule
                                {
                                    if (!WG_nodes[ii-nnx].is_PEC())connect_bins[11].push_back(&WG_nodes[ii]);//Connect K
                                }
                            }
                            else  connect_bins[12].push_back(&WG_nodes[ii]);//Connect L     //reflect into adjacent nodes



                        if ( (z==0) || (z==nnzz-1))
                        {
                            if (z==0) connect_bins[13].push_back(&WG_nodes[ii]);//Connect M  // front face rule

                            else  connect_bins[14].push_back(&WG_nodes[ii]);//Connect N      //back face rule
                        }
                        else if ((WG_nodes[ii+nxy].is_PEC())||(WG_nodes[ii-nxy].is_PEC()))
                        {
                            if ((WG_nodes[ii-nxy].is_PEC()))   // front face rule
                            {
                                if (!WG_nodes[ii+nxy].is_PEC()) connect_bins[15].push_back(&WG_nodes[ii]);//Connect O

                                else connect_bins[16].push_back(&WG_nodes[ii]);//Connect P
                            }
                            else  connect_bins[17].push_back(&WG_nodes[ii]);//back face rule //Connect Q
                        }
                        else connect_bins[18].push_back(&WG_nodes[ii]); //Connect R//normal reflection rules
                    }

                }
                ii=ii+1;
            }
        }
    }
    cout<<" TLM Connect bins created......"<<endl;
}

void WG_16::TLM_connection_optimized(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 )
{
   // variable definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = 0.5*(nnzz - nnz);
    int nxy = nny*nnx;
    double dl = this->WG_dl;
    double w = this->WG_width;
    double h = this->WG_height;
    double l = this->WG_length;

    // Boundary conditions of computational domain
    double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
    double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
    double reflct_x = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
    double reflct_y = (Zy1 - Z0) /( Zy1 + Z0);

    if( Zz1 < 0) reflct_z1 = 1;
    if( Zz2 < 0) reflct_z2 = 1;
    if( Zy1 < 0) reflct_y  = 1;
    if( Zx1 < 0) reflct_x  = 1;

    //CONNECTION PROCESS
    int ii = 0;
   // cout<<" CONNECTION : "<<endl;
    int bin=0;

        //CONNECT FOR SPECIAL NODE
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,WG_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(3,WG_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(10,WG_nodes[ii+1].get_reflected(6));
            connect_bins[bin][i]->set_incidence(6,WG_nodes[ii-1].get_reflected(10));

            connect_bins[bin][i]->set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
            connect_bins[bin][i]->set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
            connect_bins[bin][i]->set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
            connect_bins[bin][i]->set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));

            connect_bins[bin][i]->set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
            connect_bins[bin][i]->set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
            connect_bins[bin][i]->set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
            connect_bins[bin][i]->set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
        }
        //CONNECT FOR 1
         bin = bin+1;
        #pragma omp parallel
        #pragma omp for
         for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,WG_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(3,WG_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(10,WG_nodes[ii+1].get_reflected(6));
            connect_bins[bin][i]->set_incidence(6,WG_nodes[ii-1].get_reflected(10));
        }

        //CONNECT FOR 2
         bin = bin+1;
         #pragma omp parallel
        #pragma omp for
         for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(6, reflct_x*(WG_nodes[ii].get_reflected(6)));
            connect_bins[bin][i]->set_incidence(3, reflct_x*(WG_nodes[ii].get_reflected(3)));
            connect_bins[bin][i]->set_incidence(11,WG_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(10,WG_nodes[ii+1].get_reflected(6));
        }
        //CONNECT FOR 3
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,reflct_x*(WG_nodes[ii].get_reflected(11)));
            connect_bins[bin][i]->set_incidence(10,reflct_x*(WG_nodes[ii].get_reflected(10)));
            connect_bins[bin][i]->set_incidence(3,WG_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(6,WG_nodes[ii-1].get_reflected(10));
        }

        //CONNECT FOR 4
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(6, -1*(WG_nodes[ii].get_reflected(6)));
            connect_bins[bin][i]->set_incidence(3, -1*(WG_nodes[ii].get_reflected(3)));
            connect_bins[bin][i]->set_incidence(11,WG_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(10,WG_nodes[ii+1].get_reflected(6));
        }

        //CONNECT FOR 5
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(6, -1*(WG_nodes[ii].get_reflected(6)));
            connect_bins[bin][i]->set_incidence(3, -1*(WG_nodes[ii].get_reflected(3)));
            connect_bins[bin][i]->set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
            connect_bins[bin][i]->set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
        }

        //CONNECT FOR 6
         bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,-1*(WG_nodes[ii].get_reflected(11)));
            connect_bins[bin][i]->set_incidence(10,-1*(WG_nodes[ii].get_reflected(10)));
            connect_bins[bin][i]->set_incidence(3,WG_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(6,WG_nodes[ii-1].get_reflected(10));
        }

        //CONNECT FOR 7
         bin = bin+1;
         #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(1, reflct_y*WG_nodes[ii].get_reflected(1));
            connect_bins[bin][i]->set_incidence(5, reflct_y*WG_nodes[ii].get_reflected(5));
            connect_bins[bin][i]->set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
            connect_bins[bin][i]->set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));

        }

        //CONNECT FOR 8
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(12,reflct_y*WG_nodes[ii].get_reflected(12));
            connect_bins[bin][i]->set_incidence(7, reflct_y*WG_nodes[ii].get_reflected(7));
            connect_bins[bin][i]->set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
            connect_bins[bin][i]->set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
        }
        //CONNECT FOR 9
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(1, -1*WG_nodes[ii].get_reflected(1));
            connect_bins[bin][i]->set_incidence(5, -1*WG_nodes[ii].get_reflected(5));
            connect_bins[bin][i]->set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
            connect_bins[bin][i]->set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
        }

        //CONNECT FOR 10
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(1, -1*WG_nodes[ii].get_reflected(1));
            connect_bins[bin][i]->set_incidence(5, -1*WG_nodes[ii].get_reflected(5));
            connect_bins[bin][i]->set_incidence(12,-1*WG_nodes[ii].get_reflected(12));
            connect_bins[bin][i]->set_incidence(7, -1*WG_nodes[ii].get_reflected(7));
        }
        //CONNECT FOR 11
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(12,-1*WG_nodes[ii].get_reflected(12)); // It is possible to do it this way
            WG_nodes[ii].set_incidence(7, -1*WG_nodes[ii].get_reflected(7));
            WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
            WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
        }

        //CONNECT FOR 12
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(12,WG_nodes[ii+nnx].get_reflected(1));
            WG_nodes[ii].set_incidence(1, WG_nodes[ii-nnx].get_reflected(12));
            WG_nodes[ii].set_incidence(7, WG_nodes[ii+nnx].get_reflected(5));
            WG_nodes[ii].set_incidence(5, WG_nodes[ii-nnx].get_reflected(7));
        }

        //CONNECT FOR 13
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(8, reflct_z1*WG_nodes[ii].get_reflected(8));
            WG_nodes[ii].set_incidence(9, reflct_z1*WG_nodes[ii].get_reflected(9));
            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
        }

        //CONNECT FOR 14
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(4,reflct_z2*WG_nodes[ii].get_reflected(4));
            WG_nodes[ii].set_incidence(2,reflct_z2*WG_nodes[ii].get_reflected(2));
            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
        }

        //CONNECT FOR 15
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(8, -1*WG_nodes[ii].get_reflected(8));
            WG_nodes[ii].set_incidence(9, -1*WG_nodes[ii].get_reflected(9));
            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
        }

        //CONNECT FOR 16
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(8, -1*WG_nodes[ii].get_reflected(8));
            WG_nodes[ii].set_incidence(9, -1*WG_nodes[ii].get_reflected(9));
            WG_nodes[ii].set_incidence(4,-1*WG_nodes[ii].get_reflected(4));
            WG_nodes[ii].set_incidence(2,-1*WG_nodes[ii].get_reflected(2));
        }

        //CONNECT FOR 17
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            WG_nodes[ii].set_incidence(4,-1*WG_nodes[ii].get_reflected(4));
            WG_nodes[ii].set_incidence(2,-1*WG_nodes[ii].get_reflected(2));
            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
        }

        //CONNECT FOR 18
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();
            //Connect N
            WG_nodes[ii].set_incidence(4, WG_nodes[ii+nxy].get_reflected(8));
            WG_nodes[ii].set_incidence(8, WG_nodes[ii-nxy].get_reflected(4));
            WG_nodes[ii].set_incidence(2, WG_nodes[ii+nxy].get_reflected(9));
            WG_nodes[ii].set_incidence(9, WG_nodes[ii-nxy].get_reflected(2));
        }

       if( bin!=18) cout<<"Error in optimized connect algorithm!!"<<endl;//cout<<bin<<endl;cin>>bin;
}


void WG_16::print_connect_bins()
{

    for(int i = 0;i<connect_bins.size();i++)
    {
        cout<<" Bin " <<i<<" : "<<connect_bins[i].size()<<endl;
    }
    //cout<<"centre_node  "<<centre_node<< "\\"<<connect_bins[10][0]->get_iD();
cout<<endl;
}

void WG_16:: write_output_file (const string& nodefileEX,const string& nodefileEY, const string& nodefileEZ, const string& nodefileHX,const string& nodefileHY,const string& nodefileHZ)                                           //writes to node file
{

    ofstream out_fileEY;
    ofstream out_fileEX;
    ofstream out_fileEZ;

    ofstream out_Is;
    ofstream out_Vs;

    ofstream out_fileHY;
    ofstream out_fileHX;
    ofstream out_fileHZ;

    ofstream out_averageEY;
    ofstream out_additional_ey;


    //string nodefileEX = "PEC_EX.txt";
    //string nodefileEZ = "PEC_EZ.txt";

    //string nodefileHX = "HX_TE10_pml_n2_91e105.txt"; //HX_TE10_pml_n0_91e16_ref
    //string nodefileHY = "PEC_HY.txt";
    //string nodefileHZ = "PEC_HZ.txt";

    out_fileEX.open((nodefileEX+"").c_str());
    out_fileEZ.open((nodefileEZ+"").c_str());
    out_fileEY.open((nodefileEY+"").c_str());

    out_fileHX.open((nodefileHX+"").c_str());
    out_fileHY.open((nodefileHY+"").c_str());
    out_fileHZ.open((nodefileHZ+"").c_str());

    string s = nodefileEY.substr(0,nodefileEY.length() - 4);

    out_Is.open((s+"_Is.txt").c_str());
    out_Vs.open((s+"_Vs.txt").c_str());
    out_averageEY.open((s+"_Average_EY.txt").c_str());


    int end__ = this->WG_Ex.size();
    int end_= this->WG_Ey.size();
    //cout<<end_<<endl;cin>>end_;

    out_fileEY<<setprecision(4)<<this->WG_dl<<endl;
    out_fileEX<<setprecision(4)<<this->WG_dl<<endl;
    out_fileEZ<<setprecision(4)<<this->WG_dl<<endl;
    out_fileHY<<setprecision(4)<<this->WG_dl<<endl;
    out_fileHX<<setprecision(4)<<this->WG_dl<<endl;
    out_fileHZ<<setprecision(4)<<this->WG_dl<<endl;
    out_Is<<setprecision(4)<<this->WG_dl<<endl;
    out_Vs<<setprecision(4)<<this->WG_dl<<endl;
    out_averageEY<<setprecision(4)<<this->WG_dl<<endl;


    out_fileEY<<setprecision(4)<<end_<<endl;
    out_fileEX<<setprecision(4)<<end_<<endl;
    out_fileEZ<<setprecision(4)<<end_<<endl;
    out_fileHY<<setprecision(4)<<end_<<endl;
    out_fileHX<<setprecision(4)<<end_<<endl;
    out_fileHZ<<setprecision(4)<<end_<<endl;
    out_Is<<setprecision(4)<<end_<<endl;
    out_Vs<<setprecision(4)<<end_<<endl;
    out_averageEY<<setprecision(4)<<end_<<endl;

    cout<< " Writing to file...." << endl;

    for(int ii = 0 ; ii< end_; ii++)
    {

         out_fileEY<<setprecision(numeric_limits<double>::digits10+1) <<WG_Ey[ii]<<endl;
         out_fileHX<<setprecision(numeric_limits<double>::digits10+1)<<WG_Hx[ii]<<endl;
         out_fileEX<<setprecision(numeric_limits<double>::digits10+1)<<WG_Ex[ii]<<endl;
         out_fileHY<<setprecision(numeric_limits<double>::digits10+1)<<WG_Hy[ii]<<endl;
         out_fileEZ<<setprecision(numeric_limits<double>::digits10+1)<<WG_Ez[ii]<<endl;
         out_fileHZ<<setprecision(numeric_limits<double>::digits10+1)<<WG_Hz[ii]<<endl;
         out_Is<<setprecision(numeric_limits<double>::digits10+1)<<Is[ii]<<endl;
         out_Vs<<setprecision(numeric_limits<double>::digits10+1)<<Vs[ii]<<endl;
         out_averageEY<<setprecision(numeric_limits<double>::digits10+1)<<average_Ey[ii]<<endl;

    }

    if( Ey_outputs.size() == WG_Ey.size())
    {
        out_additional_ey.open ( ( s+"_Ey_other_points.txt").c_str() );
        out_additional_ey<<this->WG_dl<<endl;
        out_additional_ey<<end_<<endl;

            for(int ii = 0 ; ii< end_; ii++)
            {
                out_additional_ey<<setprecision(numeric_limits<double>::digits10+1)<<Ey_outputs[ii]<<endl;
                cout<<Ey_outputs[ii]<<endl;
            }
    }
}

void WG_16:: write_output_file_for_yz_plane (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->WG_Ey_plane_yz[0].size();

    cout<< " Writing yz plane to file...." << endl;

    out_fileEY<<this->Nyy<<","<<this->Nzz<<","<<WG_Ey_plane_yz.size()<<endl;

        for(int ii = 0 ; ii< end_; ii++)
        {
            for (int dt=0; dt< WG_Ey_plane_yz.size() ; dt++)
            {
            out_fileEY<<setprecision(numeric_limits<double>::digits10+1) <<WG_Ey_plane_yz[dt][ii]<<",";
            }
           out_fileEY<<endl;
        }
}

void WG_16:: write_output_file_for_line (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_file_line;

    string s = nodefileEY.substr(0,nodefileEY.length() - 4);

    out_file_line.open((s+"_.csv").c_str());

    int end_ = this->WG_Ey.size(); // number of time steps

    cout<< " Writing Ey fields along line to file...." << endl;

    out_file_line<<this->Nz<< "," << end_<<endl;

    for (int dt=0; dt< end_ ; dt++){
        for(int ii = 0 ; ii< Nz; ii++)
        {
            out_file_line<<setprecision(numeric_limits<double>::digits10+1) <<Ey_field_along_line[dt][ii]<<",";
            //cout<< ii<<endl;
        }
        out_file_line<<endl;
    }
}


void WG_16:: write_output_file_for_xy_plane (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->WG_Ey_plane_xy[0].size();

    cout<< " Writing xy plane to file...." << endl;

    out_fileEY<<this->Nxx<<","<<this->Nyy<<","<<WG_Ey_plane_xy.size()<<endl;

     for(int ii = 0 ; ii< end_; ii++){
       for (int dt=0; dt< WG_Ey_plane_xy.size() ; dt++)
        {
            out_fileEY<<setprecision(numeric_limits<double>::digits10+1) <<WG_Ey_plane_xy[dt][ii]<<",";
        }
        out_fileEY<<endl;
    }
}

void WG_16:: write_output_file_for_zx_plane (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->WG_Ey_plane_zx[0].size();

    cout<< " Writing zx plane to file...." << endl;

    out_fileEY<<this->Nxx<<","<<this->Nzz<<","<<WG_Ey_plane_zx.size()<<endl;

       for(int ii = 0 ; ii< end_; ii++)
        {
            for (int dt=0; dt< WG_Ey_plane_zx.size() ; dt++)
            {
            out_fileEY<<setprecision(numeric_limits<double>::digits10+1) <<WG_Ey_plane_zx[dt][ii]<<",";
            }
        out_fileEY<<endl;
        }
}

void WG_16 :: write_incident_file_for_xy_plane (const string& nodefileEY)
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->V_incident_xy[0].size();

    cout<< " Writing INCIDENT in xy plane to file...." << endl;

    out_fileEY<<this->Nxx<<","<<this->Nyy<<","<<endl;

    for(int ii = 0 ; ii< end_; ii++)
    {
        for (int dt=0; dt< V_incident_xy.size() ; dt++)
        {
            out_fileEY<<setprecision(numeric_limits<double>::digits10+1) <<V_incident_xy[dt][ii]<<",";
        }
        out_fileEY<<endl;
    }


}

void WG_16::print_Ey_output()
{
    cout<< " PRINTING Ey......." <<endl;
    for(int i= 0; i<WG_Ex.size(); i++)
        cout<<i << ": "<< setprecision(15)<<WG_Ey[i] <<endl;
}

void WG_16::print_Ez_output()
{

    cout<< " PRINTING Ez......." <<endl;
    for(int i= 0; i<WG_Ez.size(); i++)
        cout<<i << ": "<< setprecision(15)<<WG_Ez[i] <<endl;
}

bool WG_16::print_G_PML_parameters(int L_R,int PML_layer)
{
    // initiate var with node_id corresponding to 1st node in pml layer i.e id=0;
    // multiply total node per layer by PML layer id. to
    int N_pml = this->npml ;
    int node_id = 0;
    int first_node = 0;
    int xy = this->Nx*this->Ny;

    cout<<endl;

    if(N_pml == 0) return false;

    if( PML_layer > N_pml) PML_layer = N_pml;

    if (L_R == 1) // left side of waveguide
    {
        node_id = (N_pml-PML_layer)*xy;

        for(int y=0; y< this->Ny ; y++){
            for (int x = 0; x< this->Nx ; x++) { cout<<this->WG_nodes[node_id].get_G(2)<<" ";  }
        cout<<endl<<endl;
        }
    return true;
    }

    if (L_R == 2) // right side of waveguide
    {

        first_node = 1 + return_coordinates_WG (this->WG_width,this->WG_height,this->WG_length,npml,this->WG_dl,this->WG_width,this->WG_height,this->WG_length);
        node_id = first_node + (PML_layer-1)*xy;

        for(int y=0; y< this->Ny ; y++){
            for (int x = 0; x< this->Nx ; x++) { cout<<this->WG_nodes[node_id].get_G(2)<<" ";  }
        cout<<endl;
        }
    return true;
    }

return false;
}

vector<double> WG_16::compute_far_field(const vector<int> &plane_xy1)
{/*
    int length_xy1 = plane_xy1.size();
    double M_xy1_i= 0, M_xy1_j= 0 , M_xy1_k= 0;
    double M_xy2_i= 0, M_xy2_j= 0,  M_xy2_k= 0;
    double M_yz1_i= 0, M_yz1_j= 0 , M_yz1_k= 0;
    double M_yz2_i= 0, M_yz2_j= 0,  M_yz2_k= 0;
    double M_zx1_i= 0, M_zx1_j= 0 , M_zx1_k= 0;
    double M_zx2_i= 0, M_zx2_j= 0,  M_zx2_k= 0;

    double J_xy1_i= 0, J_xy1_j= 0 , J_xy1_k= 0;
    double J_xy2_i= 0, J_xy2_j= 0,  J_xy2_k= 0;
    double J_yz1_i= 0, J_yz1_j= 0 , J_yz1_k= 0;
    double J_yz2_i= 0, J_yz2_j= 0,  J_yz2_k= 0;
    double J_zx1_i= 0, J_zx1_j= 0 , J_zx1_k= 0;
    double J_zx2_i= 0, J_zx2_j= 0,  J_zx2_k= 0;

    double x_= 0, y_= 0 , z_=0;
    double A_xy1_i= 0, A_xy1_j= 0, A_xy1_k= 0;
    double A_xy2_i= 0, A_xy2_j= 0, A_xy2_k= 0;
    double A_yz1_i= 0, A_yz1_j= 0, A_yz1_k= 0;
    double A_yz2_i= 0, A_yz2_j= 0, A_yz2_k= 0;
    double A_zx1_i= 0, A_zx1_j= 0, A_zx1_k= 0;
    double A_zx2_i= 0, A_zx2_j= 0, A_zx2_k= 0;

    for(int i = 0; i<length_xy1; i++)
    {
        //computing tangential E and H fields
       M_xy1_i= -this->Ey_output_at_node( plane_xy1[i] );
       M_xy1_j= this->Ex_output_at_node( plane_xy1[i] );
       J_xy1_i= this->Hy_output_at_node( plane_xy1[i] );
       J_xy1_j= -this->Hx_output_at_node( plane_xy1[i] );


       M_xy2_i= this->Ey_output_at_node( plane_xy2[i] );
       M_xy2_j= -this->Ex_output_at_node( plane_xy2[i] );
       J_xy2_i= -this->Hy_output_at_node( plane_xy2[i] );
       J_xy2_j= this->Hx_output_at_node( plane_xy2[i] );

       M_yz1_j= -this->Ez_output_at_node( plane_yz1[i] );
       M_yz1_k= this->Ey_output_at_node( plane_yz1[i] );
       J_yz1_j= this->Hz_output_at_node( plane_yz1[i] );
       J_yz1_k= -this->Hy_output_at_node( plane_yz1[i] );

       M_yz2_j= this->Ez_output_at_node( plane_yz2[i] );
       M_yz2_k= -this->Ey_output_at_node( plane_yz2[i] );
       J_yz2_j= -this->Hz_output_at_node( plane_yz2[i] );
       J_yz2_k= this->Hy_output_at_node( plane_yz2[i] );

       M_zx1_i= this->Ez_output_at_node( plane_zx1[i] );
       M_zx1_k= -this->Ex_output_at_node( plane_zx1[i] );
       J_zx1_i= -this->Hz_output_at_node( plane_zx1[i] );
       J_zx1_k= this->Hx_output_at_node( plane_zx1[i] );

       M_zx2_i= -this->Ez_output_at_node( plane_zx2[i] );
       M_zx2_k= this->Ex_output_at_node( plane_zx2[i] );
       J_zx2_i= this->Hz_output_at_node( plane_zx2[i] );
       J_zx2_k= -this->Hx_output_at_node( plane_zx2[i] );

       // Computing the discrete surface integral
       //xy1
       x_ = this->WG_nodes[ plane_xy1[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->WG_nodes[ plane_xy1[i] ].get_coord(2) - this->npml;
       z_ = this->WG_nodes[ plane_xy1[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_xy1_i = A_xy1_i + J_xy1_i *G_;
       A_xy1_j = A_xy1_j + J_xy1_j *G_;
       A_xy1_k = A_xy1_k + J_xy1_k *G_;

       F_xy1_i = F_xy1_i + M_xy1_i *G_;
       F_xy1_j = F_xy1_j + M_xy1_j *G_;
       F_xy1_k = F_xy1_k + M_xy1_k *G_;


       //xy2
       x_ = this->WG_nodes[ plane_xy2[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->WG_nodes[ plane_xy2[i] ].get_coord(2) - this->npml;
       z_ = this->WG_nodes[ plane_xy2[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_xy2_i = A_xy2_i + J_xy2_i *G_;
       A_xy2_j = A_xy2_j + J_xy2_j *G_;
       A_xy2_k = A_xy2_k + J_xy2_k *G_;

       F_xy2_i = F_xy2_i + M_xy2_i *G_;
       F_xy2_j = F_xy2_j + M_xy2_j *G_;
       F_xy2_k = F_xy2_k + M_xy2_k *G_;

       // Computing the discrete surface integral
       //yz1
       x_ = this->WG_nodes[ plane_yz1[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->WG_nodes[ plane_yz1[i] ].get_coord(2) - this->npml;
       z_ = this->WG_nodes[ plane_yz1[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_yz1_i = A_yz1_i + J_yz1_i *G_;
       A_yz1_j = A_yz1_j + J_yz1_j *G_;
       A_yz1_k = A_yz1_k + J_yz1_k *G_;

       F_yz1_i = F_yz1_i + M_yz1_i *G_;
       F_yz1_j = F_yz1_j + M_yz1_j *G_;
       F_yz1_k = F_yz1_k + M_yz1_k *G_;

       //yz2
       x_ = this->WG_nodes[ plane_yz2[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->WG_nodes[ plane_yz2[i] ].get_coord(2) - this->npml;
       z_ = this->WG_nodes[ plane_yz2[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_yz2_i = A_yz2_i + J_yz2_i *G_;
       A_yz2_j = A_yz2_j + J_yz2_j *G_;
       A_yz2_k = A_yz2_k + J_yz2_k *G_;

       F_yz2_i = F_yz2_i + F_yz2_i *G_;
       F_yz2_j = F_yz2_j + F_yz2_j *G_;
       F_yz2_k = F_yz2_k + F_yz2_k *G_;

        // Computing the discrete surface integral
       //zx1
       x_ = this->WG_nodes[ plane_zx1[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->WG_nodes[ plane_zx1[i] ].get_coord(2) - this->npml;
       z_ = this->WG_nodes[ plane_zx1[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_zx1_i = A_zx1_i + J_zx1_i *G_;
       A_zx1_j = A_zx1_j + J_zx1_j *G_;
       A_zx1_k = A_zx1_k + J_zx1_k *G_;

       F_zx1_i = F_zx1_i + M_zx1_i *G_;
       F_zx1_j = F_zx1_j + M_zx1_j *G_;
       F_zx1_k = F_zx1_k + M_zx1_k *G_;


       //zx2
       x_ = this->WG_nodes[ plane_zx2[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->WG_nodes[ plane_zx2[i] ].get_coord(2) - this->npml;
       z_ = this->WG_nodes[ plane_zx2[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_zx2_i = A_zx2_i + J_zx2_i *G_;
       A_zx2_j = A_zx2_j + J_zx2_j *G_;
       A_zx2_k = A_zx2_k + J_zx2_k *G_;

       F_zx2_i = F_zx2_i + M_zx2_i *G_;
       F_zx2_j = F_zx2_j + M_zx2_j *G_;
       F_zx2_k = F_zx2_k + M_zx2_k *G_;

    }

    A_i_ = A_xy1_i + A_xy2_i + A_yz1_i + A_yz2_i + A_zx1_i + A_zx2_i;
    A_j_ = A_xy1_j + A_xy2_j + A_yz1_j + A_yz2_j + A_zx1_j + A_zx2_j;
    A_k_ = A_xy1_k + A_xy2_k + A_yz1_k + A_yz2_k + A_zx1_k + A_zx2_k;

    F_i_ = F_xy1_i + F_xy2_i + F_yz1_i + F_yz2_i + F_zx1_i + F_zx2_i;
    F_j_ = F_xy1_j + F_xy2_j + F_yz1_j + F_yz2_j + F_zx1_j + F_zx2_j;
    F_k_ = F_xy1_k + F_xy2_k + F_yz1_k + F_yz2_k + F_zx1_k + F_zx2_k;

    this->A_i.push_back(A_i_);
    this->A_j.push_back(A_j_);
    this->A_k.push_back(A_k_);
    this->F_i.push_back(F_i_);
    this->F_j.push_back(F_j_);
    this->F_k.push_back(F_k_);
    */
}

void WG_16::analytical_fc()
{
  cout<< " The analytical resonant frequencies:  " <<endl;  // a = width ; b = height ; c = length ;

    int W = WG_width/WG_dl;
    W = W*WG_dl;
    int L = WG_length/WG_dl;
    L= L*WG_dl;

    cout<< W<<endl;
    cout << L<<endl;

            cout<<endl<<"***********************************RECTANGULAR WAVEGUIDE MODES RESONATOR********************************************"<<endl;
            cout<<"TE10 : ";
            cout<< sqrt((pi/WG_width)*(pi/WG_width) )*300000000000/(2*pi)<<endl;

            cout<<"**********************************************CAVITY RESONATOR***********************************************************"<<endl;

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

int return_coordinates_WG(double width, double height, double length,int npml, double dl, double x, double y, double z) //only waveguide
{
    int h = int( (height/dl) +0.5);
    int l = int ( (length/dl) + 0.5);
    int w = int ( (width/ dl) + 0.5);
    int s = h*w*npml;

    int c =  int(x/dl) +0.5;
    int d = int(y/dl) +0.5;
    int e = int(z/dl)+0.5;

    if(l == e) e = e-1;
    if (h == d) d = d-1;
    if (w == c) c= c - 1;

    int node_id =  c + d*w + e*(h*w);

    //cout<<x<<" "<<y <<" "<<z <<" "<<endl;
   // cout<<c<<" "<<d <<" "<<e <<" "<<" :"<<h << " "<<l <<" "<<w <<endl;
   // cin>>c;

    return node_id+s;
}

int return_coordinates_3D(double width, double height, double length,int n_PMl,
                          double WG_dl, double x, double y, double z)
{
    if( width  < x) x=width;
    if( height < y) y=height;
    if( length < z) z=length;

    int Nx = int ( (width  / WG_dl) + 0.5) +1;                     // note the added node which indicates position where the SD is placed.
    int Ny = int ( (height / WG_dl) + 0.5) +1;
    int Nz = int ( (length / WG_dl) + 0.5) +1;

    int x_id = int ( (x / WG_dl) + 0.5) +1;                        // note the added node which indicates position where the SD is placed.
    int y_id = int ( (y / WG_dl) + 0.5) +1;
    int z_id = int ( (z / WG_dl) + 0.5) +1;                        // 0->|_ _ _ _->width|

    //Number of nodes in entire computation along each coordinate axis
    int Nxx = Nx + 2*n_PMl;
    int Nyy = Ny + 2*n_PMl;
    int Nzz = Nz + 2*n_PMl;

    int node_id = Nxx*Nyy*(n_PMl+z_id-1) + Nxx*(n_PMl+y_id-1) + n_PMl+x_id-1;


    return node_id;
}

/*
int return_coordinates_2D(double width, double height,int npml, double dl, double x, double y, double z)
{

    int h = int( (height/dl) +0.5);
    int w = int ( (width/ dl) + 0.5);

    int h_ = h+2*npmly;
    int w_ = w+2*npmlx;


    int y_ = int(y/dl + 0.5);
    int x_ = int(x/dl + 0.5);


    int    coords = x_ + w_*((npmly+y_)-1) + npmlx-1;

    return coords;
}
*/
/*
vector<vector<double> > parse2DCsvFile(int &simtype, string inputFileName, vector< string > &e_filenames, vector< string > &h_filenames)
{

    vector< vector <double> > data;

    ifstream inputFile(inputFileName.c_str());              //input file streaming object
    int l = 0;
    int countt = 0;
    double temp =0;
    string str , str1("#E"),str2("#H"), str3("#s"), str4("#w"),str5("#q"),str6("#t");

    while (inputFile) {                              //returns -1 for EOF or EOL "\n"
        l++;
        string s, f;
        if (!getline(inputFile, s)) break;          //reads line from file stream into string s - checks if valid, breaks otherwise

        str = s.substr(0,2);
        //strll = "#H";

        //cout<<str<<endl;
         if(str == str3) simtype = 1;          // short_dipole
         if(str == str4) simtype = 0;          // waveguide
         if(str == str5) simtype = 2;          // naca0015
         if(str == str6) simtype = 3;          //  shunt dipole

         if ( str == str1 ){
            istringstream ss(s);

            while (ss){
                string line;
                if (!getline(ss, line, ',') ) break;
                if (countt >=2){                      //interested in the nth column.
                    e_filenames.push_back(line);
                    //cout<<line<<endl;
                }
                countt++;                           // have to read every element of string
            }
            countt = 0;
            } //"#H"

           if  (str == str2){
            istringstream ss(s);

            while (ss){
                string line;
                if (!getline(ss, line, ',') ) break;
                if (countt >=2){                      //interested in the nth column.
                    h_filenames.push_back(line);
                    //cout<<line<<endl;
                }
                countt++;                           // have to read every element of string
            }
            countt = 0;
            }

        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> v_temp;

            while (ss) {                           // splits the string up and converts each component to a double.
                string line;
                if (!getline(ss, line, ',') ) break;
                if (countt >=2){                      //interested in the nth column.
                    //try {                           // throws an exception if invalid argument i.e temp == string
                        stringstream s_;
                        s_ << line;
                        s_ >> temp;
                        v_temp.push_back(temp);
                        //temp = (stof(line));
                   // }
                   // catch (const std::invalid_argument e) {
                        //cout << "NaN found in file " << inputFileName << " line " << l
                        //<< endl;
                        //e.what();
                    //}
                }
                countt++;                           // have to read every element of string
            }
            countt = 0;
            data.push_back(v_temp);
        }
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return data;
}

*/
