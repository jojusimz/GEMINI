
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
#include <omp.h>

#include "HSCN_node.h"
#include "mesh_handler_HSCN.h"
#include "mesh_handler_2D.h"

using namespace std;

const   double pi=acos(-1.0L);
const   double  u0 = 12566370614e-16;
const   double e0 = 88541878176e-22;
const   double c = 299792458;
const   double Z0 = sqrt(u0/e0);//
//const   double er(1),ur(1);
//.....................................................................................................
//constructor  for waveguide
mesh_handler_HSCN::mesh_handler_HSCN(int simtype,   double width,   double height,   double length,   double dl,float ery, float urx,float sigma_ey, float sigma_mx, int n_PMl, int conduct_prof,  double Refn_factor)
{
    cout<<" Creating     Waveguide ......"<<endl;
    if ( dl <= 0 || width <=0 || height <=0 || length <=0)
    {
        mesh_handler_HSCN_width =  0;
        mesh_handler_HSCN_height = 0;
        mesh_handler_HSCN_length = 0;
        mesh_handler_HSCN_dl = 0;
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
        sim_type = simtype;
        centre_node = -1;
        mesh_handler_HSCN_width = width;
        mesh_handler_HSCN_height = height;
        mesh_handler_HSCN_length = length;
        mesh_handler_HSCN_dl = dl;
        npml = n_PMl;

        Nx = int ( (mesh_handler_HSCN_width / mesh_handler_HSCN_dl) + 0.5 );
        Nxx = Nx;
        Ny = int ( ( mesh_handler_HSCN_height / mesh_handler_HSCN_dl) + 0.5);
        Nyy = Ny;
        Nz = int ( (mesh_handler_HSCN_length / mesh_handler_HSCN_dl) + 0.5);
        Nzz = Nz;
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in waveguide

        cout<<"Ntotal:"<<Ntotal<<endl;
        cout<<"Nx * Ny "<< Nx*Ny*Nz<<endl;
        //cin>>Nx;

        mesh_handler_HSCN_PML_length = n_PMl * dl;
        Nzz =  2*n_PMl + Nz;
        Ntotal_ = Nzz * Nx *Ny;            // Ntotal + nodes in PML
        cout<<" Ntotal with PML:  "<<Ntotal_<<endl;

        int PML_boundary1 = n_PMl;
        int PML_boundary2 = n_PMl+Nz;
        mesh_handler_HSCN_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nt=Nx*Ny;

        bool special_nodex = false, special_nodey = false, special_nodez=false , spec = false;
        //Type of medium
        int medium_type = 1;
        media_identifier = medium_type ;

        // Defining the timestep
        float dl_ = 1e-3*dl;
        cout<< " maximum time : "<< 0.5*dl_/c<<endl;
        cout<<dl_ <<endl;

       // ery=5;
        //urx = 5;

        timestep = 0.5*dl_/c;//0.5*(sqrt( 2*ery/( 1/(dl_*dl_) + 1/ (dl_*dl_*urx) ))) / c;     //ur = urx and er = ery;
       // timestep = 0.2*timestep;
        cout<<"timestep "<<timestep<<endl;


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
                    // Creating HSCN_nodes;
                    if((z >= PML_boundary1 )&& (z < PML_boundary2 ) )
                    {
                        mesh_handler_HSCN_nodes.push_back(HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,false));
                    }
                    else
                    {
                        if ( z < PML_boundary1) mesh_handler_HSCN_nodes.push_back(HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,true,n_PMl,Refn_factor,conduct_prof,-1,-1,PML_boundary1-z-1));

                        else mesh_handler_HSCN_nodes.push_back(HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,true,n_PMl,Refn_factor,conduct_prof,-1,-1,z-PML_boundary2)); //cout<<"total_id : "<<id<<endl;cin>>id; }

                    }
                    id = id + 1;
                }
            }

            //defining the impedance and conductance for medium 1
            float urz =1, ury =1;
            float erx =1, erz =1;

            line_Y.assign(6,0);

            line_Y[0]= 2*c*dl_*timestep/(urz*dl_*dl_);     // yxy

            line_Y[1]= 2*c*dl_*timestep/(urx*dl_*dl_);     // yzy

            line_Y[2]= 2*c*dl_*timestep/(ury*dl_*dl_);    //yxz

            line_Y[3]= 2*erx*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urz*dl_) + dl_/(ury*dl_)) / dl_;  //yox

            line_Y[4]= 2*ery*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(urz*dl_)) / dl_;   //yoy

            line_Y[5]= 2*erz*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(ury*dl_)) / dl_;   //yoz

              double y = (2*ery*dl_/(c*timestep)  - (4*c*timestep * ( 1/urx + 1/urz ) / dl_));

            cout<<" IMPEDANCE "<<endl;
            cout<< " Y :" << line_Y[0] << " "<<line_Y[1] << " "<<line_Y[2] << " "<<line_Y[3]
            << " "<<line_Y[4] << " "<<line_Y[5] << " ..";
            cout<<setprecision(20)<<y <<endl<<endl;
            //cin>>y;

            float sigma_ex =0, sigma_ez =0;
            float sigma_my =0, sigma_mz =0;

            line_G.assign(3,0);

            line_G[0]= sigma_ex*dl_;

            line_G[1]= sigma_ey*dl_;

            line_G[2]= sigma_ez*dl_;


            line_R.assign(3,0);

            line_R[0]= sigma_mx*dl_;

            line_R[1]= sigma_my*dl_;

            line_R[2]= sigma_mz*dl_;

          // cin>>id;// cout<<"id"<<id<<endl;cin>>id;

          //defining the impedance and conductance for medium 2
            float urz2 =1, ury2 =1 , urx2 = 1;
            float erx2 =1, erz2 =1;
            float ery2 = 4.56;

            line_Y2.assign(6,0);

            line_Y2[0]= 2*c*dl_*timestep/(urz2*dl_*dl_);     // yxy

            line_Y2[1]= 2*c*dl_*timestep/(urx2*dl_*dl_);     // yzy

            line_Y2[2]= 2*c*dl_*timestep/(ury2*dl_*dl_);    //yxz

            line_Y2[3]= 2*erx2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urz2*dl_) + dl_/(ury2*dl_)) / dl_;  //yox

            line_Y2[4]= 2*ery2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx2*dl_) + dl_/(urz2*dl_)) / dl_;   //yoy

            line_Y2[5]= 2*erz2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx2*dl_) + dl_/(ury2*dl_)) / dl_;   //yoz

              double y2 = (2*ery2*dl_/(c*timestep)  - (4*c*timestep * ( 1/urx2 + 1/urz2 ) / dl_));

            cout<<" IMPEDANCE "<<endl;
            cout<< " Y :" << line_Y2[0] << " "<<line_Y2[1] << " "<<line_Y2[2] << " "<<line_Y2[3]
            << " "<<line_Y2[4] << " "<<line_Y2[5] << " ..";
            cout<<setprecision(20)<<y2 <<endl<<endl;
            //cin>>y;

            float sigma_ex2 =0, sigma_ez2 =0;
            float sigma_my2 =0, sigma_mz2 =0 , sigma_mx2 = 0;
            float sigma_ey2 = 8e3;

            line_G2.assign(3,0);

            line_G2[0]= sigma_ex2*dl_;

            line_G2[1]= sigma_ey2*dl_;

            line_G2[2]= sigma_ez2*dl_;


            line_R2.assign(3,0);

            line_R2[0]= sigma_mx2*dl_;

            line_R2[1]= sigma_my2*dl_;

            line_R2[2]= sigma_mz2*dl_;

          // cin>>id;// cout<<"id"<<id<<endl;cin>>id;
    }

}

//constructor 3D domain
mesh_handler_HSCN::mesh_handler_HSCN(  double width,   double height,   double length, int simtype,   double dl,
                         float ery, float urx,float sigma_ey, float sigma_mx, int n_PMl, int conduct_prof,  double Refn_factor)
{
    cout<<" Creating Cubic domain ......"<<endl;
    cout<<" dl "<< dl<<" length "<<length<<endl;

    if ( dl <= 0 || int(length) <=0 )
    {
        mesh_handler_HSCN_width =  0;
        mesh_handler_HSCN_height = 0;
        mesh_handler_HSCN_length = 0;
        mesh_handler_HSCN_dl = 0;
        Nx = 0;
        Ny = 0;
        Nz = 0;
        Ntotal = 0;
        check_meshed=0;
        cout<< " Troublesome space discretization value chosen " <<endl;
        cin>>mesh_handler_HSCN_width;
    }
    else
    {
        neighbours_set = false;
        sim_type = simtype;
        mesh_handler_HSCN_width = width;
        mesh_handler_HSCN_height = height;
        mesh_handler_HSCN_length = length;
        mesh_handler_HSCN_dl = dl;
        npml = n_PMl;

        //Number of nodes in the problem domain a  each coordinate axis
        Nx = int ( (mesh_handler_HSCN_width / mesh_handler_HSCN_dl)  + 0.5) +1;                     // note the added node which indicates position where the SD is placed.
        Ny = int ( (mesh_handler_HSCN_height / mesh_handler_HSCN_dl) + 0.5) +1;
        Nz = int ( (mesh_handler_HSCN_length / mesh_handler_HSCN_dl) + 0.5) +1;
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in waveguide

        int mid_x = int(0.5*width/dl  + 0.5)+ n_PMl;
        int mid_y = int(0.5*height/dl + 0.5)+ n_PMl;
        int mid_z = int(0.5*length/dl + 0.5)+ n_PMl;

        //Number of nodes in entire computation a  each coordinate axis
        Nxx = Nx + 2*n_PMl;
        Nyy = Ny + 2*n_PMl;
        Nzz = Nz + 2*n_PMl;
        Ntotal_ = Nxx*Nyy*Nzz;                  // total number of nodes in waveguide

        cout<<Nxx<<": "<<Nyy<<": "<<Nzz<<endl;
        mesh_handler_HSCN_PML_length = n_PMl * dl;

        int PML_boundaryx1 = n_PMl;
        int PML_boundaryx2 = Nx + n_PMl;
        int PML_boundaryy1 = n_PMl;
        int PML_boundaryy2 = Ny + n_PMl;
        int PML_boundaryz1 = n_PMl;
        int PML_boundaryz2 = Nz + n_PMl;

        //mesh_handler_HSCN_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nxy = Nxx*Nyy;

        cout<<"total node in problem domain:"<<Ntotal<<endl;
        cout<<"Nx * Ny *Nz "<< Ntotal_<<endl;
        //cout<< " maximum vector " << vector<HSCN_node>::max_size()<<endl;
        mesh_handler_HSCN_nodes.resize(Ntotal_);

        bool special_nodex = false, special_nodey = false, special_nodez=false , spec = false;
        int medium_type = 1;
        media_identifier = medium_type ;

          // Defining the timestep
          double dl_ = 1e-3*dl;
        timestep = 0.5*dl_/c;//0.5*(sqrt( 2*ery/( 1/(dl_*dl_) + 1/ (dl_*dl_*urx) ))) / c;     //ur = urx and er = ery;

        for(int z=0; z<Nzz; z++){
            for (int y=0 ; y<Nyy ; y++)
            {
                for (int x=0 ; x<Nxx ; x++)
                {
                    //mesh_handler_HSCN_nodes.resize(id+1);
                    //Setting flag for nodes at boundaries/ edge of domain
                    special_nodex = (( x==0)||( x==Nxx-1));
                    special_nodey = (( y==0)||( y==Nyy-1));
                    special_nodez = (( z==0)||( z==Nzz-1));
                    spec = special_nodex || special_nodey || special_nodez ;
                    // Creating HSCN_nodes;

                    if(((x >= PML_boundaryx1 )&& (x < PML_boundaryx2 ))&& ((y >= PML_boundaryy1 )&& (y < PML_boundaryy2 )) &&((z >= PML_boundaryz1 )&& (z < PML_boundaryz2 )) )
                    {
                       if( (x == mid_x )&&(z == mid_z)&& (y== mid_y)) { centre_node = id ; cout<<" id : "<<centre_node<<endl; cout<<"z plane "<<mid_z<<endl; }

                       mesh_handler_HSCN_nodes[id] = HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,false);
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

                       // mesh_handler_HSCN_nodes.push_back(HSCN_node(dl,tfactor,id,x,y,z,spec,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz));
                        mesh_handler_HSCN_nodes[id] = HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz);

                    }
                    id = id + 1;
                }
            }

    }

            //defining the impedance and conductance
            float urz =1, ury =1;
            float erx =1, erz =1;

            line_Y.assign(6,0);

            line_Y[0]= 2*c*dl_*timestep/(urz*dl_*dl_);     // yxy

            line_Y[1]= 2*c*dl_*timestep/(urx*dl_*dl_);     // yzy

            line_Y[2]= 2*c*dl_*timestep/(ury*dl_*dl_);    //yxz

            line_Y[3]= 2*erx*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urz*dl_) + dl_/(ury*dl_)) / dl_;  //yox

            line_Y[4]= 2*ery*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(urz*dl_)) / dl_;   //yoy

            line_Y[5]= 2*erz*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(ury*dl_)) / dl_;   //yoz


              double y = (2*ery*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(urz*dl_)) / dl_);

            cout<<" IMPEDANCE "<<endl;
            cout<< " Y :" << line_Y[0] << " "<<line_Y[1] << " "<<line_Y[2] << " "<<line_Y[3]
            << " "<<line_Y[4] << " "<<line_Y[5] << " ..";
            cout<<setprecision(20)<<y <<endl<<endl;

            // cin>>

            float sigma_ex =0, sigma_ez =0;
            float sigma_my =0, sigma_mz =0;

            line_G.assign(3,0);

            line_G[0]= sigma_ex*dl_;

            line_G[1]= sigma_ey*dl_;

            line_G[2]= sigma_ez*dl_;


            line_R.assign(3,0);

            line_R[0]= sigma_mx*dl_;

            line_R[1]= sigma_my*dl_;

            line_R[2]= sigma_mz*dl_;

              //defining the impedance and conductance for medium 2
            float urz2 =1, ury2 =1, urx2 = 1;
            float erx2 =1, erz2 =1;
            float ery2 = 2;

            line_Y2.assign(6,0);

            line_Y2[0]= 2*c*dl_*timestep/(urz2*dl_*dl_);     // yxy

            line_Y2[1]= 2*c*dl_*timestep/(urx2*dl_*dl_);     // yzy

            line_Y2[2]= 2*c*dl_*timestep/(ury2*dl_*dl_);    //yxz

            line_Y2[3]= 2*erx2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urz2*dl_) + dl_/(ury2*dl_)) / dl_;  //yox

            line_Y2[4]= 2*ery2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx2*dl_) + dl_/(urz2*dl_)) / dl_;   //yoy

            line_Y2[5]= 2*erz2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx2*dl_) + dl_/(ury2*dl_)) / dl_;   //yoz

              double y2 = (2*ery2*dl_/(c*timestep)  - (4*c*timestep * ( 1/urx2 + 1/urz2 ) / dl_));

            cout<<" IMPEDANCE "<<endl;
            cout<< " Y :" << line_Y2[0] << " "<<line_Y2[1] << " "<<line_Y2[2] << " "<<line_Y2[3]
            << " "<<line_Y2[4] << " "<<line_Y2[5] << " ..";
            cout<<setprecision(20)<<y2 <<endl<<endl;
            //cin>>y;

            float sigma_ex2 =0, sigma_ez2 =0;
            float sigma_my2 =0, sigma_mz2 =0 , sigma_mx2=0;
            float sigma_ey2 = 10e3;

            line_G2.assign(3,0);

            line_G2[0]= sigma_ex2*dl_;

            line_G2[1]= sigma_ey2*dl_;

            line_G2[2]= sigma_ez2*dl_;


            line_R2.assign(3,0);

            line_R2[0]= sigma_mx2*dl_;

            line_R2[1]= sigma_my2*dl_;

            line_R2[2]= sigma_mz2*dl_;
    }
}

//2D_HSCN
mesh_handler_HSCN::mesh_handler_HSCN(  double width,   double height, int simtype,   double dl,float ery, float urx,float sigma_ey, float sigma_mx, int n_PMl, int conduct_prof,  double Refn_factor)
{

    cout<<" Creating Planar domain for HSCN ......"<<endl;
    cout<<" dl "<< dl<<" length "<< dl<<endl;

    if ( dl <= 0  )
    {
        mesh_handler_HSCN_width =  0;
        mesh_handler_HSCN_height = 0;
        mesh_handler_HSCN_length = 0;
        mesh_handler_HSCN_dl = 0;
        Nx = 0;
        Ny = 0;
        Nz = 0;
        Ntotal = 0;
        check_meshed=0;
        cout<< " Troublesome space discretization value chosen " <<endl;
        cin>>mesh_handler_HSCN_width;
    }
    else
    {
        neighbours_set = false;
        sim_type = simtype;         //sim type - waveguide, dipole etc.
        mesh_handler_HSCN_width = width;
        mesh_handler_HSCN_height = height;
        mesh_handler_HSCN_length = dl;
        mesh_handler_HSCN_dl = dl;
        npml = n_PMl;
        centre_node = get_coordinate_iD_2D(width,height,npml,npml,dl,width/2,height/2,dl);

        //Number of nodes in the problem domain a  each coordinate axis
        Nx = int ( (mesh_handler_HSCN_width / mesh_handler_HSCN_dl)  + 0.5);
        Ny = int ( (mesh_handler_HSCN_height / mesh_handler_HSCN_dl) + 0.5);
        Nz = int ( (mesh_handler_HSCN_length / mesh_handler_HSCN_dl) + 0.5);
        Ntotal = Nx*Ny*Nz;                  // total number of nodes in computational domain

        //Number of nodes in entire computation a  each coordinate axis
        Nxx = Nx + 2*n_PMl;
        Nyy = Ny + 2*n_PMl;
        Nzz = Nz;
        Ntotal_ = Nxx*Nyy*Nzz;                  // total number of nodes



       //Reserve memory for mesh_handler_HSCN nodes
        //mesh_handler_HSCN_nodes.reserve(Ntotal_);

        cout<<Nxx<<" : "<<Nyy<<" : "<<Nz<<endl;
        mesh_handler_HSCN_PML_length = n_PMl * dl;

        int PML_boundaryx1 = n_PMl;
        int PML_boundaryx2 = Nx + n_PMl;
        int PML_boundaryy1 = n_PMl;
        int PML_boundaryy2 = Ny + n_PMl;

        //mesh_handler_HSCN_node_neighbours.resize(Ntotal_);
        int id= 0 ;
        int Nxy = Nxx*Nyy;

        cout<<"Total node in problem domain:"<<Ntotal<<endl;
        cout<<"Nx * Ny *Nz "<< Ntotal_<<endl;
        mesh_handler_HSCN_nodes.resize(Ntotal_);

        bool special_nodex = false, special_nodey = false, spec = false;
       int medium_type = 1;
        media_identifier = medium_type ;


         // Defining the timestep
          double dl_ = 1e-3*dl;
        timestep = 0.5*(sqrt( 2*ery/( 1/(dl_*dl_) + 1/ (dl_*dl_*urx) ))) / c;     //ur = urx and er = ery;

        for(int z=0; z<Nz; z++){
            for (int y=0 ; y<Nyy ; y++)
            {
                for (int x=0 ; x<Nxx ; x++)
                {
                    //mesh_handler_HSCN_nodes.resize(id+1);
                    //Setting flag for nodes at boundaries/ edge of domain
                    special_nodex = ( ( x==0)||( x==Nxx-1) );
                    special_nodey = ( ( y==0)||( y==Nyy-1) );
                    spec = special_nodex || special_nodey ;
                    // Creating HSCN_nodes;

                    if(((x >= PML_boundaryx1 )&& (x < PML_boundaryx2 ))&& ((y >= PML_boundaryy1 )&& (y < PML_boundaryy2 )))
                    {
                        mesh_handler_HSCN_nodes[id] = HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,false);
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

                        mesh_handler_HSCN_nodes[id] = (HSCN_node(0,dl,timestep,id,x,y,z,spec,medium_type,false,true,n_PMl,Refn_factor,conduct_prof,Lx,Ly,Lz));

                    }
                    id = id + 1;
                }
            }

    }
       //defining the impedance and conductance

            float urz =1, ury =1;
            float erx =1, erz =1;

            line_Y.assign(6,0);

            line_Y[0]= 2*c*dl_*timestep/(urz*dl_*dl_);     // yxy

            line_Y[1]= 2*c*dl_*timestep/(urx*dl_*dl_);     // yzy

            line_Y[2]= 2*c*dl_*timestep/(ury*dl_*dl_);    //yxz

            line_Y[3]= 2*erx*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urz*dl_) + dl_/(ury*dl_)) / dl_;  //yox

            line_Y[4]= 2*ery*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(urz*dl_)) / dl_;   //yoy

            line_Y[5]= 2*erz*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx*dl_) + dl_/(ury*dl_)) / dl_;   //yoz

            float sigma_ex =0, sigma_ez =0;
            float sigma_my =0, sigma_mz =0;

            line_G.assign(3,0);

            line_G[0]= sigma_ex*dl_;

            line_G[1]= sigma_ey*dl_;

            line_G[2]= sigma_ez*dl_;


            line_R.assign(3,0);

            line_R[0]= sigma_mx*dl_;

            line_R[1]= sigma_my*dl_;

            line_R[2]= sigma_mz*dl_;

              //defining the impedance and conductance for medium 2
            float urz2 =1, ury2 =1, urx2 = 1;
            float erx2 =1, erz2 =1;
            float ery2 = 4.56;                              // <<<<<<<<<<------------------------------------------------------------------------------------- change this

            line_Y2.assign(6,0);

            line_Y2[0]= 2*c*dl_*timestep/(urz2*dl_*dl_);     // yxy

            line_Y2[1]= 2*c*dl_*timestep/(urx2*dl_*dl_);     // yzy

            line_Y2[2]= 2*c*dl_*timestep/(ury2*dl_*dl_);    //yxz

            line_Y2[3]= 2*erx2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urz2*dl_) + dl_/(ury2*dl_)) / dl_;  //yox

            line_Y2[4]= 2*ery2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx2*dl_) + dl_/(urz2*dl_)) / dl_;   //yoy

            line_Y2[5]= 2*erz2*dl_*dl_/(c*timestep*dl_)  - 4*c*timestep * ( dl_/(urx2*dl_) + dl_/(ury2*dl_)) / dl_;   //yoz

              double y2 = (2*ery2*dl_/(c*timestep)  - (4*c*timestep * ( 1/urx2 + 1/urz2 ) / dl_));

            cout<<" IMPEDANCE "<<endl;
            cout<< " Y :" << line_Y2[0] << " "<<line_Y2[1] << " "<<line_Y2[2] << " "<<line_Y2[3]
            << " "<<line_Y2[4] << " "<<line_Y2[5] << " ..";
            cout<<setprecision(20)<<y2 <<endl<<endl;
            //cin>>y;

            float sigma_ex2 =0, sigma_ez2 =0;
            float sigma_my2 =0, sigma_mz2 =0, sigma_mx2=0;
            float sigma_ey2 = 8e3;                // <<<<<<<<<<------------------------------------------------------------------------------------- change this

            line_G2.assign(3,0);

            line_G2[0]= sigma_ex2*dl_;

            line_G2[1]= sigma_ey2*dl_;

            line_G2[2]= sigma_ez2*dl_;


            line_R2.assign(3,0);

            line_R2[0]= sigma_mx2*dl_;

            line_R2[1]= sigma_my2*dl_;

            line_R2[2]= sigma_mz2*dl_;

    }

}


mesh_handler_HSCN::mesh_handler_HSCN(const mesh_handler_HSCN& copy_mesh_handler_HSCN)
{
    mesh_handler_HSCN_width = copy_mesh_handler_HSCN.mesh_handler_HSCN_width;
    mesh_handler_HSCN_height = copy_mesh_handler_HSCN.mesh_handler_HSCN_height;
    mesh_handler_HSCN_length = copy_mesh_handler_HSCN.mesh_handler_HSCN_length;
    mesh_handler_HSCN_dl = copy_mesh_handler_HSCN.mesh_handler_HSCN_dl;            // discretization length
    mesh_handler_HSCN_PML_length = copy_mesh_handler_HSCN.mesh_handler_HSCN_PML_length;
    sim_type = copy_mesh_handler_HSCN.sim_type;
    Nx = copy_mesh_handler_HSCN.Nx;                 //number of nodes on the x , y and z axis
    Ny = copy_mesh_handler_HSCN.Ny;
    Nz = copy_mesh_handler_HSCN.Nz;
    Nzz = copy_mesh_handler_HSCN.Nzz;
    Nxx = copy_mesh_handler_HSCN.Nxx;
    Nyy = copy_mesh_handler_HSCN.Nyy;
    npml = copy_mesh_handler_HSCN.npml;
    centre_node = copy_mesh_handler_HSCN.centre_node;         //node of excitation for short dipole - default set to -1 if waveguide
    mesh_handler_HSCN_Ex = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ex;
    mesh_handler_HSCN_Ey = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey;
    mesh_handler_HSCN_Ez = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ez;
    mesh_handler_HSCN_Hx = copy_mesh_handler_HSCN.mesh_handler_HSCN_Hx;
    mesh_handler_HSCN_Hy = copy_mesh_handler_HSCN.mesh_handler_HSCN_Hy;
    mesh_handler_HSCN_Hz = copy_mesh_handler_HSCN.mesh_handler_HSCN_Hz;
    mesh_handler_HSCN_Ey_plane_xy = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey_plane_xy;
    mesh_handler_HSCN_Ey_plane_yz = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey_plane_yz;
    mesh_handler_HSCN_Ey_plane_zx = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey_plane_zx;
    Ntotal= copy_mesh_handler_HSCN.Ntotal;
    Ntotal_ = copy_mesh_handler_HSCN.Ntotal_;
    mesh_handler_HSCN_nodes = copy_mesh_handler_HSCN.mesh_handler_HSCN_nodes;   // represents the waveguide nodes; this is set in the constructor depending on the simulation parameters.
    mesh_handler_HSCN_node_neighbours= copy_mesh_handler_HSCN.mesh_handler_HSCN_node_neighbours;
    check_meshed = copy_mesh_handler_HSCN.check_meshed;         //checks if the waveguide has been meshed / filled with cubes
    neighbours_set = copy_mesh_handler_HSCN.neighbours_set;
    timestep = copy_mesh_handler_HSCN.timestep;
    line_Y = copy_mesh_handler_HSCN.line_Y;
    line_G = copy_mesh_handler_HSCN.line_G;
    line_R = copy_mesh_handler_HSCN.line_R;
    line_Y2 = copy_mesh_handler_HSCN.line_Y2;
    line_G2 = copy_mesh_handler_HSCN.line_G2;
    line_R2 = copy_mesh_handler_HSCN.line_R2;
    media_identifier = copy_mesh_handler_HSCN.media_identifier;
    vctr_medium_iD = copy_mesh_handler_HSCN.vctr_medium_iD;
    vctr_PML_medium_iD = copy_mesh_handler_HSCN.vctr_PML_medium_iD;
}

mesh_handler_HSCN &mesh_handler_HSCN::operator=(const mesh_handler_HSCN& copy_mesh_handler_HSCN)
{
    if (this!= &copy_mesh_handler_HSCN)
    {
        mesh_handler_HSCN_width = copy_mesh_handler_HSCN.mesh_handler_HSCN_width;
        mesh_handler_HSCN_height = copy_mesh_handler_HSCN.mesh_handler_HSCN_height;
        mesh_handler_HSCN_length = copy_mesh_handler_HSCN.mesh_handler_HSCN_length;
        mesh_handler_HSCN_dl = copy_mesh_handler_HSCN.mesh_handler_HSCN_dl;            // discretization length
        mesh_handler_HSCN_PML_length = copy_mesh_handler_HSCN.mesh_handler_HSCN_PML_length;
        sim_type = copy_mesh_handler_HSCN.sim_type;
        Nx = copy_mesh_handler_HSCN.Nx;                 //number of nodes on the x , y and z axis
        Ny = copy_mesh_handler_HSCN.Ny;
        Nz = copy_mesh_handler_HSCN.Nz;
        Nzz = copy_mesh_handler_HSCN.Nzz;
        Nxx = copy_mesh_handler_HSCN.Nxx;
        Nyy = copy_mesh_handler_HSCN.Nyy;
        npml = copy_mesh_handler_HSCN.npml;
        centre_node = copy_mesh_handler_HSCN.centre_node;         //node of excitation for short dipole - default set to -1 if waveguide
        mesh_handler_HSCN_Ex = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ex;
        mesh_handler_HSCN_Ey = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey;
        mesh_handler_HSCN_Ez = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ez;
        mesh_handler_HSCN_Hx = copy_mesh_handler_HSCN.mesh_handler_HSCN_Hx;
        mesh_handler_HSCN_Hy = copy_mesh_handler_HSCN.mesh_handler_HSCN_Hy;
        mesh_handler_HSCN_Hz = copy_mesh_handler_HSCN.mesh_handler_HSCN_Hz;
        mesh_handler_HSCN_Ey_plane_xy = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey_plane_xy;
        mesh_handler_HSCN_Ey_plane_yz = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey_plane_yz;
        mesh_handler_HSCN_Ey_plane_zx = copy_mesh_handler_HSCN.mesh_handler_HSCN_Ey_plane_zx;
        Ntotal= copy_mesh_handler_HSCN.Ntotal;
        Ntotal_ = copy_mesh_handler_HSCN.Ntotal_;
        mesh_handler_HSCN_nodes = copy_mesh_handler_HSCN.mesh_handler_HSCN_nodes;   // represents the waveguide nodes; this is set in the constructor depending on the simulation parameters.
        mesh_handler_HSCN_node_neighbours= copy_mesh_handler_HSCN.mesh_handler_HSCN_node_neighbours;
        check_meshed = copy_mesh_handler_HSCN.check_meshed;         //checks if the waveguide has been meshed / filled with cubes
        neighbours_set = copy_mesh_handler_HSCN.neighbours_set;
        timestep = copy_mesh_handler_HSCN.timestep;
        line_Y = copy_mesh_handler_HSCN.line_Y;
        line_G = copy_mesh_handler_HSCN.line_G;
        line_R = copy_mesh_handler_HSCN.line_R;
        line_Y2 = copy_mesh_handler_HSCN.line_Y2;
        line_G2 = copy_mesh_handler_HSCN.line_G2;
        line_R2 = copy_mesh_handler_HSCN.line_R2;
        media_identifier = copy_mesh_handler_HSCN.media_identifier;
        vctr_medium_iD = copy_mesh_handler_HSCN.vctr_medium_iD;
        vctr_PML_medium_iD = copy_mesh_handler_HSCN.vctr_PML_medium_iD;
    }
    return (*this);
}

 mesh_handler_HSCN::~mesh_handler_HSCN()
{
  cout<<"destroyed"<<endl;//destructor
}

void mesh_handler_HSCN::E_output_at_node(int output_node)
{
      double Vx (0), Vy(0) , Vz(0);

    //if (mesh_handler_HSCN_nodes[output_node].is_PML_HSCN()==1) { cout<<"CAUTION!!: output node is in PML"<<endl<<endl;}

    if (mesh_handler_HSCN_nodes[output_node].is_PML_HSCN()==1)
    {

        if (media_identifier == 1)
        {
            mesh_handler_HSCN_nodes[output_node].compute_complex_Vj( line_Y, line_G , Vx, Vy ,Vz );
        }

        else if (media_identifier == 2)
        {
            mesh_handler_HSCN_nodes[output_node].compute_complex_Vj( line_Y2, line_G2 , Vx, Vy ,Vz );
        }
    }
    else
    {
        if (media_identifier == 1)
        {
            mesh_handler_HSCN_nodes[output_node].compute_standard_Vj( line_Y, line_G , Vx, Vy ,Vz );
        }

        else if (media_identifier == 2)
        {
            mesh_handler_HSCN_nodes[output_node].compute_standard_Vj( line_Y2, line_G2 , Vx, Vy ,Vz );
        }
    }

    mesh_handler_HSCN_Ex.push_back(Ex_output_at_node(output_node , Vx));

    mesh_handler_HSCN_Ey.push_back(Ey_output_at_node(output_node , Vy));

    mesh_handler_HSCN_Ez.push_back(Ez_output_at_node(output_node , Vz));

}

void mesh_handler_HSCN::H_output_at_node(int output_node)
{

    mesh_handler_HSCN_Hx.push_back(Hx_output_at_node(output_node ));

    mesh_handler_HSCN_Hy.push_back(Hy_output_at_node(output_node ));

    mesh_handler_HSCN_Hz.push_back(Hz_output_at_node(output_node ));
}


  double mesh_handler_HSCN::Ex_output_at_node(int node_id,   double Vx)
{
    return (-Vx /(1e-3*this->mesh_handler_HSCN_dl) );
}

  double mesh_handler_HSCN::Ey_output_at_node(int node_id,   double Vy)
{
    return (-Vy /(1e-3*this->mesh_handler_HSCN_dl) );
}

  double mesh_handler_HSCN::Ez_output_at_node(int node_id,   double Vz)
{
    return (-Vz /(1e-3*this->mesh_handler_HSCN_dl) );
}

  double mesh_handler_HSCN::Hx_output_at_node(int node_id,   double ix)
{
    return 0;
}

  double mesh_handler_HSCN::Hy_output_at_node(int node_id ,   double iy)
{
    return 0;
}

  double mesh_handler_HSCN::Hz_output_at_node(int node_id ,   double iz)
{
    return 0;
}

void mesh_handler_HSCN::set_special_nodes()
{
    if(neighbours_set == false)
    {
        int temp;
        cout<<" neighbours have not been set!!" <<endl;
        cout<<" press any number to proceed...neighbours will now be set"<<endl;
        cin>>temp;
        this->set_mesh_handler_HSCN_neighbour();               // potential problems could arise here when parallelization of the code
    }

    for (int mesh_handler_HSCN_i = 0; mesh_handler_HSCN_i<this->mesh_handler_HSCN_nodes.size(); mesh_handler_HSCN_i++)
    {
    int numbr_neigh = mesh_handler_HSCN_node_neighbours[mesh_handler_HSCN_i].size();

        for (int i =0; i<numbr_neigh;i++)
        {
            if ( (  this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[mesh_handler_HSCN_i][i] ].is_PEC() ) &&   (   this->mesh_handler_HSCN_nodes[mesh_handler_HSCN_i].check_special_node())  )
            {
                int temp;
                //cout<<" |[][x] A PEC node is one node away from boundary!!"<<endl;
                //  cout<<"This will create problems in the connection process!!"<<endl;
                // cout<<"Suggestion is to extend the domain in the " << i <<"direction to create |[][][x]"<<endl;

            }
            if ( (  this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[mesh_handler_HSCN_i][i] ].is_PEC() ) &&   (  !this->mesh_handler_HSCN_nodes[mesh_handler_HSCN_i].check_special_node())  )
            {
                this->mesh_handler_HSCN_nodes[mesh_handler_HSCN_i].set_special_HSCN_node(true);
                break;
            }
        }
    }
}


void mesh_handler_HSCN::place_dipole_antenna(  double dipole_length, int x, int y, int z )
{
    cout<<"dipole length: "<<dipole_length<<endl<<endl;

    if( (x==-1)||(y==-1)||(z==-1) )
    {
        int sd_node = this->centre_node;
        int antenna_node_no = dipole_length/this->mesh_handler_HSCN_dl + 0.5;
        int max_l = ( Ny-1 );

        if (antenna_node_no>max_l)
        {
            cout<<"CAUTION: Antenna length bigger than geometry!!"<<endl;
            cout<<"Setting to full length mode"<<endl;
            antenna_node_no = max_l;
        }

        for(int i = 1; i<=int(antenna_node_no*0.5+0.5) ; i++)
        {
            this->mesh_handler_HSCN_nodes[sd_node +i*Nxx].set_PEC();
            this->mesh_handler_HSCN_nodes[sd_node -i*Nxx].set_PEC();
            cout<<i;
        }
    }
}


void mesh_handler_HSCN::set_mesh_handler_HSCN_neighbour()
{
    int id = 0;
    int Nt = this->Nxx*this->Nyy;

    mesh_handler_HSCN_node_neighbours.resize(this->Ntotal_);

    for(int z=0; z<this->Nzz; z++)
    {
        //cout<<id<<endl;
        for (int y=0 ; y<this->Nyy ; y++)
        {
            for (int x=0 ; x<this->Nxx ; x++)
            {
                if (x!=0)
                {   //cout<<mesh_handler_HSCN_node_neighbours.size()<<endl;
                //cin>>id;
                    mesh_handler_HSCN_node_neighbours[id].push_back(mesh_handler_HSCN_nodes[id-1].get_iD());
                    mesh_handler_HSCN_node_neighbours[id-1].push_back(mesh_handler_HSCN_nodes[id].get_iD());
                }
                if (y!=0)
                {
                    mesh_handler_HSCN_node_neighbours[id].push_back(mesh_handler_HSCN_nodes[id-this->Nxx].get_iD());
                    mesh_handler_HSCN_node_neighbours[id-this->Nxx].push_back(mesh_handler_HSCN_nodes[id].get_iD());
                }
                if(z!=0)
                {
                    mesh_handler_HSCN_node_neighbours[id].push_back(mesh_handler_HSCN_nodes[id-Nt].get_iD());
                    mesh_handler_HSCN_node_neighbours[id-Nt].push_back(mesh_handler_HSCN_nodes[id].get_iD());
                }
                id = id+1;
            }
        }
    }
    neighbours_set = true;
    cout<<"neighbours set"<<endl;

}

int::mesh_handler_HSCN::get_node_coord(int node_iD, int coord_x_y_z)
{
    return( this->mesh_handler_HSCN_nodes[node_iD].get_coord(coord_x_y_z) );
}
void mesh_handler_HSCN::print_HSCN_neighbour(int node_id)
{
    int countt = mesh_handler_HSCN_node_neighbours[node_id].size();
    for ( int i = 0; i<countt; i++)
        cout<<mesh_handler_HSCN_node_neighbours[node_id][i];
}

void mesh_handler_HSCN::print_mesh_handler_HSCN_nodes()
{
    for( int i=0; i< this->mesh_handler_HSCN_nodes.size(); i++)
    {
        cout<<this->mesh_handler_HSCN_nodes[i]<<endl;
    }
}


void mesh_handler_HSCN::print_mesh_handler_HSCN_node(int _node)
{
    cout<<" Node : " <<_node << endl;
    for(int i=0; i<12; i++)
    {
        cout<<this->mesh_handler_HSCN_nodes[_node].get_incidence(i+1)<<endl;
    }

}

void mesh_handler_HSCN::print_mesh_handler_HSCN_plane_par(int z_plane, int par )
{
    cout<<" Plane : " << z_plane << endl;

 int node_id = (z_plane-1) * this->Nxx*Nyy;

 for( int y = 0; y< Nyy ; y++){
    for ( int x = 0; x< Nxx; x++){
        if( par == 1) cout<<this->mesh_handler_HSCN_nodes[node_id].is_PML_HSCN()<<" :  ";
        else if( par == 2) cout<<this->mesh_handler_HSCN_nodes[node_id].get_iD()<<" :  ";
        else if( par == 3) cout<<this->mesh_handler_HSCN_nodes[node_id].is_PEC()<<" :  ";
        else if ( par ==4 ) cout<< this->mesh_handler_HSCN_nodes[node_id].get_medium_type()<<" :  ";
        node_id = node_id+1;
    }
    cout<<endl<<endl;
 }
}

void mesh_handler_HSCN::insert_iris(int thickness , float z_plane)
{
    if (2*thickness < Nyy)
    {
        cout<<" Inserting Iris......in ";
        cout<< " Plane : "<< z_plane <<endl;
          double dl = mesh_handler_HSCN_dl;
        int start_node = get_coordinate_iD_WG (mesh_handler_HSCN_width,mesh_handler_HSCN_height,mesh_handler_HSCN_length,npml,dl,0,0,z_plane); // on the left
        cout<<"start node bottom"<<start_node<<endl;
        int width_node = start_node + Nxx*Nyy-1; // on the right
        cout<<"start node top "<<width_node<<endl;

        for(int i = 0; i<thickness;i++)
        {
            int ii=start_node + i*Nxx;
            int jj=width_node - i*Nxx;
            for(int counter =0 ; counter<Nx; counter++)
            {
                mesh_handler_HSCN_nodes[ii].set_PEC();
                mesh_handler_HSCN_nodes[jj].set_PEC();
                ii = ii + 1;
                jj = jj - 1;
            }
        }
    }

    else cout<<"Cannot insert Iris "<<endl;
}

void mesh_handler_HSCN::insert_perfect_cube(int cube_l,int nx0,int ny0, int nz0)
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
                    mesh_handler_HSCN_nodes[loc2 + i].set_PEC();   //
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

void mesh_handler_HSCN::insert_dielectric_sheet(double sheet_thickness, double y_loc)
{
    int debug = 0;
    if( (y_loc > (mesh_handler_HSCN_height + mesh_handler_HSCN_dl*npml)) || ( y_loc + sheet_thickness > (mesh_handler_HSCN_height + mesh_handler_HSCN_dl*npml )) )
    {
        cout<<" FAILED TO INSERT "<<endl;

        return ;
    }


    int Ny_loc = int ( (y_loc / mesh_handler_HSCN_dl) + 0.5);
    int start_node = Ny_loc*Nxx;
    int node_id = 0;
 //cout<<debug;cin>>debug;

    for ( int i = 0; i< int(sheet_thickness/mesh_handler_HSCN_dl); i++)
        for( int j = 0; j<Nzz ; j++)
        {
            node_id = start_node + i*Nxx + Nxx*Nyy*j;

            for ( int k = 0 ; k< Nxx; k++)
            {

                this->mesh_handler_HSCN_nodes[node_id].set_material_type(2);
                //cout<<node_id<<" ";//cin>>debug;
                node_id = node_id + 1;
            }
           // cout<<endl;
        }
}


void mesh_handler_HSCN::breakpoint_(int i)
{
    cout<<" breakpoint " << i <<endl;
    cin>>i;
}

void mesh_handler_HSCN::print_this_neighbour(int excited_node)
{
    int ii = this->mesh_handler_HSCN_node_neighbours[excited_node].size();

    cout<<" Node : " << excited_node << endl;
    for( int i=1; i<=12; i++)  cout<<this->mesh_handler_HSCN_nodes[excited_node].get_reflected(i)<<endl;
    for (int j =0 ; j<ii ; j++)
    {
        cout<<" Neighbour "<< this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_iD() <<" : "<<endl;
        cout<<"1: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(1)<<endl;
        cout<<"2: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(2)<<endl;
        cout<<"3: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(3)<<endl;
        cout<<"4: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(4)<<endl;
        cout<<"5: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(5)<<endl;
        cout<<"6: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(6)<<endl;
        cout<<"7: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(7)<<endl;
        cout<<"8: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(8)<<endl;
        cout<<"9: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(9)<<endl;
        cout<<"10: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(10)<<endl;
        cout<<"11: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(11)<<endl;
        cout<<"12: "<<this->mesh_handler_HSCN_nodes[ mesh_handler_HSCN_node_neighbours[excited_node][j] ].get_incidence(12)<<" ? "<<endl;

        //cout<<this->mesh_handler_HSCN_node_neighbours[excited_node][j].get_reflected()<<endl;
    }
}

void mesh_handler_HSCN::TLM_excitation_te10(int node ,   double v_inc)
{
          double start_node = 0;
          double v = 0, v_ = 0;
        int i=0;
        for(int j = 0; j<Nyy; j++)
        {
            start_node = node + j*Nxx;

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(3);
            v_ = v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(3,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(11);
            v_ = v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(11,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(4);
            v_= v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(4,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(8);
            v_= v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(8,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(14);
            v_ = v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(14,v_);

        }

        //cout<<"here:"<<v_inc<<"node"<<endl<<mesh_handler_HSCN_nodes[28].get_incidence(3)<<endl;


}

void mesh_handler_HSCN::TLM_excitation_single_node(int node ,   double v_inc, int axis)
{
      double v = 0, v_ = 0;
    int start_node = node;

    //EXCITE EY field component
    if( axis ==2){
    v = mesh_handler_HSCN_nodes[start_node].get_incidence(3);
    v_ = v_inc + v;
    mesh_handler_HSCN_nodes[start_node].set_incidence(3,v_);

    v = mesh_handler_HSCN_nodes[start_node].get_incidence(11);
    v_ = v_inc + v;
    mesh_handler_HSCN_nodes[start_node].set_incidence(11,v_);

    v = mesh_handler_HSCN_nodes[start_node].get_incidence(4);
    v_= v_inc + v;
    mesh_handler_HSCN_nodes[start_node].set_incidence(4,v_);

    v = mesh_handler_HSCN_nodes[start_node].get_incidence(8);
    v_= v_inc + v;
    mesh_handler_HSCN_nodes[start_node].set_incidence(8,v_);

    v = mesh_handler_HSCN_nodes[start_node].get_incidence(14);
    v_ = v_inc + v;
    mesh_handler_HSCN_nodes[start_node].set_incidence(14,v_);

    }
    else if(axis == 3)
    {
        v = mesh_handler_HSCN_nodes[start_node].get_incidence(7);
        v_ = v_inc + v;
        mesh_handler_HSCN_nodes[start_node].set_incidence(7,v_);

        v = mesh_handler_HSCN_nodes[start_node].get_incidence(5);
        v_ = v_inc + v;
        mesh_handler_HSCN_nodes[start_node].set_incidence(5,v_);

        v = mesh_handler_HSCN_nodes[start_node].get_incidence(6);
        v_= v_inc + v;
        mesh_handler_HSCN_nodes[start_node].set_incidence(6,v_);

        v = mesh_handler_HSCN_nodes[start_node].get_incidence(10);
        v_= v_inc + v;
        mesh_handler_HSCN_nodes[start_node].set_incidence(10,v_);

        v = mesh_handler_HSCN_nodes[start_node].get_incidence(15);
        v_ = v_inc + v;
        mesh_handler_HSCN_nodes[start_node].set_incidence(15,v_);

        }
}

void mesh_handler_HSCN::TLM_excitation_plane_node(int node,   double v_inc, int exct_plane )
{

   int nxy = Nyy*Nxx;
     double v = 0, v_ = 0;
   int start_node = node;

    if( exct_plane ==1) //exciting xy plane i.e. z direction
    {

        for ( int j = 0; j < nxy ; j++){
            //EXCITE EY field component
            v = mesh_handler_HSCN_nodes[start_node].get_incidence(3);
            v_ = v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(3,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(11);
            v_ = v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(11,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(4);
            v_= v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(4,v_);

            v = mesh_handler_HSCN_nodes[start_node].get_incidence(8);
            v_= v_inc + v;
            mesh_handler_HSCN_nodes[start_node].set_incidence(8,v_);

                v = mesh_handler_HSCN_nodes[start_node].get_incidence(14);
                v_ = v_inc + v;
                mesh_handler_HSCN_nodes[start_node].set_incidence(14,v_);

            //cout<<"plane excitation "<<endl;
            start_node = start_node + 1;
        }
    }

}

void mesh_handler_HSCN::TLM_simulation1(int dt_total, int output_node,int excited_node,  double Zz1,   double Zz2,
                                    double Zy1,   double Zy2,  double Zx1,   double Zx2, f_gaussian func,   double exctn_type)
{
    int pause;
    int nxy = Nyy*Nxx;
      double dl = this->mesh_handler_HSCN_dl;
      double w = this->mesh_handler_HSCN_width;
      double h = this->mesh_handler_HSCN_height;
      double l = this->mesh_handler_HSCN_length;

// Guassian Pulse Parameters
      double f = 0;           //frequency
      double ff = 0;          //dummy variable
      double te = 0;          //dummy variable
      double t = timestep;

// Compute the start node in desired excitation plane
    int plane_start_node= excited_node;
// TE10 mode
    int te10_nd = excited_node;//get_coordinate_iD_WG ( w,h,l,npml,dl,0,0,dl );
    cout<< " excitation node " << excited_node<<endl;
    cout<<"output_node  "<< output_node <<endl;

// TIME LOOP
    for( int dt=0; dt<dt_total; dt++ )
    {
//EXCITATION
        if( (dt<2e3) )
        {
            // single node excitation
            if ( exctn_type==-1 )
            {
                if( !mesh_handler_HSCN_nodes[excited_node].is_PEC())
                {
                    f = func(dt)*(-1e-3*dl);
                    this->TLM_excitation_single_node(excited_node, f);
                }
            }
            // te10 excitation
            else if(exctn_type==-2)
            {
                int xy_end = te10_nd + Nxx;

                for (int i = te10_nd; i< xy_end ;i++ )
                {
                    f = sin((i-te10_nd)*pi/(Nxx-1))*func(dt)*1e-3*dl;
                    this->TLM_excitation_te10(i, f);
                }
            }
            //plane excitation
            else if ( exctn_type >=0 )
            {
               f = func(dt)*(-1e-3)*dl;
               this->TLM_excitation_plane_node(plane_start_node, f );
            }
        }

//OUTPUT FIELDS
    int ii= 0;


    if( mesh_handler_HSCN_nodes[output_node].is_PEC() )
        {
            cout<<" error!! "<<endl;
            cout<<output_node<<endl;
            cin>>pause;
        }

    E_output_at_node(output_node);

    // H_output_at_node(output_node);
    bool check = false;

//Scatter
    //Standard nodes
    #pragma omp parallel
    #pragma omp for
    for(int i=0; i< vctr_medium_iD[0].size(); i++)
    {
        this->mesh_handler_HSCN_nodes[this->vctr_medium_iD[0][i]].scattering_standard(line_Y,line_G,line_R);
    }
    #pragma omp parallel
    #pragma omp for
    for(int i=0; i< vctr_medium_iD[1].size(); i++)
    {
         this->mesh_handler_HSCN_nodes[this->vctr_medium_iD[1][i]].scattering_standard(line_Y2,line_G2,line_R2);
    }


    //PML nodes
    #pragma omp parallel
    #pragma omp for
    for(int i=0; i< vctr_PML_medium_iD[0].size(); i++)
    {
        this->mesh_handler_HSCN_nodes[this->vctr_PML_medium_iD[0][i]].scattering(line_Y,line_G,line_R);
    }

    #pragma omp parallel
    #pragma omp for
    for(int i=0; i< vctr_PML_medium_iD[1].size(); i++)
    {
        this->mesh_handler_HSCN_nodes[this->vctr_PML_medium_iD[1][i]].scattering(line_Y2,line_G2,line_R2);
    }
    //cout<<"s"<<endl<<endl;


//Connect
   if( sim_type == 4) // waveguide
    {
        if (Ny > 1 )this->TLM_connection_optimized(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2); //this->TLM_connection_sd(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);
        else if( Ny==1) this->TLM_connection_2D_HSCN_mesh_handler_HSCN(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);
    }
    // cubic domain
   else if( sim_type ==5 ) this->TLM_connection_optimized(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2); //this->TLM_connection_sd(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);

    //cout << "c"<<endl<<endl;

    //cout<< " dt = " <<dt <<" of "<<dt_total<<endl;


    }



}


void mesh_handler_HSCN::TLM_simulation_2D_HSCN(int dt_total, int output_node,int excited_node,  double Zz1,   double Zz2,   double Zy1,   double Zy2,  double Zx1,   double Zx2, f_gaussian func,   double z)
{

    int nxy = Nyy*Nxx;
      double dl = this->mesh_handler_HSCN_dl;
      double w = this->mesh_handler_HSCN_width;
      double h = this->mesh_handler_HSCN_height;
      double l = this->mesh_handler_HSCN_length;

    output_node = centre_node;
    excited_node  = centre_node;

// Guassian Pulse Parameters
      double f = 0;           //frequency
      double t = timestep;     //time factor

// Compute the start node in desired excitation plan

// TIME LOOP
    for( int dt=0; dt<dt_total; dt++ )
    {
        //EXCITATION
        if( (dt<300) )
        {
            // single node excitation
                if( !mesh_handler_HSCN_nodes[excited_node].is_PEC())
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

       E_output_at_node(output_node);
       H_output_at_node(output_node);

//cout<<Ez_output_at_node(output_node)<<endl;
        //Far Field Computation variables
        //vector<HSCN_node> v_plane_xy1,v_plane_xy2;

        bool check = false;

    #pragma omp parallel
    #pragma omp for
    for(int i=0; i< vctr_medium_iD[0].size(); i++)
    {
        this->mesh_handler_HSCN_nodes[this->vctr_medium_iD[0][i]].scattering(line_Y,line_G,line_R);
    }

    #pragma omp parallel
    #pragma omp for
    for(int i=0; i< vctr_medium_iD[1].size(); i++)
    {
         this->mesh_handler_HSCN_nodes[this->vctr_medium_iD[1][i]].scattering(line_Y2,line_G2,line_R2);
    }


    this->TLM_connection_2D_HSCN(Zz1,Zz2,Zy1,Zy2,Zx1,Zx2);

    }
}


void mesh_handler_HSCN::create_connect_bins()
{
     cout<<"Creating Bins for TLM Connect process ...."<<endl;
    // variable definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml; //0.5*(nnzz - nnz);
    int nxy = nny*nnx;
    double dl = this->mesh_handler_HSCN_dl;
    double w = this->mesh_handler_HSCN_width;
    double h = this->mesh_handler_HSCN_height;
    double l = this->mesh_handler_HSCN_length;

    connect_bins.resize(19);  //creates 19 different connect categories
    int ii=0;

    for(int z=0 ; z<this->Nzz; z++)
        {
            for (int y=0 ; y<this->Nyy ; y++)
            {
                for (int x=0 ; x<this->Nxx ; x++)
                {
                     // cout<<"here : "<<endl;
                    if( !mesh_handler_HSCN_nodes[ii].is_PEC() )
                    {
                        if(!mesh_handler_HSCN_nodes[ii].check_special_node()) connect_bins[0].push_back(&mesh_handler_HSCN_nodes[ii]);// STORE POINTER TO THE REGULAR NODES

                        else
                        {
                            if ( (x==0)||(x==nnx-1) )
                            {
                                if ( (x==0) ) connect_bins[2].push_back(&mesh_handler_HSCN_nodes[ii]);//left wall rule.   //Connect B

                                else  connect_bins[3].push_back(&mesh_handler_HSCN_nodes[ii]);  //right wall rule  //Connect C
                            }
                            else if ((mesh_handler_HSCN_nodes[ii+1].is_PEC())||(mesh_handler_HSCN_nodes[ii-1].is_PEC()) )
                            {
                                if (mesh_handler_HSCN_nodes[ii-1].is_PEC()) //left wall rule.
                                   {
                                    if (!mesh_handler_HSCN_nodes[ii+1].is_PEC()) connect_bins[4].push_back(&mesh_handler_HSCN_nodes[ii]); //Connect D

                                    else connect_bins[5].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect E
                                   }
                                else connect_bins[6].push_back(&mesh_handler_HSCN_nodes[ii]);// Connect F
                            }
                            else connect_bins[1].push_back(&mesh_handler_HSCN_nodes[ii]);////Connect A

                            if (( y==0) || (y==nny-1))
                            {
                                if (y==0)connect_bins[7].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect G   //bottom wall rule

                                else  connect_bins[8].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect H      //top wall rule
                            }

                            else if ((mesh_handler_HSCN_nodes[ii+nnx].is_PEC())||(mesh_handler_HSCN_nodes[ii-nnx].is_PEC()) )
                            {
                               if (mesh_handler_HSCN_nodes[ii-nnx].is_PEC())   //bottom wall rule
                                {
                                    if (!mesh_handler_HSCN_nodes[ii+nnx].is_PEC()) connect_bins[9].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect I

                                    else connect_bins[10].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect J
                                }
                                else        //top wall rule
                                {
                                    if (!mesh_handler_HSCN_nodes[ii-nnx].is_PEC())connect_bins[11].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect K
                                }
                            }
                            else  connect_bins[12].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect L     //reflect into adjacent nodes

                        if ( (z==0) || (z==nnzz-1))
                        {
                            if (z==0) connect_bins[13].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect M  // front face rule

                            else  connect_bins[14].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect N      //back face rule
                        }
                        else if ((mesh_handler_HSCN_nodes[ii+nxy].is_PEC())||(mesh_handler_HSCN_nodes[ii-nxy].is_PEC()))
                        {
                            if ((mesh_handler_HSCN_nodes[ii-nxy].is_PEC()))   // front face rule
                            {
                                if (!mesh_handler_HSCN_nodes[ii+nxy].is_PEC()) connect_bins[15].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect O

                                else connect_bins[16].push_back(&mesh_handler_HSCN_nodes[ii]);//Connect P
                            }
                            else  connect_bins[17].push_back(&mesh_handler_HSCN_nodes[ii]);//back face rule //Connect Q
                        }
                        else connect_bins[18].push_back(&mesh_handler_HSCN_nodes[ii]); //Connect R//normal reflection rules
                    }

                }
                ii=ii+1;
            }
        }
    }
    cout<<" TLM Connect bins created......"<<endl;
}

void mesh_handler_HSCN::TLM_connection_optimized(double Zz1 , double Zz2, double Zy1, double Zy2, double Zx1, double Zx2 )
{
   // variable definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml; //0.5*(nnzz - nnz);
    int nxy = nny*nnx;
    double dl = this->mesh_handler_HSCN_dl;
    double w = this->mesh_handler_HSCN_width;
    double h = this->mesh_handler_HSCN_height;
    double l = this->mesh_handler_HSCN_length;

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

            connect_bins[bin][i]->set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
            connect_bins[bin][i]->set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));

            connect_bins[bin][i]->set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
            connect_bins[bin][i]->set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
            connect_bins[bin][i]->set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
            connect_bins[bin][i]->set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));

            connect_bins[bin][i]->set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
            connect_bins[bin][i]->set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
            connect_bins[bin][i]->set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
            connect_bins[bin][i]->set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
        }
        //CONNECT FOR 1
         bin = bin+1;
        #pragma omp parallel
        #pragma omp for
         for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
            connect_bins[bin][i]->set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
        }

        //CONNECT FOR 2
         bin = bin+1;
         #pragma omp parallel
        #pragma omp for
         for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(6, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
            connect_bins[bin][i]->set_incidence(3, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
            connect_bins[bin][i]->set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
        }
        //CONNECT FOR 3
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
            connect_bins[bin][i]->set_incidence(10,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
            connect_bins[bin][i]->set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
        }

        //CONNECT FOR 4
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(6, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
            connect_bins[bin][i]->set_incidence(3, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
            connect_bins[bin][i]->set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
            connect_bins[bin][i]->set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
        }

        //CONNECT FOR 5
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(6, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
            connect_bins[bin][i]->set_incidence(3, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
            connect_bins[bin][i]->set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
            connect_bins[bin][i]->set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
        }

        //CONNECT FOR 6
         bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
            connect_bins[bin][i]->set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
            connect_bins[bin][i]->set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
            connect_bins[bin][i]->set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
        }

        //CONNECT FOR 7
         bin = bin+1;
         #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(1, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(1));
            connect_bins[bin][i]->set_incidence(5, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(5));
            connect_bins[bin][i]->set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
            connect_bins[bin][i]->set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));

        }

        //CONNECT FOR 8
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(12,reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(12));
            connect_bins[bin][i]->set_incidence(7, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(7));
            connect_bins[bin][i]->set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
            connect_bins[bin][i]->set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
        }
        //CONNECT FOR 9
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(1, -1*mesh_handler_HSCN_nodes[ii].get_reflected(1));
            connect_bins[bin][i]->set_incidence(5, -1*mesh_handler_HSCN_nodes[ii].get_reflected(5));
            connect_bins[bin][i]->set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
            connect_bins[bin][i]->set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
        }

        //CONNECT FOR 10
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            connect_bins[bin][i]->set_incidence(1, -1*mesh_handler_HSCN_nodes[ii].get_reflected(1));
            connect_bins[bin][i]->set_incidence(5, -1*mesh_handler_HSCN_nodes[ii].get_reflected(5));
            connect_bins[bin][i]->set_incidence(12,-1*mesh_handler_HSCN_nodes[ii].get_reflected(12));
            connect_bins[bin][i]->set_incidence(7, -1*mesh_handler_HSCN_nodes[ii].get_reflected(7));
        }
        //CONNECT FOR 11
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(12,-1*mesh_handler_HSCN_nodes[ii].get_reflected(12)); // It is possible to do it this way
            mesh_handler_HSCN_nodes[ii].set_incidence(7, -1*mesh_handler_HSCN_nodes[ii].get_reflected(7));
            mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
            mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
        }

        //CONNECT FOR 12
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
            mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
            mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
            mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
        }

        //CONNECT FOR 13
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(8, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
            mesh_handler_HSCN_nodes[ii].set_incidence(9, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
        }

        //CONNECT FOR 14
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(4,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(4));
            mesh_handler_HSCN_nodes[ii].set_incidence(2,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(2));
            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
        }

        //CONNECT FOR 15
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(8, -1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
            mesh_handler_HSCN_nodes[ii].set_incidence(9, -1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
        }

        //CONNECT FOR 16
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(8,-1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
            mesh_handler_HSCN_nodes[ii].set_incidence(9,-1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
            mesh_handler_HSCN_nodes[ii].set_incidence(4,-1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
            mesh_handler_HSCN_nodes[ii].set_incidence(2,-1*mesh_handler_HSCN_nodes[ii].get_reflected(2));
        }

        //CONNECT FOR 17
        bin = bin+1;
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();

            mesh_handler_HSCN_nodes[ii].set_incidence(4,-1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
            mesh_handler_HSCN_nodes[ii].set_incidence(2,-1*mesh_handler_HSCN_nodes[ii].get_reflected(2));
            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
        }

        //CONNECT FOR 18
        bin = bin+1;
        #pragma omp parallel
        #pragma omp for
        for(int i=0; i<connect_bins[bin].size();i++)
        {
            ii = connect_bins[bin][i]->get_iD();
            //Connect N
            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
        }

       if( bin!=18) cout<<"Error in optimized connect algorithm!!"<<endl;//cout<<bin<<endl;cin>>bin;
}

void mesh_handler_HSCN::TLM_scatter_optimized()
{
     vctr_medium_iD.resize(2);
     vctr_PML_medium_iD.resize(2);

    // SORT MEDIUM INTO DIFFERENT REGIONS
    for(int i=0; i< Ntotal_; i++)
    {
        if(!mesh_handler_HSCN_nodes[i].is_PEC())
            {
                if(mesh_handler_HSCN_nodes[i].get_medium_type() == 1)
                    {
                        if ( mesh_handler_HSCN_nodes[i].is_PML_HSCN() == 0)
                            vctr_medium_iD[0].push_back(    mesh_handler_HSCN_nodes[i].get_iD()   );

                        else
                            vctr_PML_medium_iD[0].push_back(   mesh_handler_HSCN_nodes[i].get_iD()   );
                    }

                else if(mesh_handler_HSCN_nodes[i].get_medium_type() == 2)
                    {
                        if ( mesh_handler_HSCN_nodes[i].is_PML_HSCN() == 0)
                            vctr_medium_iD[1].push_back(    mesh_handler_HSCN_nodes[i].get_iD()   );
                        else
                            vctr_PML_medium_iD[1].push_back(   mesh_handler_HSCN_nodes[i].get_iD()   );
                    }

            }
    }

}

void mesh_handler_HSCN::print_connect_bins()
{

    for(int i = 0;i<connect_bins.size();i++)
    {
        cout<<" Bin " <<i<<" : "<<connect_bins[i].size()<<endl;
    }
    //cout<<"centre_node  "<<centre_node<< "\\"<<connect_bins[10][0]->get_iD();
cout<<endl;
}

void mesh_handler_HSCN::TLM_connection(  double Zz1 ,   double Zz2,   double Zy1,   double Zy2,   double Zx1,   double Zx2 )
{
   // Variable Definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml;
    int nxy = nny*nnx;
      double dl = this->mesh_handler_HSCN_dl;
      double w = this->mesh_handler_HSCN_width;
      double h = this->mesh_handler_HSCN_height;
      double l = this->mesh_handler_HSCN_length;

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

                    if(!mesh_handler_HSCN_nodes[ii].check_special_node()){

                        mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                        mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                        mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                        mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));

                        mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));

                        mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                        mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                        mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                        mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));

                    }

                    else {
                    if ( (x==0)||(x==nnx-1) )
                    {

                            if ( (x==0) )//left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                            }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                            }
                    }
                    else
                    {
                        // reflect into adjacent nodes
                        mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                        mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                        mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                        mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                        //cout<<" Position : "<<" 5 "<<endl;
                    }

                    if (( y==0) || (y==nny-1))
                    {
                            if (y==0)   //bottom wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                                //cout<<" Position : "<<" 9 "<<endl;
                            }
                            else       //top wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                            }
                    }
                    else       //reflect into adjacent nodes
                    {
                        mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                        //cout<<" Position : "<<"12"<<endl;
                    }

                    if ( (z==0) || (z==nnzz-1))
                    {
                        if (z==0)   // front face rule
                        {
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
                            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                        }
                        else        //back face rule
                        {
                            mesh_handler_HSCN_nodes[ii].set_incidence(4,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                        }
                    }
                    else        //normal reflection rules
                    {
                        // reflect into adjacent nodes
                        mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                        mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                        mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                        mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                        //cout<<" Position : "<<"17"<<endl;
                    }
                }
                    //STUB CONNECTION
                        mesh_handler_HSCN_nodes[ii].set_incidence(13,mesh_handler_HSCN_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        mesh_handler_HSCN_nodes[ii].set_incidence(14,mesh_handler_HSCN_nodes[ii].get_reflected(14));
                        mesh_handler_HSCN_nodes[ii].set_incidence(15,mesh_handler_HSCN_nodes[ii].get_reflected(15));

                     ii=ii+1;
                    }
                }
            }
}


void mesh_handler_HSCN::TLM_connection_2D_HSCN_mesh_handler_HSCN(  double Zz1 ,   double Zz2,   double Zy1,   double Zy2,   double Zx1,   double Zx2 )
{
     if( Ny >1)
     {
        int a;
        cout<<" Error!! calling 2D HSCN connection for 3D HSCN node" <<endl;
        cin>>a;
     }

     // Variable Definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml;
    int nxy = nny*nnx;
      double dl = this->mesh_handler_HSCN_dl;
      double w = this->mesh_handler_HSCN_width;
      double h = this->mesh_handler_HSCN_height;
      double l = this->mesh_handler_HSCN_length;

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

                    if( !mesh_handler_HSCN_nodes[ii].is_PEC() )
                    {

                        if(!mesh_handler_HSCN_nodes[ii].check_special_node())
                        {

                            mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                            mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                            mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                            mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));

                            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));

                        }

                    else
                    {

                        if ( (x==0)||(x==nnx-1) )
                        {
                            if ( (x==0) )//left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                                }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                            }
                        }

                        else if ((mesh_handler_HSCN_nodes[ii+1].is_PEC())||(mesh_handler_HSCN_nodes[ii-1].is_PEC()) )
                        {
                            if (mesh_handler_HSCN_nodes[ii-1].is_PEC()) //left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                if (!mesh_handler_HSCN_nodes[ii+1].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                               // cout<<" Position : "<<" 6a"<<endl;
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                //cout<<" Position : "<<" 6 "<<endl;
                                }
                            }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                                //cout<<" Position : "<<" 6b"<<endl;
                            }
                        }
                        else
                        {
                            // reflect into adjacent nodes
                            mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                            mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                            mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                            mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                            //cout<<" Position : "<<" 5 "<<endl;
                        }

                        if ( (z==0) || (z==nnzz-1))
                        {
                            if (z==0)   // front face rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(8, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
                                mesh_handler_HSCN_nodes[ii].set_incidence(9, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
                                mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                                mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                            }
                            else        //back face rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(4,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                                mesh_handler_HSCN_nodes[ii].set_incidence(2,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                                mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                                mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                            }
                        }

                        else if ((mesh_handler_HSCN_nodes[ii+nxy].is_PEC())||(mesh_handler_HSCN_nodes[ii-nxy].is_PEC()))
                        {

                            if ((mesh_handler_HSCN_nodes[ii-nxy].is_PEC()))   // front face rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(8, -1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
                                mesh_handler_HSCN_nodes[ii].set_incidence(9, -1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
                                if (!mesh_handler_HSCN_nodes[ii+nxy].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                                mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                                //cout<<" Position : "<<"18a "<<endl;
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(4,-1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                                mesh_handler_HSCN_nodes[ii].set_incidence(2,-1*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                                //cout<<" Position : "<<"18"<<endl;
                                }
                            }
                            else       //back face rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(4,-1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                                mesh_handler_HSCN_nodes[ii].set_incidence(2,-1*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                                mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                                mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                                //cout<<" Position : "<<"18b"<<endl;
                            }
                        }
                        else        //normal reflection rules
                        {
                            // reflect into adjacent nodes
                            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                            //cout<<" Position : "<<"17"<<endl;
                        }
                        }
                     //Connection a  y axis
                        mesh_handler_HSCN_nodes[ii].set_incidence(12,-mesh_handler_HSCN_nodes[ii].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, -mesh_handler_HSCN_nodes[ii].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, -mesh_handler_HSCN_nodes[ii].get_reflected(7));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, -mesh_handler_HSCN_nodes[ii].get_reflected(5));

                    //STUB CONNECTION
                        mesh_handler_HSCN_nodes[ii].set_incidence(13,mesh_handler_HSCN_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        mesh_handler_HSCN_nodes[ii].set_incidence(14,mesh_handler_HSCN_nodes[ii].get_reflected(14));
                        mesh_handler_HSCN_nodes[ii].set_incidence(15,mesh_handler_HSCN_nodes[ii].get_reflected(15));
                    }
                    ii=ii+1;
                }
        }
        //cout<<"here"<<endl;
}

void mesh_handler_HSCN::TLM_connection_2D_HSCN(  double Zz1 ,   double Zz2,   double Zy1,   double Zy2,   double Zx1,   double Zx2 )
{
    // Variable Definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = this->npml;
    int nxy = nny*nnx;
      double dl = this->mesh_handler_HSCN_dl;
      double w = this->mesh_handler_HSCN_width;
      double h = this->mesh_handler_HSCN_height;
      double l = this->mesh_handler_HSCN_length;

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

                    if( !mesh_handler_HSCN_nodes[ii].is_PEC() ) {

                    if(!mesh_handler_HSCN_nodes[ii].check_special_node()){

                        mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                        mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                        mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                        mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));

                        mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));

                    }
                    else{

                    if ( (x==0)||(x==nnx-1) )
                    {

                            if ( (x==0) )//left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                            }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                            }
                    }

                    else if ((mesh_handler_HSCN_nodes[ii+1].is_PEC())||(mesh_handler_HSCN_nodes[ii-1].is_PEC()) )
                    {
                            if (mesh_handler_HSCN_nodes[ii-1].is_PEC()) //left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                if (!mesh_handler_HSCN_nodes[ii+1].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                               // cout<<" Position : "<<" 6a"<<endl;
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                //cout<<" Position : "<<" 6 "<<endl;
                                }
                            }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                                //cout<<" Position : "<<" 6b"<<endl;
                            }
                    }
                    else
                    {
                        // reflect into adjacent nodes
                        mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                        mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                        mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                        mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                        //cout<<" Position : "<<" 5 "<<endl;
                    }

                    if (( y==0) || (y==nny-1))
                    {
                            if (y==0)   //bottom wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                                //cout<<" Position : "<<" 9 "<<endl;
                            }
                            else       //top wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                            }
                    }

                    else if ((mesh_handler_HSCN_nodes[ii+nnx].is_PEC())||(mesh_handler_HSCN_nodes[ii-nnx].is_PEC()) )
                    {
                       if (mesh_handler_HSCN_nodes[ii-nnx].is_PEC())   //bottom wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, -1*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, -1*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                if (!mesh_handler_HSCN_nodes[ii+nnx].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,-1*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, -1*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                }
                            }
                        else        //top wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,-1*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, -1*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                if (!mesh_handler_HSCN_nodes[ii-nnx].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                                //cout<<" Position : "<<" 11b"<<endl;
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, -1*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, -1*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                //cout<<" Position : "<<"11"<<endl;

                                }
                            }
                    }

                    else       //reflect into adjacent nodes
                    {
                        mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                        //cout<<" Position : "<<"12"<<endl;
                    }

                    }
                    //connection in z axis
                        mesh_handler_HSCN_nodes[ii].set_incidence(8, -1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
                        mesh_handler_HSCN_nodes[ii].set_incidence(9, -1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
                        mesh_handler_HSCN_nodes[ii].set_incidence(4, -1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                        mesh_handler_HSCN_nodes[ii].set_incidence(2, -1*mesh_handler_HSCN_nodes[ii].get_reflected(2));

                    //STUB CONNECTION
                        mesh_handler_HSCN_nodes[ii].set_incidence(13,mesh_handler_HSCN_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        mesh_handler_HSCN_nodes[ii].set_incidence(14,mesh_handler_HSCN_nodes[ii].get_reflected(14));
                        mesh_handler_HSCN_nodes[ii].set_incidence(15,mesh_handler_HSCN_nodes[ii].get_reflected(15));

                    }
                     ii=ii+1;

                }
            }

}

void mesh_handler_HSCN::TLM_connection_sd(  double Zz1 ,   double Zz2,   double Zy1,   double Zy2,   double Zx1,   double Zx2 )
{

   // variable definitions
    int nny = this->Nyy;
    int nnx = this->Nxx;
    int nnz = this->Nz;
    int nnzz = this->Nzz;
    int n_PML = 0.5*(nnzz - nnz);
    int nxy = nny*nnx;
      double dl = this->mesh_handler_HSCN_dl;
      double w = this->mesh_handler_HSCN_width;
      double h = this->mesh_handler_HSCN_height;
      double l = this->mesh_handler_HSCN_length;

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

    for(int z=0 ; z<this->Nzz; z++)
        {
            for (int y=0 ; y<this->Nyy ; y++)
            {
                for (int x=0 ; x<this->Nxx ; x++)
                {
                    // cout<<"here : "<<endl;
                    if( !mesh_handler_HSCN_nodes[ii].is_PEC() ) {

                    if(!mesh_handler_HSCN_nodes[ii].check_special_node()){

                        mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                        mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                        mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                        mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));

                        mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));

                        mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                        mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                        mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                        mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));

                    }
                    else{

                    if ( (x==0)||(x==nnx-1) )
                    {

                            if ( (x==0) )//left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                            }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,reflct_x*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                            }
                    }
                    else if ((mesh_handler_HSCN_nodes[ii+1].is_PEC())||(mesh_handler_HSCN_nodes[ii-1].is_PEC()) )
                    {
                            if (mesh_handler_HSCN_nodes[ii-1].is_PEC()) //left wall rule.
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(6, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(6)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3, -1*(mesh_handler_HSCN_nodes[ii].get_reflected(3)));
                                if (!mesh_handler_HSCN_nodes[ii+1].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                               // cout<<" Position : "<<" 6a"<<endl;
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                //cout<<" Position : "<<" 6 "<<endl;
                                }
                            }
                            else    //right wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(11,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(11)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(10,-1*(mesh_handler_HSCN_nodes[ii].get_reflected(10)));
                                mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                                mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                                //cout<<" Position : "<<" 6b"<<endl;
                            }
                    }
                    else
                    {
                        // reflect into adjacent nodes
                        mesh_handler_HSCN_nodes[ii].set_incidence(11,mesh_handler_HSCN_nodes[ii+1].get_reflected(3));
                        mesh_handler_HSCN_nodes[ii].set_incidence(3,mesh_handler_HSCN_nodes[ii-1].get_reflected(11));
                        mesh_handler_HSCN_nodes[ii].set_incidence(10,mesh_handler_HSCN_nodes[ii+1].get_reflected(6));
                        mesh_handler_HSCN_nodes[ii].set_incidence(6,mesh_handler_HSCN_nodes[ii-1].get_reflected(10));
                        //cout<<" Position : "<<" 5 "<<endl;
                    }

                    if (( y==0) || (y==nny-1))
                    {
                            if (y==0)   //bottom wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                                //cout<<" Position : "<<" 9 "<<endl;
                            }
                            else       //top wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, reflct_y*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                            }
                    }

                    else if ((mesh_handler_HSCN_nodes[ii+nnx].is_PEC())||(mesh_handler_HSCN_nodes[ii-nnx].is_PEC()) )
                    {
                       if (mesh_handler_HSCN_nodes[ii-nnx].is_PEC())   //bottom wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, -1*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, -1*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                if (!mesh_handler_HSCN_nodes[ii+nnx].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,-1*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, -1*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                }
                            }
                        else        //top wall rule
                            {
                                mesh_handler_HSCN_nodes[ii].set_incidence(12,-1*mesh_handler_HSCN_nodes[ii].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(7, -1*mesh_handler_HSCN_nodes[ii].get_reflected(7));
                                if (!mesh_handler_HSCN_nodes[ii-nnx].is_PEC())
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                                //cout<<" Position : "<<" 11b"<<endl;
                                }
                                else
                                {
                                mesh_handler_HSCN_nodes[ii].set_incidence(1, -1*mesh_handler_HSCN_nodes[ii].get_reflected(1));
                                mesh_handler_HSCN_nodes[ii].set_incidence(5, -1*mesh_handler_HSCN_nodes[ii].get_reflected(5));
                                //cout<<" Position : "<<"11"<<endl;

                                }
                            }
                    }

                    else       //reflect into adjacent nodes
                    {
                        mesh_handler_HSCN_nodes[ii].set_incidence(12,mesh_handler_HSCN_nodes[ii+nnx].get_reflected(1));
                        mesh_handler_HSCN_nodes[ii].set_incidence(1, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(7, mesh_handler_HSCN_nodes[ii+nnx].get_reflected(5));
                        mesh_handler_HSCN_nodes[ii].set_incidence(5, mesh_handler_HSCN_nodes[ii-nnx].get_reflected(7));
                        //cout<<" Position : "<<"12"<<endl;
                    }

                    if ( (z==0) || (z==nnzz-1))
                    {
                        if (z==0)   // front face rule
                        {
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, reflct_z1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
                            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                        }
                        else        //back face rule
                        {
                            mesh_handler_HSCN_nodes[ii].set_incidence(4,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2,reflct_z2*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                        }
                    }
                    else if ((mesh_handler_HSCN_nodes[ii+nxy].is_PEC())||(mesh_handler_HSCN_nodes[ii-nxy].is_PEC()))
                    {
                        if ((mesh_handler_HSCN_nodes[ii-nxy].is_PEC()))   // front face rule
                        {
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, -1*mesh_handler_HSCN_nodes[ii].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, -1*mesh_handler_HSCN_nodes[ii].get_reflected(9));
                            if (!mesh_handler_HSCN_nodes[ii+nxy].is_PEC())
                            {
                            mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                            //cout<<" Position : "<<"18a "<<endl;
                            }
                            else
                            {
                            mesh_handler_HSCN_nodes[ii].set_incidence(4,-1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2,-1*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                            //cout<<" Position : "<<"18"<<endl;
                            }
                        }
                        else       //back face rule
                        {
                            mesh_handler_HSCN_nodes[ii].set_incidence(4,-1*mesh_handler_HSCN_nodes[ii].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(2,-1*mesh_handler_HSCN_nodes[ii].get_reflected(2));
                            mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                            mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                            //cout<<" Position : "<<"18b"<<endl;
                        }
                    }
                    else        //normal reflection rules
                    {
                        // reflect into adjacent nodes
                        mesh_handler_HSCN_nodes[ii].set_incidence(4, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(8));
                        mesh_handler_HSCN_nodes[ii].set_incidence(8, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(4));
                        mesh_handler_HSCN_nodes[ii].set_incidence(2, mesh_handler_HSCN_nodes[ii+nxy].get_reflected(9));
                        mesh_handler_HSCN_nodes[ii].set_incidence(9, mesh_handler_HSCN_nodes[ii-nxy].get_reflected(2));
                        //cout<<" Position : "<<"17"<<endl;

                    }
                }
                    //STUB CONNECTION
                    //set_incidence(12,get_reflected(12));
                        mesh_handler_HSCN_nodes[ii].set_incidence(13,mesh_handler_HSCN_nodes[ii].get_reflected(13));            // capacitive stubs become  open circuits
                        mesh_handler_HSCN_nodes[ii].set_incidence(14,mesh_handler_HSCN_nodes[ii].get_reflected(14));
                        mesh_handler_HSCN_nodes[ii].set_incidence(15,mesh_handler_HSCN_nodes[ii].get_reflected(15));
                    }
                    ii=ii+1;
                }
            }
        }
}

void mesh_handler_HSCN:: write_output_file (const string& nodefileEX,const string& nodefileEY, const string& nodefileEZ, const string& nodefileHX,const string& nodefileHY,const string& nodefileHZ)                                           //writes to node file
{
    ofstream out_fileEY;
    ofstream out_fileEX;
    ofstream out_fileEZ;

    ofstream out_fileHY;
    ofstream out_fileHX;
    ofstream out_fileHZ;

    out_fileEX.open((nodefileEX+"").c_str());
    out_fileEZ.open((nodefileEZ+"").c_str());
    out_fileEY.open((nodefileEY+"").c_str());

    out_fileHX.open((nodefileHX+"").c_str());
    out_fileHY.open((nodefileHY+"").c_str());
    out_fileHZ.open((nodefileHZ+"").c_str());

    int end__ = this->mesh_handler_HSCN_Ex.size();
    int end_= this->mesh_handler_HSCN_Ey.size();

    out_fileEY<<setprecision(4)<<this->mesh_handler_HSCN_dl<<endl;
    out_fileEX<<setprecision(4)<<this->mesh_handler_HSCN_dl<<endl;
    out_fileEZ<<setprecision(4)<<this->mesh_handler_HSCN_dl<<endl;
    out_fileHY<<setprecision(4)<<this->mesh_handler_HSCN_dl<<endl;
    out_fileHX<<setprecision(4)<<this->mesh_handler_HSCN_dl<<endl;
    out_fileHZ<<setprecision(4)<<this->mesh_handler_HSCN_dl<<endl;

    out_fileEY<<setprecision(4)<<end_<<endl;
    out_fileEX<<setprecision(4)<<end_<<endl;
    out_fileEZ<<setprecision(4)<<end_<<endl;
    out_fileHY<<setprecision(4)<<end_<<endl;
    out_fileHX<<setprecision(4)<<end_<<endl;
    out_fileHZ<<setprecision(4)<<end_<<endl;

    cout<< " Writing to file...." << endl;
    for(int ii = 0 ; ii< end_; ii++)
    {
     out_fileEY<<setprecision(numeric_limits<  double>::digits10+1) <<mesh_handler_HSCN_Ey[ii]<<endl;
     //cout<<mesh_handler_HSCN_Ey[ii]<<endl;
     //out_fileHX<<setprecision(numeric_limits<  double>::digits10+1)<<mesh_handler_HSCN_Hx[ii]<<endl;
     out_fileEX<<setprecision(numeric_limits<  double>::digits10+1)<<mesh_handler_HSCN_Ex[ii]<<endl;
    // out_fileHY<<setprecision(numeric_limits<  double>::digits10+1)<<mesh_handler_HSCN_Hy[ii]<<endl;
     out_fileEZ<<setprecision(numeric_limits<  double>::digits10+1)<<mesh_handler_HSCN_Ez[ii]<<endl;
     //out_fileHZ<<setprecision(numeric_limits<  double>::digits10+1)<<mesh_handler_HSCN_Hz[ii]<<endl;
    }
}

void mesh_handler_HSCN:: write_output_file_for_yz_plane (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->mesh_handler_HSCN_Ey_plane_yz[0].size();

    cout<< " Writing yz plane to file...." << endl;

    out_fileEY<<this->Nyy<<","<<this->Nzz<<","<<mesh_handler_HSCN_Ey_plane_yz.size()<<endl;

    for (int dt=0; dt< mesh_handler_HSCN_Ey_plane_yz.size() ; dt++){
        for(int ii = 0 ; ii< end_; ii++)
        {
            out_fileEY<<setprecision(numeric_limits<  double>::digits10+1) <<mesh_handler_HSCN_Ey_plane_yz[dt][ii]<<",";
        }
        out_fileEY<<endl;
    }
}

void mesh_handler_HSCN:: write_output_file_for_xy_plane (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->mesh_handler_HSCN_Ey_plane_xy[0].size();

    cout<< " Writing xy plane to file...." << endl;

    out_fileEY<<this->Nxx<<","<<this->Nyy<<","<<mesh_handler_HSCN_Ey_plane_xy.size()<<endl;

    for (int dt=0; dt< mesh_handler_HSCN_Ey_plane_xy.size() ; dt++){
        for(int ii = 0 ; ii< end_; ii++)
        {
            out_fileEY<<setprecision(numeric_limits<  double>::digits10+1) <<mesh_handler_HSCN_Ey_plane_xy[dt][ii]<<",";
        }
        out_fileEY<<endl;
    }
}

void mesh_handler_HSCN:: write_output_file_for_zx_plane (const string& nodefileEY)                                           //writes to node file
{
    ofstream out_fileEY;

    out_fileEY.open((nodefileEY+"").c_str());

    int end_ = this->mesh_handler_HSCN_Ey_plane_zx[0].size();

    cout<< " Writing zx plane to file...." << endl;

    out_fileEY<<this->Nxx<<","<<this->Nzz<<","<<mesh_handler_HSCN_Ey_plane_zx.size()<<endl;

    for (int dt=0; dt< mesh_handler_HSCN_Ey_plane_zx.size() ; dt++){
        for(int ii = 0 ; ii< end_; ii++)
        {
            out_fileEY<<setprecision(numeric_limits<  double>::digits10+1) <<mesh_handler_HSCN_Ey_plane_zx[dt][ii]<<",";
        }
        out_fileEY<<endl;
    }
}


void mesh_handler_HSCN::print_Ey_output()
{
    cout<< " PRINTING Ey......." <<endl;
    for(int i= 0; i<mesh_handler_HSCN_Ey.size(); i++)
        cout<<i << ": "<< setprecision(15)<<mesh_handler_HSCN_Ey[i] <<endl;
}

void mesh_handler_HSCN::print_Ez_output()
{
    cout<< " PRINTING Ez......." <<endl;
    for(int i= 0; i<mesh_handler_HSCN_Ez.size(); i++)
        cout<<i << ": "<< setprecision(15)<<mesh_handler_HSCN_Ez[i] <<endl;
}

bool mesh_handler_HSCN::print_G_PML_parameters(int L_R,int PML_layer)
{

}

vector<  double> mesh_handler_HSCN::compute_far_field(const vector<int> &plane_xy1)
{
    /*
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
       x_ = this->mesh_handler_HSCN_nodes[ plane_xy1[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->mesh_handler_HSCN_nodes[ plane_xy1[i] ].get_coord(2) - this->npml;
       z_ = this->mesh_handler_HSCN_nodes[ plane_xy1[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_xy1_i = A_xy1_i + J_xy1_i *G_;
       A_xy1_j = A_xy1_j + J_xy1_j *G_;
       A_xy1_k = A_xy1_k + J_xy1_k *G_;

       F_xy1_i = F_xy1_i + M_xy1_i *G_;
       F_xy1_j = F_xy1_j + M_xy1_j *G_;
       F_xy1_k = F_xy1_k + M_xy1_k *G_;


       //xy2
       x_ = this->mesh_handler_HSCN_nodes[ plane_xy2[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->mesh_handler_HSCN_nodes[ plane_xy2[i] ].get_coord(2) - this->npml;
       z_ = this->mesh_handler_HSCN_nodes[ plane_xy2[i] ].get_coord(3) - this->npml;
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
       x_ = this->mesh_handler_HSCN_nodes[ plane_yz1[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->mesh_handler_HSCN_nodes[ plane_yz1[i] ].get_coord(2) - this->npml;
       z_ = this->mesh_handler_HSCN_nodes[ plane_yz1[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_yz1_i = A_yz1_i + J_yz1_i *G_;
       A_yz1_j = A_yz1_j + J_yz1_j *G_;
       A_yz1_k = A_yz1_k + J_yz1_k *G_;

       F_yz1_i = F_yz1_i + M_yz1_i *G_;
       F_yz1_j = F_yz1_j + M_yz1_j *G_;
       F_yz1_k = F_yz1_k + M_yz1_k *G_;

       //yz2
       x_ = this->mesh_handler_HSCN_nodes[ plane_yz2[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->mesh_handler_HSCN_nodes[ plane_yz2[i] ].get_coord(2) - this->npml;
       z_ = this->mesh_handler_HSCN_nodes[ plane_yz2[i] ].get_coord(3) - this->npml;
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
       x_ = this->mesh_handler_HSCN_nodes[ plane_zx1[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->mesh_handler_HSCN_nodes[ plane_zx1[i] ].get_coord(2) - this->npml;
       z_ = this->mesh_handler_HSCN_nodes[ plane_zx1[i] ].get_coord(3) - this->npml;
       r_abs  = sqrt( (x_ - x )*(x_ - x) + (y_ - y )*(y_ - y) + (z_ - z )*(z_ - z) ) ;  // absolute value
       G_= cos(r_abs*k0)/r_abs  // surface multiplier
       A_zx1_i = A_zx1_i + J_zx1_i *G_;
       A_zx1_j = A_zx1_j + J_zx1_j *G_;
       A_zx1_k = A_zx1_k + J_zx1_k *G_;

       F_zx1_i = F_zx1_i + M_zx1_i *G_;
       F_zx1_j = F_zx1_j + M_zx1_j *G_;
       F_zx1_k = F_zx1_k + M_zx1_k *G_;


       //zx2
       x_ = this->mesh_handler_HSCN_nodes[ plane_zx2[i] ].get_coord(1) - this->npml;        // coordinates on the plane_xy1
       y_ = this->mesh_handler_HSCN_nodes[ plane_zx2[i] ].get_coord(2) - this->npml;
       z_ = this->mesh_handler_HSCN_nodes[ plane_zx2[i] ].get_coord(3) - this->npml;
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




