

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

#include "mesh_handler_2D.h"
#include "shunt_node.h"
#include "gemini_utility.h"

const double pi = 3.14159265;
const double  u0 = 12566370614e-16;
const double e0 = 88541878176e-22;
const double c = 299792458.;
double Z0 = sqrt(u0/e0);//
double Z0_shunt = Z0*sqrt(2);

using namespace std;

mesh_handler_2D::mesh_handler_2D( bool freespace,float er, float sigma_e, float h, float w,float dl_,int PMLx, int PMLy, double Refn_fctr, int cndct_prof)
{

// SETTING MEMBER VARIABLES
    height = h;
    width = w;
    dl = dl_;
    Ny = int((  h/dl    )+  0.5);
    Nx = int((  w/dl    )+  0.5);
    Nxy = Ny*Nx;
    npmlx = PMLx;
    npmly = PMLy;
    Nxx  = Nx + 2*PMLx;
    Nyy  = Ny + 2*PMLy;
    Nxxyy = Nxx * Nyy;
    nodetype = 0;                       //not set

    long double dt = dl*1e-3/(c*sqrt(2));

// DEFINING PML
    int PML_boundaryx1 = PMLx;
    int PML_boundaryx2 = Nx + PMLx;
    int PML_boundaryy1 = PMLy;
    int PML_boundaryy2 = Ny + PMLy;

    //float g = 3.9;
    double sigmax = -(   e0*0.5*c *log(Refn_fctr)   )/(1e-3*dl_*pow(npmlx,cndct_prof+1)    );
    double sigmay = -(   e0*0.5*c *log(Refn_fctr)   )/(1e-3*dl_*pow(npmly,cndct_prof+1)    );

    int id=0;               //id

    for (int y=0 ; y<Nyy ; y++)
    {
        for (int x=0 ; x<Nxx ; x++)
        {
            // Creating Shunt_nodes;
            if( ((x >= PML_boundaryx1 )&& (x < PML_boundaryx2 ))&& ((y >= PML_boundaryy1 )&& (y < PML_boundaryy2 )) )
            {
                prob_domain2D.push_back(shunt_node(freespace,er,sigma_e,id,dt,x,y,dl));
            }
            else
            {
                int Lx = -1;
                int Ly = -1;

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

                prob_domain2D.push_back(    shunt_node(freespace,er,sigma_e,id,dt,x,y,dl,true,Lx,Ly,sigmax,sigmay,cndct_prof) );
            }

            id = id + 1;
        }
    }
}

mesh_handler_2D::mesh_handler_2D(const mesh_handler_2D& copy_mesh_handler_2D)
{
    height = copy_mesh_handler_2D.height;
    width = copy_mesh_handler_2D.width;
    dl = copy_mesh_handler_2D.dl;
    Ny = copy_mesh_handler_2D.Ny;
    Nx = copy_mesh_handler_2D.Nx;
    Nxy = copy_mesh_handler_2D.Nxy;
    npmlx = copy_mesh_handler_2D.npmlx;
    npmly = copy_mesh_handler_2D.npmly;
    Nxx = copy_mesh_handler_2D.Nxx;
    Nyy = copy_mesh_handler_2D.Nyy;
    Nxxyy = copy_mesh_handler_2D.Nxxyy;
    nodetype = copy_mesh_handler_2D.nodetype;
    Ez_plane_xy = copy_mesh_handler_2D.Ez_plane_xy;
    Ez_output = copy_mesh_handler_2D.Ez_output;
    Ez_far_field = copy_mesh_handler_2D.Ez_far_field;
    Hx_output = copy_mesh_handler_2D.Hx_output;
    Is = copy_mesh_handler_2D.Is;
    Vs = copy_mesh_handler_2D.Vs;
}



mesh_handler_2D &mesh_handler_2D::operator=(const mesh_handler_2D& c_mesh_handler_2D)
{
    if (this!= &c_mesh_handler_2D)
    {
        height = c_mesh_handler_2D.height;
        width = c_mesh_handler_2D.width;
        dl = c_mesh_handler_2D.dl;
        Ny = c_mesh_handler_2D.Ny;
        Nx = c_mesh_handler_2D.Nx;
        Nxy = c_mesh_handler_2D.Nxy;
        npmlx = c_mesh_handler_2D.npmlx;
        npmly = c_mesh_handler_2D.npmly;
        Nxx = c_mesh_handler_2D.Nxx;
        Nyy = c_mesh_handler_2D.Nyy;
        Nxxyy = c_mesh_handler_2D.Nxxyy;
        nodetype = c_mesh_handler_2D.nodetype;
        Ez_plane_xy = c_mesh_handler_2D.Ez_plane_xy;
        Ez_output = c_mesh_handler_2D.Ez_output;
        Ez_far_field = c_mesh_handler_2D.Ez_far_field;
        Hx_output = c_mesh_handler_2D.Hx_output;
        Is = c_mesh_handler_2D.Is;
        Vs = c_mesh_handler_2D.Vs;
    }
    return (*this);
}

void mesh_handler_2D::Simulate_2D(int ttltimestep, vector<int> &bdry_cdn, int excit_iD, int excit_lngth, int obsrv_iD,f_gaussian func, int exct_type, int d, int source_res)
{

    Ez_output.reserve(ttltimestep);
    float Vz = 0, Vz_prime=0;
    double Isrc=0;
    double sigma (10e8), f0(6e9);
    double dt (sqrt(2)*1e-3*dl/c);
    double f (0);

    //double start_x = abs(width-1000)/2 - d*dl;           //d nodes from scatterer....needed for plane wave propagation and computing reflection coefficient

//Setting the Excitation in L
    double te10_nd =excit_iD;

// TLM ALGORITHM
    for(int t=0; t<ttltimestep; t++)
    {
        if( t < excit_lngth )
        {
            if(exct_type ==-1) // point source excitation
            {
                f = func(t)*(-1e-3*dl);
                if( source_res>0)
                {
                    this->prob_domain2D[excit_iD].compute_source_voltages(f,source_res,Vz,Vz_prime);
                    Isrc = (Vz_prime - f)/source_res;
                }

                this->Geometry_excitation_Ez(f,excit_iD,source_res);

            }

            else if( exct_type ==0)      //line source - naca0015
            {
                for ( int j=0; j<Ny; j++)
                {
                    f = func(t)*(-1e-3*dl);
                    this->Geometry_excitation_Ez(f,excit_iD+j*Nxx,source_res);
                }
            }

            else if (exct_type ==-2)      //TE10 - wr28
            {
                for ( int j = te10_nd; j< Nx+te10_nd ; j++)
                {
                    f = sin((j-te10_nd)*pi/(Nx-1))*func(t)*(-1e-3*dl);
                    this->Geometry_excitation_Ez(f,j,source_res);
                }
            }

        }

//TLM SCATTER
        this->TLM_scatter();

        int aa;
// OUTPUT FIELDS
        Ez_output.push_back(  this->Geometry_observation_Ez(obsrv_iD) );     //observation must be had after the scattering because v_total is calculated in scatter member function
        Hx_output.push_back(  this->Geometry_observation_Hx(obsrv_iD) );

        //Ez_far_field.push_back(far_field_RCS(3,100000));
        //cin>>aa;

        Is.push_back(Isrc);
        Vs.push_back(f);

//OUTPUT FIELDS ACROSS PLANE
        if( t % 30 == 0)
        {
            vector< double > xy_plane;
            for(int i=0; i<this->Nxxyy ; i++) xy_plane.push_back(  this->Geometry_observation_Ez(i) );     //observation must be had after the scattering because v_total is calculated in scatter member function

            //store field value for each node in mesh plane
            Ez_plane_xy.push_back(xy_plane);

        }

// TLM CONNECT
        this->TLM_connect(bdry_cdn);
    }
}

void mesh_handler_2D:: write_output_file_2D (const string& nodefileEZ, const string& nodefileHX,const string& nodefileIs,const string& nodefileVs)
{

    ofstream out_fileEZ;
    ofstream out_fileHX;
    ofstream out_fileIs;
    ofstream out_fileVs;
    ofstream out_fileEz_far;

    out_fileEZ.open((nodefileEZ+"").c_str());
    out_fileHX.open((nodefileHX+"").c_str());
    out_fileIs.open((nodefileIs+"").c_str());
    out_fileVs.open((nodefileVs+"").c_str());

    string s = nodefileEZ.substr(0,nodefileEZ.length() - 4);

    out_fileEz_far.open((s+"_far_field_.txt").c_str());

    int end_ = this->Ez_output.size();

    out_fileEZ<<setprecision(4)<<this->dl<<endl;
    out_fileEZ<<setprecision(4)<<end_<<endl;

    out_fileHX<<setprecision(4)<<this->dl<<endl;
    out_fileHX<<setprecision(4)<<end_<<endl;

    out_fileIs<<setprecision(4)<<this->dl<<endl;
    out_fileIs<<setprecision(4)<<end_<<endl;

    out_fileVs<<setprecision(4)<<this->dl<<endl;
    out_fileVs<<setprecision(4)<<end_<<endl;

    //out_fileEz_far<<setprecision(4)<<this->dl<<endl;
    //out_fileEz_far<<setprecision(4)<<end_<<endl;

    cout<< " Writing to file...." << endl;

    for(int ii = 0 ; ii< end_; ii++)
    {
        out_fileEZ<<setprecision(numeric_limits<double>::digits10+1) <<Ez_output[ii]<<endl;
        out_fileHX<<setprecision(numeric_limits<double>::digits10+1) <<Hx_output[ii]<<endl;
        out_fileIs<<setprecision(numeric_limits<double>::digits10+1) <<Is[ii]<<endl;
        out_fileVs<<setprecision(numeric_limits<double>::digits10+1) <<Vs[ii]<<endl;
        //out_fileEz_far<<setprecision(numeric_limits<double>::digits10+1) <<Ez_far_field[ii]<<endl;
    }

}

void mesh_handler_2D::write_output_file_for_xy_plane_shuntnode(const string& nodefileEZ)
{

    ofstream out_fileEZ;

    out_fileEZ.open((nodefileEZ+"").c_str());

    int end_ = this->Ez_plane_xy[0].size();
    cout<< end_<<endl;

    int a = 0;
    //cin>>a;

    cout<< " Writing xy plane to file...." << endl;

    out_fileEZ<<this->Nxx<<","<<this->Nyy<<","<<Ez_plane_xy.size()<<endl;
    /*
        for (int t=0; t< Ez_plane_xy.size() ; t++){
            for(int ii = 0 ; ii< end_; ii++)
            {
                out_fileEZ<<setprecision(numeric_limits<double>::digits10+1) <<Ez_plane_xy[t][ii]<<",";
            }
            out_fileEZ<<endl;
        }

    */

    int total_time = Ez_plane_xy.size();

    for(int ii = 0 ; ii< end_; ii++)
    {
        for (int t=0; t< total_time ; t++)
        {
            out_fileEZ<<Ez_plane_xy[t][ii]<<",";
        }
        out_fileEZ<<endl;
    }


}

void mesh_handler_2D::Geometry_excitation_Ez(float ez_inc, int excit_iD, int source_res)
{

    int a;
    float Vz = 0,  Vz_prime = 0;
    float delta_Vz = 0;
    // Check if excitation node is PML / PEC
    if((this->prob_domain2D[excit_iD].is_PML_shunt( ) )||(this->prob_domain2D[excit_iD].is_PEC( ) ))
    {
        cout<<" excitation iD : "<< excit_iD<<endl;
        if(this->prob_domain2D[excit_iD].is_PML_shunt( ))
        {
            cout<<" PML "<<endl;
        }
        else
        {
            cout<<" PEC "<<endl;
        }
        cout<<" exciting in pml layer OR Exciting PEC...."<<endl;
        cin>>a;
    }

    if( source_res > 0)
    {
        this->prob_domain2D[excit_iD].compute_source_voltages(ez_inc,source_res,Vz,Vz_prime);
        delta_Vz = Vz_prime - Vz;
        ez_inc = delta_Vz;
    }

    this->prob_domain2D[excit_iD].shunt_excitation(ez_inc );

}

double mesh_handler_2D::Geometry_observation_Ez(int obsrv)
{
    //cout<<"observing at "<< obsrv<<endl;
    return (this->prob_domain2D[obsrv].get_Ez());
}

double mesh_handler_2D::Geometry_observation_Hx(int obsrv)
{
    //cout<<"observing at "<< obsrv<<endl;
    return (this->prob_domain2D[obsrv].get_Hx());
}

void mesh_handler_2D::TLM_scatter()
{
    for(int i=0 ; i<this->prob_domain2D.size(); i++)
        if(!this->prob_domain2D[i].is_PEC()) this->prob_domain2D[i].shunt_scatter();
}

void mesh_handler_2D::TLM_connect(vector<int> &bdry_cdn)
{
    int i(0);
    double bdry_cdntn1 = bdry_cdn[0];
    double bdry_cdntn2 = bdry_cdn[1];


    for (int y=0 ; y<this->Nyy ; y++)
        for (int x=0 ; x<this->Nxx ; x++)
        {
            double sigmax = this->prob_domain2D[i].get_cnductvty_x();
            double sigmay = this->prob_domain2D[i].get_cnductvty_y();
            double dt_ = this->prob_domain2D[i].get_dt();
            double dt_scale_x = 1;//exp(-sigmax*dt_);
            double dt_scale_y = 1;//exp(-sigmay*dt_);

            if(!this->prob_domain2D[i].is_PEC())
            {

                if ( (x==0)||(x==this->Nxx-1))
                {
                    if ( (x==0) )//left wall rule.
                    {
                        this->prob_domain2D[i].set_incidence( 1, dt_scale_x*bdry_cdntn1*this->prob_domain2D[i].get_reflected(1) );
                        this->prob_domain2D[i].set_incidence(3, dt_scale_x*this->prob_domain2D[i+1].get_reflected(1) );
                    }
                    else    //right wall rule
                    {
                        this->prob_domain2D[i].set_incidence(3, dt_scale_x*bdry_cdntn1*this->prob_domain2D[i].get_reflected(3) );
                        this->prob_domain2D[i].set_incidence(1, dt_scale_x*this->prob_domain2D[i-1].get_reflected(3)  );
                    }

                }
                else if ((this->prob_domain2D[i-1].is_PEC()) ||(this->prob_domain2D[i+1].is_PEC()))
                {
                    if( this->prob_domain2D[i-1].is_PEC())
                    {
                        this->prob_domain2D[i].set_incidence( 1, -1*dt_scale_x*this->prob_domain2D[i].get_reflected(1) );
                        this->prob_domain2D[i].set_incidence(3, dt_scale_x*this->prob_domain2D[i+1].get_reflected(1) );
                    }
                    else if ( this->prob_domain2D[i+1].is_PEC())
                    {
                        this->prob_domain2D[i].set_incidence(3, -1*dt_scale_x*this->prob_domain2D[i].get_reflected(3) );
                        this->prob_domain2D[i].set_incidence(1, dt_scale_x*this->prob_domain2D[i-1].get_reflected(3)  );

                    }
                }

                else
                {
                    this->prob_domain2D[i].set_incidence( 1, dt_scale_x*this->prob_domain2D[i-1].get_reflected(3) );
                    this->prob_domain2D[i].set_incidence( 3,dt_scale_x*this->prob_domain2D[i+1].get_reflected(1) );
                }

                if ( (y==0)||(y==this->Nyy-1) )
                {
                    if ( (y==0) )//bottom wall rule.
                    {
                        this->prob_domain2D[i].set_incidence( 2, dt_scale_y*bdry_cdntn2*this->prob_domain2D[i].get_reflected(2) );
                        this->prob_domain2D[i].set_incidence(0, dt_scale_y*this->prob_domain2D[i+Nxx].get_reflected(2) );
                    }
                    else    //right wall rule
                    {
                        this->prob_domain2D[i].set_incidence(0, dt_scale_y*bdry_cdntn2*this->prob_domain2D[i].get_reflected(0) );
                        this->prob_domain2D[i].set_incidence(2, dt_scale_y*this->prob_domain2D[i-this->Nxx].get_reflected(0) );

                    }
                }

                else if ((this->prob_domain2D[i+Nxx].is_PEC()) || (this->prob_domain2D[i-Nxx].is_PEC()))
                {

                    if((this->prob_domain2D[i-Nxx].is_PEC()) &&(this->prob_domain2D[i+Nxx].is_PEC()))
                    {
                        this->prob_domain2D[i].set_incidence( 2, -1*dt_scale_y*this->prob_domain2D[i].get_reflected(2) );
                        this->prob_domain2D[i].set_incidence(0, -1*dt_scale_y*this->prob_domain2D[i].get_reflected(0) );
                    }
                    else
                    {
                        if (this->prob_domain2D[i-Nxx].is_PEC())
                        {
                            this->prob_domain2D[i].set_incidence( 2, -1*dt_scale_y*this->prob_domain2D[i].get_reflected(2) );
                            this->prob_domain2D[i].set_incidence(0, dt_scale_y*this->prob_domain2D[i+Nxx].get_reflected(2) );
                        }
                        else if (this->prob_domain2D[i+Nxx].is_PEC())
                        {
                            this->prob_domain2D[i].set_incidence(0, -1*dt_scale_y*this->prob_domain2D[i].get_reflected(0) );
                            this->prob_domain2D[i].set_incidence(2, dt_scale_y*this->prob_domain2D[i-this->Nxx].get_reflected(0) );
                        }
                    }
                }
                else
                {
                    this->prob_domain2D[i].set_incidence(0, dt_scale_y*this->prob_domain2D[i+this->Nxx].get_reflected(2) );
                    this->prob_domain2D[i].set_incidence(2,dt_scale_y*this->prob_domain2D[i-this->Nxx].get_reflected(0) );
                }
            }
            //stub
            this->prob_domain2D[i].set_incidence(4,this->prob_domain2D[i].get_reflected(4) );
            i=i+1;
        }

}

void mesh_handler_2D::insert_structure_Naca_x_axis(int start_node,int end_node,naca_f func)
{
    double delta_x = 1/double(((end_node-start_node)));
    double x=0;

    cout<<"start_node"<<start_node<<endl;
    cout<<"end_node"<<end_node<<endl;

    for(int i = start_node; i <=end_node; i++)
    {

        //cout<<"node "<<x<<endl;wdw
        int f = func(x,this->dl);
        //cout<<" f: "<<f<<endl;
        x= x+delta_x;

        for(int j = 0 ; j<=f; j++)
        {
            this->make_PEC(i+j*Nxx);
            this->make_PEC(i-j*Nxx);
        }

    }

}

void mesh_handler_2D::insert_structure_dipole_y_axis(int mid_point,int length)
{
    cout<<" Inserting Dipole Antenna in.... "<< mid_point<<" "<<endl;
    cout<<" Length of dipole  " << 2*length*dl<<endl;
    for(int j = 1; j<=length; j++)
    {
        this->make_PEC(mid_point+j*Nxx);
        this->make_PEC(mid_point-j*Nxx);
    }

}

void:: mesh_handler_2D::insert_iris_2D(int thickness, float pos)
{
    int start_node = get_coordinate_iD_2D (width,height,npmlx,npmly,dl,dl,pos,0);
    int end_node = start_node+ Nx - 1;

    if( 2*thickness < Nx)
    {
        for( int i = 0; i<thickness; i++)
        {
            this->make_PEC(start_node+i);
            this->make_PEC(end_node-i);
        }
    }

    else cout <<" Cannot insert iris!... check geometry.....!! "<<endl;
}


void mesh_handler_2D::make_PEC(int node_id)
{
    this->prob_domain2D[node_id].make_PEC_node();
}

void mesh_handler_2D::print_nodes(int propty)
{
    int i(0);

    if( propty == 1)  cout<< "..... Displaying the iD of nodes in the computational domain.... "<<endl;
    if( propty == 2)  cout<< "..... Displaying the PML '1' nodes in the computational domain.... "<<endl;
    if( propty == 3)  cout<< "..... Displaying the PEC '1' nodes in the computational domain.... "<<endl;

    for( int y = 0; y< Nyy ; y++)
    {
        for ( int x = 0; x< Nxx; x++)
        {

            if( propty == 1) cout<<this->prob_domain2D[i].get_iD()<<"  ";
            if( propty == 2)  cout<<this->prob_domain2D[i].is_PML()<<"  ";
            if( propty == 3) cout<<this->prob_domain2D[i].is_PEC()<<"  ";
            if( propty == 4) cout<<" [ "<<this->prob_domain2D[i].get_cnductvty_x()<<"   :  "<<this->prob_domain2D[i].get_cnductvty_y()<<" ] ";
            i = i+1;
        }
        cout<<endl<<endl;
    }

}

void mesh_handler_2D:: print_all_PEC()
{
    int i(0);

    for( int y = 0; y< Nyy ; y++)
    {
        for ( int x = 0; x< Nxx; x++)
        {
            cout<<this->prob_domain2D[i].is_PEC()<<"  ";
            i = i+1;
        }
        cout<<endl<<endl;
    }


}

double mesh_handler_2D:: far_field_RCS(int bdry, double r )
{
    double d_surg = bdry*dl;
    if (bdry <= 0) return 0;

    int vertex1 = get_coordinate_iD_2D (width,height,npmlx,npmly,dl,dl+d_surg,dl+d_surg,0);  // bottom left
    int vertex2 = get_coordinate_iD_2D (width,height,npmlx,npmly,dl,width-d_surg,dl+d_surg,0);  // bottom right
    int vertex3 = get_coordinate_iD_2D (width,height,npmlx,npmly,dl,dl+d_surg,height-d_surg,0);  //top left
    int vertex4 = get_coordinate_iD_2D (width,height,npmlx,npmly,dl,width-d_surg,height-d_surg,0); // top right

    double ez1(0), hx1(0), hy1(0);
    double ez2(0), hx2(0), hy2(0);
    double Gz1(0), Gz2(0);
    double Hz1(0), Hz2(0);

    int j(0);

    double A= 0, B=0, C=0, D=0, E=0, F=0;

    // x directed surface piece A->B
    for( int i =vertex1; i<=vertex2; i++)
    {
        if( this->prob_domain2D[i].is_PML_shunt() )
        {
            cout<< " PML node!! " <<endl;
            cin>>i;
        }

        ez1 = this->prob_domain2D[i].get_Ez();
        hx1 = this->prob_domain2D[i].get_Hx();
        //hy1 = this->prob_domain2D[i].get_Hy();

        ez2 = this->prob_domain2D[vertex3+j].get_Ez();
        hx2 = this->prob_domain2D[vertex3+j].get_Hx();
        //hy2 = this->prob_domain2D[vertex3+j].get_Hy();

        //cout<<" id :  A -> B "<< i <<endl;
        //cout<<" id :  A -> B "<<vertex3+j<<endl;

        j = j+1;

        //surface piece 1
        A = A + hx1;  //
        B = B + ez1;  //

        //surface piece 2
        C  = C + (-1*hx2); //negative
        D  = D + (-1*ez2); //negative

    }

    double ez3(0), hx3(0), hy3(0);
    double ez4(0), hx4(0), hy4(0);
    j = 0;

    // y directed surface piece C->D
    for( int i =vertex1; i<=vertex3; i=i+Nxx)
    {
        // ez3 = this->prob_domain2D[i].get_Ez();
        // hx3 = this->prob_domain2D[i].get_Hx();
        hy3 = this->prob_domain2D[i].get_Hy();

        //ez4 = this->prob_domain2D[vertex3+j].get_Ez();
        //hx4 = this->prob_domain2D[vertex3+j].get_Hx();
        hy4 = this->prob_domain2D[vertex2+j*Nxx].get_Hy();

        //cout<<" id :  C -> D "<< i <<endl;
        //cout<<" id :  C -> D "<<vertex2+j*Nxx<<endl;

        j = j+1;

        //surface piece 3
        E = E + (-1*hy3);  //negative

        //surface piece 4
        F  = F + hy4; //
    }

    double coef1 =  u0/(4*pi*r);
    double coef2 =  1/(4*pi*r*c);
    double Ez_far_fld = coef1*(A+C+E+F) - (B+D)*coef2;

    return (Ez_far_fld*dl*1e-3);

}

mesh_handler_2D::~mesh_handler_2D()
{
    //dtor
}
/*
int get_coordinate_iD_2D(double width, double height,int npmlx,int npmly, double dl, double x, double y, double z)
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
