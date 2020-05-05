



/*
 * This file is part of GEMINI.
 *
 * Author: Jomiloju Odeyemi <simi.odeyemi@gmail.com>
 *
 */

#include <string>
#include <utility>      // std::pair
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <functional>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "SCN_node.h"
#include "HSCN_node.h"
#include "shunt_node.h"
#include "sim_handler.h"
#include "mesh_handler_2D.h"
#include "mesh_handler_HSCN.h"
#include "mesh_handler_SCN.h"
#include "gemini_utility.h"

using namespace std;

const double pi = 3.14159265;
const double  u0 = 12566370614e-16;
const double e0 = 88541878176e-22;
const double c = 299792458;
const double Z0 = sqrt(u0/e0);//

sim_handler::sim_handler( string simFILE, int no_sims)//constructor
{
    int simtype(-1);
    vector <string> E_filenames, H_filenames;

    v_sim_param = parse2DCsvFile(simtype, simFILE, E_filenames, H_filenames);
    Ey_file_names = E_filenames;
    Hz_file_names = H_filenames;
    sim_status = 0;
    sim_type = simtype;

}

sim_handler::~sim_handler()
{
    //dtor
}

void sim_handler::Model_definition_for_sim_batch(int i)
{

//Geonmetry
    model.height=v_sim_param[5][i], model.width=v_sim_param[6][i],model.length=v_sim_param[7][i];

//Simulation parameters
    model.total_steps = int(v_sim_param[8][i]);
    model.dl = v_sim_param[9][i],    model.tfactor = v_sim_param[10][i];
    model.excitation_type = int(v_sim_param[11][i]);
    model.sin_freq = v_sim_param[12][i], model.gauss_bw = v_sim_param[13][i];

//medium parameters
    model.er = v_sim_param[14][i], model.ur = v_sim_param[15][i];

//boundary conditions
    model.Zz1 = int(v_sim_param[16][i])*Z0, model.Zz2 = int(v_sim_param[17][i])*Z0, model.Zy1 = int(v_sim_param[18][i])*Z0, model.Zy2 = int(v_sim_param[19][i])*Z0, model.Zx1 = int(v_sim_param[20][i])*Z0, model.Zx2 = int(v_sim_param[21][i])*Z0;
    model.reflct_z1 = (model.Zz1 - Z0) /( model.Zz1 + Z0);
    model.reflct_z2 = (model.Zz2 - Z0) /( model.Zz2 + Z0);
    model.reflct_x = (model.Zx1 - Z0) /( model.Zx1 + Z0);    //condition for PMC
    model.reflct_y = (model.Zy1 - Z0) /( model.Zy1 + Z0);

    if( model.Zz1 < 0) model.reflct_z1 = 1;
    if( model.Zz2 < 0) model.reflct_z2 = 1;
    if( model.Zy1 < 0) model.reflct_y  = 1;
    if( model.Zx1 < 0) model.reflct_x  = 1;


//Pml parameters
    model.npmlx = int(v_sim_param[22][i]), model.npmly = int(v_sim_param[23][i]),model.conduct_prof = int(v_sim_param[25][i]);
    model.R_factor= v_sim_param[24][i] ;

//Excitation coordinates
    model.x1 = v_sim_param[26][i], model.y1 = v_sim_param[27][i], model.z1 = v_sim_param[28][i];

//Observation coordinates
    model.x2 = v_sim_param[29][i], model.y2 = v_sim_param[30][i], model.z2 = v_sim_param[31][i];

//Antenna parameters
    model.lambda = v_sim_param[33][i], model.antenna_length = v_sim_param[32][i];

//SCN PML node
    model.PML_type = v_sim_param[34][i];
    model.freespace = v_sim_param[35][i];
    model.sigma_e = v_sim_param[36][i];
    //Structures
    model.iris = v_sim_param[37][i];
    model.dipole = v_sim_param[38][i];
    model.pec_cube = v_sim_param[39][i];
    model.naca0015 = v_sim_param[37][i];
    model.square_loop = v_sim_param[40][i];
    model.FSS = v_sim_param[41][i];
    model.FSS_JC = v_sim_param[42][i];
    model.dielectric_sheet = v_sim_param[43][i];

//Hard coded changes - simulation constants
    model.tfactor = 0.9999999999998;

    model.dt =  model.tfactor*0.5*1e-3*model.dl/c;

//................................................................................................................

    unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);

    if( sim_type ==waveguide)
    {
        model.excitation_node_3D = return_coordinates_WG (model.width,model.height,model.length,model.npmlx,model.dl,model.x1,model.y1,model.z1);
        model.output_node_3D = return_coordinates_WG (model.width,model.height,model.length,model.npmlx,model.dl,model.x2,model.y2,model.z2);
        model.centre_node_3D = return_coordinates_WG (model.width,model.height,model.length,model.npmlx,model.dl,model.width/2,model.height/2,model.length/2);
    }

//Definitions for 3D cubic geometry
    else if ( sim_type ==cube_3D)
    {
         model.excitation_node_3D = return_coordinates_3D (model.width,model.height,model.length,model.npmlx,model.dl,model.x1,model.y1,model.z1);
         model.output_node_3D = return_coordinates_3D (model.width,model.height,model.length,model.npmlx,model.dl,model.x2,model.y2,model.z2);
         model.centre_node_3D = return_coordinates_3D (model.width,model.height,model.length,model.npmlx,model.dl,model.width/2,model.height/2,model.length/2);
    }

    if(( sim_type ==plane_2D) ||( sim_type ==waveguide_2D)  )
    {
        model.excitation_node_2D = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,model.x1,model.y1,model.dl);
        model.output_node_2D = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,model.x2,model.y2,model.dl);
        model.centre_node_2D = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,model.width/2,model.height/2,model.dl);

//Setting up Excitation
        int exct_type = model.excitation_type; // -1 = single node; 0 = line ;-2= te10

        if( exct_type == 0 )
        {
            model.excitation_node_2D = return_coordinates_2D (model.width, model.height, model.npmlx, model.npmly, model.dl, model.x1, model.dl, model.dl) ; // exciting from  left to right
            model.output_node_2D  = model.excitation_node_2D;    // line source

            cout<< " Excitation on line "<<model.excitation_node_2D<<endl;
            cout<< " output node observed in "<< model.output_node_2D <<endl;

        }

        else if( exct_type == -2)
        {
            model.excitation_node_2D = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,model.dl,model.y1,model.dl);
            model.output_node_2D = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,model.dl,model.y2,model.dl);

            cout<< " TE10  excitation in "<<model.excitation_node_2D<<endl;
            cout<< " output node observed in "<< model.output_node_2D <<endl;
        }

    }
}

void sim_handler::Display_simulation_details_3D(int i)
{

//Prints to screen
    cout<<" ................................................................"<<endl;
    cout<<" ......Display Simulation Parameters for Batch [ " << i<<" ]....."<<endl;
    cout<<" ................................................................"<<endl;
    cout<<" 3D dimensions width x height x length : " <<model.width << " x " << model.height <<" x "<<model.length <<endl;
    cout<<" Centre node : " << model.centre_node_3D<<endl;
    cout<<" Total steps : " << model.total_steps<<endl;
    cout<<" dl : "<<model.dl<<endl;
    cout<<" dt : "<< model.dt<<endl;
    cout<<" Excitation type : "<<model.excitation_type<<endl;
    cout<<" Sine frequency : " <<model.sin_freq<<endl;
    cout<<" Gaussian bandwidth : "<<model.gauss_bw<<endl;
    cout<<" Excitation node : " << model.excitation_node_3D<<endl;
    cout<<" Centre node : " << model.centre_node_3D<<endl;
    cout<<" Observation node : " << model.output_node_3D<<endl;
    cout<<" Reflection Factor : " <<model.R_factor <<endl;
    cout<<" Number of PML : "<<model.npmlx<<endl;
    cout<<" Conductivity profile : "<< model.conduct_prof<<endl;
    cout<<" Medium ? : ";
    if(model.freespace) cout<<" Freespace "<<endl;
    else cout<< " Dielectric space"<<endl;
    cout<<" relative permittivity : "<<model.ur<<endl;
    cout<<" relative permeability : "<<model.er<<endl;
    cout<<" electric conductivity : " << model.sigma_e <<" S/m "<<endl;
    cout<<" Impedances : "<< model.Zz1 <<" "<<model.Zz2<<" "<< model.Zy1 <<" "<<model.Zy2<< " "<<model.Zx1 <<" "<<model.Zx2<<endl;
    cout<<" PML type : "<< model.PML_type <<endl<<endl;
    cout<<" Boundary Reflection in x : " <<model.reflct_x <<endl;
    cout<<" Boundary Reflection in y : " <<model.reflct_y <<endl;
    cout<<" Boundary Reflection in z : " <<model.reflct_z2 <<endl;
    cout<<" ................................................................"<<endl;
    cout<<" ................................................................"<<endl;
    cout<<" ................................................................"<<endl;
//.......................................................................................................................................................................
//.............................................................................GEOMETRY DEFINTION AND MESHING............................................................
//.......................................................................................................................................................................

}


void sim_handler::Display_simulation_details_2D(int i)
{
    cout<<" ................................................................"<<endl;
    cout<<" ......Display Simulation Parameters for Batch [ " << i<<" ]....."<<endl;
    cout<<" ................................................................"<<endl;

// PRINT SIMULATION PARAMETERS TO SCREEN
    cout<<" 2D Geometry width x height x length : " <<model.width << " x " << model.height <<"  y " << endl;
    cout<<" Total steps " << model.total_steps<<endl;
    cout<<" dl : "<<model.dl<<endl;
    cout<<" dt : "<< model.dt ;
    cout<<" Excitation type: "<<model.excitation_type<<endl;
    cout<<" Sine frequency :" <<model.sin_freq<<endl;
    cout<<" Gaussian bandwidth : " <<model.gauss_bw<<endl;
    cout<<" Excitation node iD : " << model.excitation_node_2D<<endl;
    cout<<" Observation node iD:" << model.excitation_node_2D<<endl;
    cout<<" Centre node iD:" << model.centre_node_2D<<endl;
    cout<<" Reflection Factor :" <<model.R_factor <<endl;
    cout<<" Number of PML in x plane : "<<model.npmlx<<endl;
    cout<<" Number of PML in y plane: "<<model.npmly<<endl;
    cout<<" Conductivity profile :"<< model.conduct_prof<<endl;
    cout<<" Space ? :  ";
    if(model.freespace) cout<<" Freespace "<<endl;
    else cout<< " Dielectric space"<<endl;
    cout<<" ur :"<<model.ur<<endl;
    cout<<" er :"<<model.er<<endl;
    cout<<" electric conductivity " << model.sigma_e <<" S/m "<<endl;
    cout<<" Impedances "<< model.Zy1 <<" "<<model.Zy2<< " "<<model.Zx1 <<" "<<model.Zx2<<endl;
    cout<<" SCN PML node "<< model.PML_type <<endl;
    cout<<" Antenna size "<<model.antenna_length<<endl;
    cout<<"Boundary Reflection in x " <<model.reflct_x <<endl;
    cout<<"Boundary Reflection in y " <<model.reflct_y <<endl;
    cout<<" ................................................................"<<endl;
    cout<<" ................................................................"<<endl;
    cout<<" ................................................................"<<endl;

}

 void sim_handler::Output_results_to_file_3D(int i, string directory_string, mesh_handler_SCN * SCN_mesh)
 {

     unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);

     cout<<" Write output to file "<<endl;

     ostringstream sii;
     sii << i;
     string s_i = sii.str();

     if(sim_type == cube_3D)  // 3D cubic structure
        {
            SCN_mesh->write_output_file( "./"+directory_string+"/"+s_i+"_Ex_3D.txt",   "./"+directory_string+"/"+s_i+"_"+Ey_file_names[i]+".txt",       "./"+directory_string+"/"+s_i+"_Ez_3D.txt",
                                            "./"+directory_string+"/"+s_i+"_Hx_3D.txt",   "./"+directory_string+"/"+s_i+"_Hy_3D.txt",  "./"+directory_string+"/"+s_i+"_"+Hz_file_names[i]+".txt"         );   // note that the folder_string is added to the file name
        }

        else // Waveguide
        {
            SCN_mesh->write_output_file("./"+directory_string+"/"+s_i+"_Ex_WG.txt",   "./"+directory_string+"/"+s_i+"_"+Ey_file_names[i]+".txt",       "./"+directory_string+"/"+s_i+"_Ez_WG.txt",
                                            "./"+directory_string+"/"+s_i+"_"+Hz_file_names[i]+".txt",   "./"+directory_string+"/"+s_i+"_Hy_WG.txt", "./"+directory_string+"/"+s_i+"_Hz_WG.txt") ;    // note that the folder_string is added to the file name
        }

//UNCOMMENT the output options
            //SCN_mesh->write_output_file_for_zx_plane ("./"+directory_string+"/"+s_i+"_EY_zx_plane.csv") ;
            // SCN_mesh->write_output_file_for_xy_plane ("./"+directory_string+"/"+"EY_xy_plane.csv") ;
            // SCN_mesh->write_incident_file_for_xy_plane ("./"+directory_string+"/"+"V_incident_wg_xy.csv") ;
            //SCN_mesh->write_output_file_for_line( "./"+directory_string+"/"+s_i+"_Ey_along_line.csv" ) ;

 }


 void sim_handler::Insert_structures_in_3D_domain( mesh_handler_SCN * SCN_mesh)
 {
    unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);

    int Nz = (model.length/model.dl +0.5);
//Number of nodes in entire computation along each coordinate axis
    int Nx_ = int ( (model.width  / model.dl) + 0.5) +1;       // note the added node which indicates position where the SD is placed.
    int Ny_ = int ( (model.height / model.dl) + 0.5) +1;
    int Ny = Ny_-1;
    int Nz_ = int ( (model.length / model.dl) + 0.5) +1;

//Number of nodes in entire computation along each coordinate axis including PML
    int Nxx_ = Nx_ + 2*model.npmlx;
    int Nyy_ = Ny_ + 2*model.npmlx;
    int Nzz_ = Nz_ + 2*model.npmlx;

    int npml = model.npmlx;
//inserting into the computational domain of cubic geometry
        if( (sim_type == cube_3D)&&(Nz>1) )
        {
            model.excitation_node_3D = SCN_mesh->get_true_centre_node();  // setting the excitation at the centre node
            model.output_node_3D = model.excitation_node_3D;
            //excitation_node_3D = return_coordinates_3D (width,height,length,npml,dl,x1,y1,z1);
            // output_node_3D = return_coordinates_3D (width,height,length,npml,dl,x2,y2,z2);

            //SCN_mesh->print_WG_plane_par(int(length/dl),2);
            if (model.dipole)
            {
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
                model.excitation_type = -1;
                SCN_mesh->place_dipole_antenna(model.antenna_length);
                model.excitation_node_3D = SCN_mesh->get_true_centre_node();  // setting the excitation at the centre node
                model.output_node_3D = model.excitation_node_3D;                // setting the output at the centre node
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
            }

            else if( model.square_loop)
            {
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
                model.excitation_type = -1;
                SCN_mesh->insert_square_loop(model.antenna_length);
                model.excitation_node_3D = SCN_mesh->get_true_centre_node();  // setting the excitation at the centre node
                model.output_node_3D = model.excitation_node_3D;                // setting the output at the centre node
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
            }

            else if (model.pec_cube)
            {
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
                SCN_mesh->insert_perfect_cube(5,10,10,15);
                //SCN_mesh->print_WG_plane_par(npml+24,3);

                int nz1=12;  // 12 nodes in
                // exciting at plane
                if ((model.excitation_type ==-2)||( model.excitation_type >=0))
                {
                    model.excitation_node_3D = Nxx_*Nyy_*(npml+nz1) + Nxx_*npml + npml;
                    model.output_node_3D =  Nxx_*Nyy_*(npml+35) + Nxx_*(npml+10) + npml+10;
                }

                //excite domain at single node
                if( model.excitation_type ==-1)
                {
                    model.excitation_node_3D = return_coordinates_3D (model.width, model.height, model.length, npml, model.dl, model.x1, model.y1, model.z1);
                    model.output_node_3D = return_coordinates_3D (model.width, model.height, model.length, npml, model.dl, model.x2, model.y2, model.z2);
                }
            }
        }

// inserting into the waveguide geometry
        if( (sim_type ==waveguide)&&(Ny>1) )
        {

            if( model.iris )
            {
                float z_plane = 11 - 10*model.dl; // position of iris relative to the far end of waveguide
                int node_id = z_plane/model.dl + 0.5; // node id
                int Ny = model.height/(3*model.dl) + 0.5;   // height of the iris

                SCN_mesh->insert_iris(Ny,z_plane);
            }


            else if( model.FSS )
            {
                SCN_mesh->insert_FSS_square(2.5,model.length/2);
                model.excitation_node_3D = return_coordinates_WG (model.width, model.height, model.length, npml, model.dl,0,0, model.length/2 - 1*model.dl); // fixed point of excitation always 10 cells from FSS
                model.output_node_3D = return_coordinates_WG (model.width, model.height, model.length, npml, model.dl,0,0,model.length/2 +1*model.dl);   //fixed point of observation always 10 cells from FSS - other side
            }

            else if( model.FSS_JC )
            {
                SCN_mesh->insert_FSS_jerusalem_cross(5, model.length/2);
                model.excitation_node_3D = return_coordinates_WG (model.width, model.height, model.length,npml, model.dl,0,0, model.length/2 - 1*model.dl); // fixed point of excitation always 10 cells from FSS
                model.output_node_3D = return_coordinates_WG (model.width, model.height, model.length,npml, model.dl,0,0,model.length/2 +1*model.dl);   //fixed point of observation always 10 cells from FSS - other side
            }

        }

 }

 void sim_handler::Insert_structures_in_2D_domain( mesh_handler_2D * shunt_node_mesh)
 {

    if(model.naca0015)  //inserting the NACA0015
    {
        float  start_x = abs(model.width-1000)/2;     // position of the naca0015 along the x axis
        //Constructing the geometry for the naca0015
        int start_node = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,start_x,model.height/2,0);      // coordinate iD for naca0015
        int end_node = return_coordinates_2D (model.width,model.height,model.npmlx,model.npmly,model.dl,start_x+1000,model.height/2,0);   // coordinate iD for naca0015
        shunt_node_mesh->insert_structure_Naca_x_axis(start_node,end_node,naca_f(0.15));

        // Can change this to observe at the opposite end.
        //obsrv_nd = end_node + d
        }

 }

void sim_handler::Begin_TLM_simulation_SCN(int i,string directory_string)
{
    unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);

//..........................................................................................................................
// SIMULATION PARAMETERS ARE DEFINED BELOW. IF EDITED, THIS MUST BE DONE IN LINE WITH MODIFICATIONS IN THE SIM.FILE SYSTEM
//............................................................................................................................

    Model_definition_for_sim_batch(i);

// Other Model Definitions

    model.dt =  model.tfactor*0.5*1e-3*model.dl/c;

    int npml = model.npmlx;

    int Nz = (model.length/model.dl +0.5);

//............................................................................................................................
//............................................PRINT TO SCREEN..............................................................
//............................................................................................................................

    Display_simulation_details_3D( i);

//.......................................................................................................................................................................
//.............................................................................GEOMETRY DEFINTION AND MESHING............................................................
//.......................................................................................................................................................................

//3D waveguide geometry definition

    mesh_handler_SCN *SCN_mesh;

    if ( sim_type == waveguide )    SCN_mesh = new mesh_handler_SCN(sim_type, model.width, model.height, model.length, model.dl,
                                                                    model.tfactor, model.er, model.sigma_e, model.PML_type, npml, model.conduct_prof, model.R_factor);    // waveguide


    else if ( sim_type == cube_3D )
    {
//2D SCN waveguide geometry definition
        if( Nz==1 ) SCN_mesh = new mesh_handler_SCN(model.width, model.height, sim_type, model.dl, model.tfactor, model.er, model.sigma_e,
                                                    model.PML_type, npml, model.conduct_prof, model.R_factor);  //

//Cubic domain geometry definition
        else  SCN_mesh = new mesh_handler_SCN(model.width, model.height, model.length, sim_type, model.dl, model.tfactor,
                                                    model.er, model.sigma_e, model.PML_type, npml, model.conduct_prof, model.R_factor);
    }

//............................................................................................................................
//...............................INSERTING STRUCTURES INTO THE COMPUTATIONAL DOMAIN......................................................
//............................................................................................................................

    this->Insert_structures_in_3D_domain( SCN_mesh );

//............................................................................................................................................
//........................................................SET UP FOR TIME STEPPING ALGORITHM..................................................
//............................................................................................................................................


//Setting neighbours
    cout<<endl<<"Setting neighbours"<<endl;
    SCN_mesh->set_WG_neighbour();

//Set special nodes
    cout<<endl<<"Setting special nodes "<<endl<<endl;
    SCN_mesh->set_special_nodes();

//Setting the connect bins
    SCN_mesh->create_connect_bins();
    SCN_mesh->print_connect_bins();

//........................................................................................................................
//........................................................TLM TIME STEPPING ALGORITHM..................................................
//........................................................................................................................

    cout<< ".................THE TLM ALGORITHM: EXCITE->SCATTER->CONNECT->OUTPUT............"<<endl;

//2D SCN simulation
    if(Nz==1)        cout<<" THIS PORTION OF CODE REQUIRES YOUR ATTENTION AND HAS BEEN COMMENTED OUT " <<endl;

//3D SCN
    else   //3D SCN
    {

//EXCITE->SCATTER->CONNECT->OUTPUT
        SCN_mesh->TLM_simulation1( model.total_steps, model.tfactor, model.output_node_3D, model.excitation_node_3D, model.Zz1, model.Zz2, model.Zy1, model.Zy2, model.Zx1,
                                   model.Zx2, f_gaussian(model.sin_freq,model.gauss_bw,model.dt), model.excitation_type,0,true ); //dl = 0.5
//Print output to screen
        cout<< ".......Print to EY field values to screen......"<<endl;
        SCN_mesh->print_Ey_output();

// writing field values to file for output
        Output_results_to_file_3D(i, directory_string, SCN_mesh);
    }

//...........................................................................................................................
//......................................................END OF SIMULATION.....................................................................
//...........................................................................................................................
    delete SCN_mesh;

    cout<< ".........................................................................................................................."<<endl;
    cout<< ".....................................................END OF SIMULATION OF " << i << "..................................................."<<endl;
    cout<< ".........................................................................................................................." <<endl;

}


void sim_handler::Begin_TLM_simulation_2D(int i,string directory_string)
{
   unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);

//..........................................................................................................................
// SIMULATION PARAMETERS ARE DEFINED BELOW. IF EDITED, THIS MUST BE DONE IN LINE WITH MODIFICATIONS IN THE SIM.FILE SYSTEM
//............................................................................................................................

    this->Model_definition_for_sim_batch(i);

// Additional Model Definitions

    model.dt = 1e-3*model.dl/(sqrt(2)*c);

    int Nx = model.width/model.dl+0.5;

    //excitation distance for naca0015 is set by variable d
    int d=10; // distance from Naca0015

    if (( int (model.width/model.dl +0.5) - int(1000/model.dl + 0.5) ) <= 20 )  // setting a fixed point of excitation according to the dimension of the naca0015
    {
        d = 0.5*( int (model.width/model.dl +0.5) - int(1000/model.dl + 0.5) )-1;     // changes the variable d because the distance from the airfoil to the boundary is less than 10nodes away.
    }

    //boundary conditions for domain
    vector <int> bdry_cdn;
    bdry_cdn.push_back(model.reflct_x);
    bdry_cdn.push_back(model.reflct_y);

//............................................................................................................................
//............................................PRINT TO SCREEN..............................................................
//............................................................................................................................

    this->Display_simulation_details_2D( i);

//.......................................................................................................................................................................
//.............................................................................GEOMETRY DEFINTION AND MESHING............................................................
//.......................................................................................................................................................................

//3D waveguide geometry definition

    mesh_handler_2D *shunt_node_mesh;

    shunt_node_mesh = new mesh_handler_2D(model.freespace, model.er, model.sigma_e, model.height, model.width,
                                          model.dl, model.npmlx, model.npmly, model.R_factor, model.conduct_prof);

//............................................................................................................................
//...............................INSERTING STRUCTURES INTO THE COMPUTATIONAL DOMAIN......................................................
//............................................................................................................................

    this->Insert_structures_in_2D_domain( shunt_node_mesh );

//........................................................................................................................
//........................................................DISPLAY COMPUTAITONAL DOMAIN DETAILS TO SCREEN.................................................
//........................................................................................................................

    shunt_node_mesh->print_nodes(1);
    cout<<endl<<endl;
    shunt_node_mesh->print_all_PEC();

//........................................................................................................................
//........................................................TLM TIME STEPPING ALGORITHM..................................................
//........................................................................................................................


    cout<< ".................THE TLM ALGORITHM: EXCITE->SCATTER->CONNECT->OUTPUT............"<<endl;


//EXCITE->SCATTER->CONNECT->OUTPUT

     ostringstream sii;
     sii << i;
     string s_i = sii.str();

     shunt_node_mesh->Simulate_2D(model.total_steps, bdry_cdn, model.excitation_node_2D, 5800, model.output_node_2D,
                                     f_gaussian(model.sin_freq, model.gauss_bw, model.dt), model.excitation_type, d,0);  //void Simulate_2D(int ttltimestep,vector<int> &bdry_cdn,int excit_iD,int excit_lngth,int obsrv_iD);

     shunt_node_mesh->write_output_file_2D("./"+directory_string+"/"+s_i+"_"+Ey_file_names[i]+".txt","./"+directory_string+"/"+s_i+"_Hx_.txt","./"+directory_string+"/"+"Is2.txt","./"+directory_string+"/"+"Vs2.txt");
     shunt_node_mesh->write_output_file_for_xy_plane_shuntnode("./"+directory_string+"/"+s_i+"2D__plane_animation_test_for_stub.csv");

//...........................................................................................................................
//......................................................END OF SIMULATION.....................................................................
//...........................................................................................................................
    delete shunt_node_mesh;

    cout<<" End of Simulation " << i <<endl<<endl<<endl;
}



void sim_handler::Begin_TLM_simulation(int i,string directory_string)
{
    unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);



    if (( sim_type  == waveguide) || ( sim_type == cube_3D))
    {
        this->Begin_TLM_simulation_SCN(i, directory_string);
    }


    else if (( sim_type  == waveguide_2D) || ( sim_type == plane_2D))
    {
cout<<" here "<<sim_type<<endl;
        this->Begin_TLM_simulation_2D(i, directory_string);
    }

}

/*

    int debug;
    int simtype = sim_type;
    ostringstream sii;
    sii << i;
    string s_i = sii.str();



//..........................................................................................................................
// SIMULATION PARAMETERS ARE DEFINED BELOW. IF EDITED, THIS MUST BE DONE IN LINE WITH MODIFICATIONS IN THE SIM.FILE SYSTEM
//............................................................................................................................

    Model_definition_for_sim_batch(i);

/*
//Geonmetry
    double height=sim_param[5][i], width=sim_param[6][i],length=sim_param[7][i];

//Simulation parameters
    int total_steps = int(sim_param[8][i]);
    double dl = sim_param[9][i], tfactor = sim_param[10][i];
    int excitation_type = int(sim_param[11][i]);
    double sin_freq = sim_param[12][i], gauss_bw = sim_param[13][i];

//medium parameters
    float er = sim_param[14][i], ur = sim_param[15][i];
    double Zz1 = int(sim_param[16][i])*Z0, Zz2 = int(sim_param[17][i])*Z0, Zy1 = int(sim_param[18][i])*Z0, Zy2 = int(sim_param[19][i])*Z0, Zx1 = int(sim_param[20][i])*Z0, Zx2 = int(sim_param[21][i])*Z0;

//Pml parameters
    int npmlx = int(sim_param[22][i]), npmly = int(sim_param[23][i]),conduct_prof = int(sim_param[25][i]);
    double R_factor= sim_param[24][i] ;

//Excitation coordinates
    double x1 = sim_param[26][i], y1 = sim_param[27][i], z1 = sim_param[28][i];

//Observation coordinates
    double x2 = sim_param[29][i], y2 = sim_param[30][i], z2 =sim_param[31][i];

//Antenna parameters
    double lambda = sim_param[33][i], antenna_length = sim_param[32][i];

//SCN PML node
    int PML_type = sim_param[34][i];
    bool freespace = sim_param[35][i];
    float sigma_e = sim_param[36][i];
    //Structures
    bool iris = sim_param[37][i];
    bool dipole = sim_param[38][i];
    bool pec_cube = sim_param[39][i];
    bool naca0015 = sim_param[37][i];
    bool square_loop = sim_param[40][i];
    bool FSS = sim_param[41][i];
    bool FSS_JC = sim_param[42][i];
    bool dielectric_sheet = sim_param[43][i];

//Hard coded changes - simulation constants
    tfactor = 0.9999999999998;



    double dt =  tfactor*0.5*1e-3*dl/c;

// Node Definitions
    int Nx_ = int ( (width  / dl) + 0.5) +1;       // note the added node which indicates position where the SD is placed.
    int Ny_ = int ( (height / dl) + 0.5) +1;
    int Ny = Ny_-1;
    int Nz_ = int ( (length / dl) + 0.5) +1;

//Number of nodes in entire computation along each coordinate axis including PML
    int Nxx_ = Nx_ + 2*npmlx;
    int Nyy_ = Ny_ + 2*npmlx;
    int Nzz_ = Nz_ + 2*npmlx;

    int excitation_node_3D = 0, output_node_3D = 0, centre_node = 0;
    unsigned int waveguide(0),cube_3D(1),waveguide_2D(2),plane_2D(3),waveguide_HSCN(4),cube_3D_HSCN(5);

//............................................................................................................................
//............................................PRINT TO SCREEN..............................................................
//............................................................................................................................

/*
    cout<<" ...............Display Simulation Parameters for Batch [ " << i<<" ]..........."<<endl;
    cout<<" .............................................................................."<<endl;
    cout<<" .............................................................................."<<endl;


    if( simtype == waveguide)     cout<<" ..................SIMULATING WAVEGUIDE using SCN...................................."<<endl<<endl;
    else if (simtype == cube_3D)      cout<<" .................SIMULATING 3D CUBIC DOMAIN...................................."<<endl<<endl;
    else if (simtype == waveguide_2D)  cout<<" ..................2D PROBLEM WAVEGUIDE................................."<<endl<<endl;
    else if (simtype == plane_2D)  cout<<" ..................2D DIPOLE PROBLEM.................................."<<endl<<endl;
    else if (simtype == waveguide_HSCN)  cout<<" ..................SIMULATING WAVEGUIDE using HSCN...................................."<<endl<<endl;
    else if (simtype == cube_3D_HSCN)  cout<<" ..................SIMULATING 3D CUBIC DOMAIN using HSCN...................................."<<endl<<endl;

    int Nz = (length/dl +0.5);                                            // Number of cells in the z direction
    if( Nz == 1)  cout<<" ..................2D SCN.................................."<<endl<<endl;


// GEOMETRY

// 3D Waveguide
    if(( simtype ==waveguide) || (simtype ==cube_3D))
    {
        int npml = npmlx;

//Definitions for waveguide geometry
        if( simtype ==waveguide)
        {
            // cout<<" here "<<endl;
            excitation_node_3D = return_coordinates_WG (width,height,length,npml,dl,x1,y1,z1);
            output_node_3D = return_coordinates_WG (width,height,length,npml,dl,x2,y2,z2);
            centre_node = return_coordinates_WG (width,height,length,npml,dl,width/2,height/2,length/2);
        }
//Definitions for 3D cubic geometry
        if ( simtype ==cube_3D)
        {

            excitation_node_3D = return_coordinates_3D (width,height,length,npml,dl,x1,y1,z1);
            output_node_3D = return_coordinates_3D (width,height,length,npml,dl,x2,y2,z2);
            centre_node = return_coordinates_3D (width,height,length,npml,dl,width/2,height/2,length/2);
        }

/*
//Prints to screen
        cout<<" 3D dimensions width x height x length : " <<width << " x " << height <<" x "<<length <<endl;
        cout<<" Centre node " << centre_node<<endl;
        cout<<" Total steps " << total_steps<<endl;
        cout<<" dl : "<<dl<<endl;
        cout<<" dt : "<< dt<<endl;
        cout<<" Excitation type: "<<excitation_type<<endl;
        cout<<" Sine frequency :" <<sin_freq<<endl;
        cout<<" Gaussian bandwidth : "<<gauss_bw<<endl;
        cout<<" Excitation node : " << excitation_node_3D<<endl;
        cout<<" Observation node :" << output_node_3D<<endl;
        cout<<" Reflection Factor :" <<R_factor <<endl;
        cout<<" Number of PML : "<<npml<<endl;
        cout<<" Conductivity profile :"<< conduct_prof<<endl;
        cout<<" Medium ? :  ";
        if(freespace) cout<<" Freespace "<<endl;
        else cout<< " Dielectric space"<<endl;
        cout<<" relative permittivity :"<<ur<<endl;
        cout<<" relative permeability :"<<er<<endl;
        cout<<" electric conductivity " << sigma_e <<" S/m "<<endl;
        cout<<" Impedances "<< Zz1 <<" "<<Zz2<<" "<< Zy1 <<" "<<Zy2<< " "<<Zx1 <<" "<<Zx2<<endl;
        cout<<" PML type "<< PML_type <<endl<<endl;

// Boundary conditions of computational domain
        double reflct_z1 = (Zz1 - Z0) /( Zz1 + Z0);
        double reflct_z2 = (Zz2 - Z0) /( Zz2 + Z0);
        double reflct_x = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
        double reflct_y = (Zy1 - Z0) /( Zy1 + Z0);

        if( Zz1 < 0) reflct_z1 = 1;
        if( Zz2 < 0) reflct_z2 = 1;
        if( Zy1 < 0) reflct_y  = 1;
        if( Zx1 < 0) reflct_x  = 1;

        cout<<"Reflection in x " <<reflct_x <<endl;
        cout<<"Reflection in y " <<reflct_y <<endl;
        cout<<"Reflection in z " <<reflct_z2 <<endl;
        cout<<" ..............................................................................."<<endl;
        cout<<" ..............................................................................."<<endl<<endl;
        cout<<" ..................................SIMULATION............................................."<<endl;


//.......................................................................................................................................................................
//.............................................................................GEOMETRY DEFINTION AND MESHING............................................................
//.......................................................................................................................................................................

//3D waveguide geometry definition

        mesh_handler_SCN *SCN_mesh;

        if ( simtype == waveguide )
        {
            SCN_mesh = new mesh_handler_SCN(simtype,width,height,length,dl,tfactor,er,sigma_e, PML_type,npml,conduct_prof,R_factor);    // waveguide
        }

        else if ( simtype == cube_3D )
        {
//2D SCN waveguide geometry definition
            if( Nz==1 )
            {
                SCN_mesh = new mesh_handler_SCN(width,height,simtype,dl,tfactor,er,sigma_e,PML_type, npml,conduct_prof,R_factor);  //
            }

//Cubic domain geometry definition
            else
            {
                SCN_mesh = new mesh_handler_SCN(width,height,length,simtype,dl,tfactor,er,sigma_e,PML_type, npml,conduct_prof,R_factor);
            }

        //excitation_node_3D = SCN_mesh->get_true_centre_node();
        //cout<<centre_node <<endl;
        //SCN_mesh->get_true_centre_node();

        }



//............................................................................................................................
//...............................INSERTING STRUCTURES INTO THE COMPUTATIONAL DOMAIN......................................................
//............................................................................................................................

//inserting into the computational domain of cubic geometry
        if( (simtype == cube_3D)&&(Nz>1) )
        {
            excitation_node_3D = SCN_mesh->get_true_centre_node();  // setting the excitation at the centre node
            output_node_3D = excitation_node_3D;
            //excitation_node_3D = return_coordinates_3D (width,height,length,npml,dl,x1,y1,z1);
            // output_node_3D = return_coordinates_3D (width,height,length,npml,dl,x2,y2,z2);

            //SCN_mesh->print_WG_plane_par(int(length/dl),2);
            if (dipole)
            {
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
                excitation_type = -1;
                SCN_mesh->place_dipole_antenna(antenna_length);
                excitation_node_3D = SCN_mesh->get_true_centre_node();  // setting the excitation at the centre node
                output_node_3D = excitation_node_3D;                // setting the output at the centre node
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
            }

            else if( square_loop)
            {
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
                excitation_type = -1;
                SCN_mesh->insert_square_loop(antenna_length);
                excitation_node_3D = SCN_mesh->get_true_centre_node();  // setting the excitation at the centre node
                output_node_3D = excitation_node_3D;                // setting the output at the centre node
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
            }

            else if (pec_cube)
            {
                //SCN_mesh->print_WG_plane_par(int(length/dl),2);
                SCN_mesh->insert_perfect_cube(5,10,10,15);
                //SCN_mesh->print_WG_plane_par(npml+24,3);

                int nz1=12;  // 12 nodes in
                // exciting at plane
                if ((excitation_type ==-2)||( excitation_type >=0))
                {
                    excitation_node_3D = Nxx_*Nyy_*(npml+nz1) + Nxx_*npml + npml;
                    output_node_3D =  Nxx_*Nyy_*(npml+35) + Nxx_*(npml+10) + npml+10;
                }

                //excite domain at single node
                if( excitation_type ==-1)
                {
                    excitation_node_3D = return_coordinates_3D (width,height,length,npml,dl,x1,y1,z1);
                    output_node_3D = return_coordinates_3D (width,height,length,npml,dl,x2,y2,z2);
                }
            }
        }

// inserting into the waveguide geometry
        if( (simtype ==waveguide)&&(Ny>1) )
        {

            if( iris )
            {
                float z_plane = 11 - 10*dl; // position of iris relative to the far end of waveguide
                int node_id = z_plane/dl + 0.5; // node id
                int Ny = height/(3*dl) + 0.5;   // height of the iris

                SCN_mesh->insert_iris(Ny,z_plane);
            }


            else if( FSS )
            {
                SCN_mesh->insert_FSS_square(2.5,length/2);
                excitation_node_3D = return_coordinates_WG (width,height,length,npml,dl,0,0,length/2 - 1*dl); // fixed point of excitation always 10 cells from FSS
                output_node_3D = return_coordinates_WG (width,height,length,npml,dl,0,0,length/2 +1*dl);   //fixed point of observation always 10 cells from FSS - other side
            }

            else if( FSS_JC )
            {
                SCN_mesh->insert_FSS_jerusalem_cross(5,length/2);
                excitation_node_3D = return_coordinates_WG (width,height,length,npml,dl,0,0,length/2 - 1*dl); // fixed point of excitation always 10 cells from FSS
                output_node_3D = return_coordinates_WG (width,height,length,npml,dl,0,0,length/2 +1*dl);   //fixed point of observation always 10 cells from FSS - other side
            }

        }

//............................................................................................................................................
//........................................................SET UP FOR TIME STEPPING ALGORITHM..................................................
//............................................................................................................................................


//Setting neighbours
        cout<<endl<<"Setting neighbours"<<endl;
        SCN_mesh->set_WG_neighbour();

//Set special nodes
        cout<<endl<<"Setting special nodes "<<endl<<endl;
        SCN_mesh->set_special_nodes();

//Setting the connect bins
        SCN_mesh->create_connect_bins();
        SCN_mesh->print_connect_bins();


//........................................................................................................................
//........................................................TLM TIME STEPPING ALGORITHM..................................................
//........................................................................................................................

        cout<< ".................THE TLM ALGORITHM: EXCITE->SCATTER->CONNECT->OUTPUT............"<<endl;


//2D SCN simulation
        if(Nz==1)  //2D SCN only appropriate in single node excitations!!
        {
            cout<<" THIS PORTION OF CODE REQUIRES YOUR ATTENTION AND HAS BEEN COMMENTED OUT " <<endl;
            cout<<" THIS PORTION OF CODE REQUIRES YOUR ATTENTION AND HAS BEEN COMMENTED OUT " <<endl;
            cout<<" THIS PORTION OF CODE REQUIRES YOUR ATTENTION AND HAS BEEN COMMENTED OUT " <<endl;
            cout<<" THIS PORTION OF CODE REQUIRES YOUR ATTENTION AND HAS BEEN COMMENTED OUT " <<endl;
        }

//3D SCN
        else   //3D SCN
        {

//EXCITE->SCATTER->CONNECT->OUTPUT
            SCN_mesh->TLM_simulation1( total_steps, tfactor, output_node_3D, excitation_node_3D, Zz1, Zz2, Zy1, Zy2, Zx1,
                                          Zx2, f_gaussian(sin_freq,gauss_bw,dt), excitation_type,0,true ); //dl = 0.5

            cout<< ".......Print to EY field values to screen......"<<endl;
            SCN_mesh->print_Ey_output();

// writing field values to file for output
            cout<<" Write output to file "<<endl;
            if(simtype == cube_3D)  // 3D cubic structure
            {
                SCN_mesh->write_output_file( "./"+directory_string+"/"+s_i+"_Ex_3D.txt",   "./"+directory_string+"/"+s_i+"_"+Ey_file_names[i]+".txt",       "./"+directory_string+"/"+s_i+"_Ez_3D.txt",
                                                "./"+directory_string+"/"+s_i+"_Hx_3D.txt",   "./"+directory_string+"/"+s_i+"_Hy_3D.txt",  "./"+directory_string+"/"+s_i+"_"+Hz_file_names[i]+".txt"         );   // note that the folder_string is added to the file name
            }

            else // Waveguide
            {
                SCN_mesh->write_output_file("./"+directory_string+"/"+s_i+"_Ex_WG.txt",   "./"+directory_string+"/"+s_i+"_"+Ey_file_names[i]+".txt",       "./"+directory_string+"/"+s_i+"_Ez_WG.txt",
                                               "./"+directory_string+"/"+s_i+"_"+Hz_file_names[i]+".txt",   "./"+directory_string+"/"+s_i+"_Hy_WG.txt", "./"+directory_string+"/"+s_i+"_Hz_WG.txt") ;    // note that the folder_string is added to the file name
            }

//UNCOMMENT the output options
            //SCN_mesh->write_output_file_for_zx_plane ("./"+directory_string+"/"+s_i+"_EY_zx_plane.csv") ;
            // SCN_mesh->write_output_file_for_xy_plane ("./"+directory_string+"/"+"EY_xy_plane.csv") ;
            // SCN_mesh->write_incident_file_for_xy_plane ("./"+directory_string+"/"+"V_incident_wg_xy.csv") ;
            //SCN_mesh->write_output_file_for_line( "./"+directory_string+"/"+s_i+"_Ey_along_line.csv" ) ;
        }

//...........................................................................................................................
//...........................................................................................................................
//...........................................................................................................................

delete SCN_mesh;

//...........................................................................................................................
//......................................................END OF SIMULATION.....................................................................
//...........................................................................................................................


    cout<<" End of Simulation " << i <<endl<<endl<<endl;
    cout<<" ............................................Warnings.................................... "<<endl;
    cout<<" .........................................................................................."<<endl<<endl;
    if ( Ny==1) cout<< " Simulating the waveguide in 2D  mode" <<endl<< " Only suited for empty waveguide "<<endl;
    if ( Nz==1) cout<< " Simulating the cubic domain in 2D mode." <<endl<<" Only suited for Single point excitation " <<endl;
    if ( dipole) cout << " Inserted Dipole" <<endl;
    if ( pec_cube) cout << " Inserted PEC CUBE " <<endl;
    if ( iris ) cout<< " Inserted Iris in waveguide" <<endl;
    if ( square_loop) cout<< " Inserting Square Loop"<<endl;

}

//2D shunt node
    else if( (simtype==2) || (simtype == 3))
    {
        int exctn_nd = 0;
        int obsrv_nd = 0;

        exctn_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,x1,y1,dl);
        obsrv_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,x2,y2,dl);

        //Construct Geometry handler
        mesh_handler_2D *Automate_2D;
        Automate_2D = new mesh_handler_2D(freespace,er,sigma_e,height,width,dl,npmlx,npmly,R_factor,conduct_prof);


        if(dipole)  // to ensure symmetry when placing the dipole in the computational domain
        {
            if( int(width / dl +0.5)%2 == 0 ) width = width + dl;
            if( int(height / dl + 0.5)%2 == 0 )height = height + dl;
            int mid_node = return_coordinates_2D (width,height,npmlx,npmly,dl,width/2,height/2,0);
            exctn_nd = mid_node;
            obsrv_nd = mid_node;
            int ant_len_no =  0.25*antenna_length /dl + 0.5-5 ;
            Automate_2D->insert_structure_dipole_y_axis(mid_node,ant_len_no);
            excitation_type = -1;

        }

        dt = 1e-3*dl/(sqrt(2)*c);

//Variables for the Naca0015
        float  start_x = abs(width-1000)/2;     //width distance to naca0015
        //Constructing the geometry for the naca0015
        int start_node = return_coordinates_2D (width,height,npmlx,npmly,dl,start_x,height/2,0);      // for naca0015
        int end_node = return_coordinates_2D (width,height,npmlx,npmly,dl,start_x+1000,height/2,0);   //for naca0015

//FIXES THE POINT OF EXCITATION IN THE NACA0015 - variable d denotes the distance from the object.
        int d = 10;
        if(naca0015)
        {
            float  start_x = abs(width-1000)/2;     //width distance to naca0015
            //Constructing the geometry for the naca0015
            int start_node = return_coordinates_2D (width,height,npmlx,npmly,dl,start_x,height/2,0);      // for naca0015
            int end_node = return_coordinates_2D (width,height,npmlx,npmly,dl,start_x+1000,height/2,0);   //for naca0015
            Automate_2D->insert_structure_Naca_x_axis(start_node,end_node,naca_f(0.15));

            d=10; // distance from Naca0015

            if (( int (width/dl +0.5) - int(1000/dl + 0.5) ) <= 20 )  // setting a fixed point of excitation according to the dimension of the naca0015
            {
                cout<< " changing excitation point to....";
                d = 0.5*( int (width/dl +0.5) - int(1000/dl + 0.5) )-1;     // changes the variable d because the distance from the airfoil to the boundary is less than 10nodes away.
            }

            if(excitation_type != -1 ) exctn_nd = return_coordinates_2D (width,height,npmlx,npmly,x1,dl,dl,0); //start_node - d;
            obsrv_nd = exctn_nd;
            // Can change this to observe at the opposite end.
            //obsrv_nd = end_node + d

        }

// PRINT SIMULATION PARAMETERS TO SCREEN
        cout<<" 2D Geometry width x height x length : " <<width << " x " << height <<"  y " << endl;
        cout<<" Total steps " << total_steps<<endl;
        cout<<" dl : "<<dl<<endl;
        cout<<" dt : "<< dt ;
        cout<<" Excitation type: "<<excitation_type<<endl;
        cout<<" Sine frequency :" <<sin_freq<<endl;
        cout<<" Gaussian bandwidth : "<<gauss_bw<<endl;
        cout<<" Excitation node : " << exctn_nd<<endl;
        cout<<" Observation node :" << obsrv_nd<<endl;
        cout<<" Reflection Factor :" <<R_factor <<endl;
        cout<<" Number of PML in x : "<<npmlx<<endl;
        cout<<" Number of PML in y : "<<npmly<<endl;
        cout<<" Conductivity profile :"<< conduct_prof<<endl;
        cout<<" Space ? :  ";
        if(freespace) cout<<" Freespace "<<endl;
        else cout<< " Dielectric space"<<endl;
        cout<<" ur :"<<ur<<endl;
        cout<<" er :"<<er<<endl;
        cout<<" electric conductivity " << sigma_e <<" S/m "<<endl;
        cout<<" Impedances "<< Zy1 <<" "<<Zy2<< " "<<Zx1 <<" "<<Zx2<<endl;
        cout<<" SCN PML node "<< PML_type <<endl;
        cout<<" Antenna size "<<antenna_length<<endl<<endl;

// Boundary conditions of computational domain
        int reflct_x = (Zx1 - Z0) /( Zx1 + Z0);    //condition for PMC
        int reflct_y = (Zy1 - Z0) /( Zy1 + Z0);

        if( Zy1 < 0) reflct_y  = 1;             // Magnetic wall Z
        if( Zx1 < 0) reflct_x  = 1;

        vector <int> bdry_cdn;
        bdry_cdn.push_back(reflct_x);
        bdry_cdn.push_back(reflct_y);

        int Nx = width/dl+0.5;

        cout<<"Boundary Reflection in x " <<reflct_x <<endl;
        cout<<"Boundary Reflection in y " <<reflct_y <<endl;
        cout<<" ..............................................................................."<<endl;
        cout<<" ..............................................................................."<<endl<<endl;
        cout<<" ..................................SIMULATION............................................."<<endl;


//Setting up Excitation
        int exct_type = excitation_type; // -1 = single node; 0 = line ;-2= te10

        if( exct_type == 0 )
        {
            if( naca0015 )
            {
                start_x = start_x - 5*dl;
                exctn_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,x1,dl,dl) ; // exciting from  left to right
            }
            else exctn_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,2*dl,dl,dl) ;
            obsrv_nd = exctn_nd;    // line source
        }

        if( exct_type == -2)
        {
            exctn_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,dl,y1,dl);
            obsrv_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,dl,y2,dl);
            cout<< " TE10  excitation in "<<exctn_nd<<endl;
            cout<< " output node observed in "<<obsrv_nd <<endl;
        }
//Print
        // Automate_2D->print_nodes(1);
        cout<<endl<<endl;
        // Automate_2D->print_all_PEC();
//Simulate
        Automate_2D->Simulate_2D(total_steps, bdry_cdn, exctn_nd, 5800, obsrv_nd,f_gaussian(sin_freq,gauss_bw,dt),exct_type,d,0);  //void Simulate_2D(int ttltimestep,vector<int> &bdry_cdn,int excit_iD,int excit_lngth,int obsrv_iD);
        Automate_2D->write_output_file_2D("./"+directory_string+"/"+s_i+"_"+Ey_file_names[i]+".txt","./"+directory_string+"/"+s_i+"_Hx_.txt","./"+directory_string+"/"+"Is2.txt","./"+directory_string+"/"+"Vs2.txt");
        Automate_2D->write_output_file_for_xy_plane_shuntnode("./"+directory_string+"/"+s_i+"2D__plane_animation_test_for_stub.csv");

        delete Automate_2D;
    }

//HSCN
    else if(( simtype ==5) || (simtype ==4))
    {

    }

}


*/


/*sim_handler::sim_handler(const sim_handler& other)
{
    //copy ctor
}

sim_handler& sim_handler::operator=(const sim_handler& rhs)
{
    if (this == &rhs) return *this; // handle self assignment
    //assignment operator
    return *this;
}
*/
