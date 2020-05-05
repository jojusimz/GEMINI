#ifndef  SIM_PARAM_H
#define  SIM_PARAM_H


struct sim_param
{


//Geometry
    double height, width,length ;

//Simulation parameters
    int total_steps;
    double dl, dt, tfactor;

//Excitation and output
    int excitation_type;

    double sin_freq, gauss_bw;

//medium parameters
    float er, ur ;

//boundary conditions
    double Zz1, Zz2, Zy1, Zy2, Zx1, Zx2;
    double reflct_z1, reflct_z2, reflct_x, reflct_y;

//Pml parameters
    int npmlx, npmly,conduct_prof ;
    double R_factor ;

//Excitation coordinates
    double x1, y1, z1;

//Observation coordinates
    double x2, y2, z2;

   int excitation_node_3D, output_node_3D, centre_node_3D;
   int excitation_node_2D, output_node_2D, centre_node_2D;

//Antenna parameters
    double lambda, antenna_length ;

//SCN PML node
    int PML_type;
    bool freespace;
    float sigma_e;

//Structures
    bool iris ;
    bool dipole;
    bool pec_cube ;
    bool naca0015 ;
    bool square_loop ;
    bool FSS ;
    bool FSS_JC ;
    bool dielectric_sheet ;

};

#endif //  _H
/*

 int exctn_nd = 0;
int obsrv_nd = 0;

exctn_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,x1,y1,dl);
obsrv_nd = return_coordinates_2D (width,height,npmlx,npmly,dl,x2,y2,dl);


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



*/




