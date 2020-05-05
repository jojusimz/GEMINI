/*
 * This file is part of GEMINI software package
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
#include <iterator>
#include <functional>
#include <cmath>
#include <iomanip>
#include <limits>

//#include<windows.h>       //Uncomment for windows
//#include <ctime>          //uncomment for windows

#include "sim_handler.h"
#include "gemini_utility.h"

using namespace std;

const double pi = 3.14159265;
const double  u0 = 12566370614e-16;
const double e0 = 88541878176e-22;
const double c = 299792458;
const double Z0 = sqrt(u0/e0);
const double er(1),ur(1);

int main(int argc, char* argv[])
{

//Begin TLM SIMULATIONS
    string sim_file;
    string results_directory;

    int num_of_sim = simulation_start_up(sim_file, results_directory);

    cout<<endl<< " Simulation file:    " << sim_file<<endl;
    cout<< " Results directory: " << results_directory<<endl<<endl;

    cout<<" ..............................................................."<<endl<<endl;
    cout<< " Running ' " << num_of_sim << " ' TLM simulations from  "<< sim_file<< " "<<endl<<endl;

    sim_handler TLM_sims(sim_file,num_of_sim);

// Runs the Batch of TLM simulations
    for( int i=0; i<num_of_sim; ++i)
    {
        TLM_sims.Begin_TLM_simulation(i,results_directory);
    }

    return 0;

}
