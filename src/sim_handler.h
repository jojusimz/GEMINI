#ifndef sim_handler_H
#define sim_handler_H

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
#include <limits>   // double precision7*89


#include "SCN_node.h"
#include "mesh_handler_SCN.h"
#include "mesh_handler_2D.h"
#include "sim_param.h"

using namespace std;

class sim_handler
{
    public:
        sim_handler(string simFILE, int no_sims=1);
        virtual ~sim_handler();
        void Model_definition_for_sim_batch(int i);
        void Begin_TLM_simulation(int i, string folder_str);
        void Begin_TLM_simulation_SCN(int i, string folder_str);
        void Begin_TLM_simulation_2D(int i, string folder_str);
        void Insert_structures_in_3D_domain( mesh_handler_SCN * SCN_mesh);
        void Insert_structures_in_2D_domain( mesh_handler_2D * shunt_node_mesh);
        void Output_results_to_file_3D(int i, string _results_directory, mesh_handler_SCN * SCN_mesh);
        void Output_results_to_file_2D(mesh_handler_2D *shunt_node_mesh);
        void Display_simulation_details_3D(int i);
        void Display_simulation_details_2D(int i);
        //sim_handler(const sim_handler& other);
        //sim_handler& operator=(const sim_handler& other);
    private:
        vector <vector< double> > v_sim_param;
        int sim_status;
        int sim_type;
        vector < string > Ey_file_names;
        vector < string > Hz_file_names;
        sim_param model;
};

#endif // sim_handler_H
