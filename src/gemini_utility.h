#ifndef GEMINI_UTILITY_H
#define GEMINI_UTILITY_H

#include <string>
#include <vector>
using namespace std;

int simulation_start_up(string & _sim_file, string & _results_directory);
void set_sim_file(string & _sim_file);
void set_results_directory (string & _results_directory);
vector<vector<double> > parse2DCsvFile(int &simtype, string inputFileName, vector< string >& e_filemames, vector< string >& h_filenames);
vector<vector<double> > parseIncident_2File(string inputFileName);

int return_coordinates_WG(double width, double height, double length,int npml, double dl, double x, double y, double zz);
int return_coordinates_3D(double width, double height, double length,int npml, double dl, double x, double y, double zz);
int return_coordinates_2D(double width, double height,int npmlx,int npmly, double dl, double x, double y, double z);
//vector<vector<double> > parse2DCsvFile(int &simtype, string inputFileName, vector< string >& e_filemames, vector< string >& h_filenames);

class f_gaussian
{
public:
    f_gaussian(double f0_,double bw, double dt_) : f0(f0_),bandwidth(bw),dt(dt_){  }
    inline double operator () (int ti) const { return ( (dt*ti-(dt))*bandwidth* exp( -pow( bandwidth*(dt*ti-(dt*10) ),2))*sin(2*3.14159265*f0*dt*ti)  ); }

private:
    double f0;
    double bandwidth;
    double dt;
};
#endif //
