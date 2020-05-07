
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

#include "gemini_utility.h"
//#include<windows.h>       //Uncomment for windows
//#include <ctime>          //uncomment for windows


using namespace std;

int simulation_start_up(string & _sim_file, string & _results_directory)
{
    cout<<endl<<endl;
    cout<< " ---------------------------------------------------------------"<<endl;
    cout<< " G E M I N I : C++ TLM ELECTROMAGNETIC FIELD SIMULATOR "<<endl;
    cout<< " ---------------------------------------------------------------"<<endl;
    cout<< " Created by: J.Odeyemi (2017)  "<<endl;
    cout<< " ---------------------------------------------------------------"<<endl<<endl;

    cout<< " Enter '1' for the default simulation setting " <<endl;
    cout<< " Enter '0' for the user defined simulation setting "<<endl<<endl;

    unsigned int usr_input(1);
    unsigned int num_sim_to_run (1);

    _sim_file = "SIMfile.csv";                    //default simulation file
    _results_directory = "Results";                    //default name of results directory

    cout<<"....";
    cin>> usr_input;

    switch (usr_input)
    {
        case 1:
            return (num_sim_to_run);
        case 0:
            cout<< "Enter the number of simulations you wish to run....."<<endl;
            cin>> num_sim_to_run;
            return (num_sim_to_run);
        default:
            return (1);
    }

}

void set_sim_file_name(string & _sim_file_name)
{
    cout<< " Enter the name of the .csv simulation file" <<endl;
    cin>>_sim_file_name;
}

void set_results_folder_name (string & _results_folder_name)
{
    // create folder

    cout<< " Enter the name of the user defined results folder" <<endl;
    cin>>_results_folder_name;

    string foldername(""+_results_folder_name);
    //CreateDirectory(foldername.c_str(),NULL); // windows specific utility

}



vector<vector<double> > sim_file_parser(int &simtype, string inputFileName, vector< string > &e_filenames, vector< string > &h_filenames)
{

    vector< vector <double> > data;

    ifstream inputFile(inputFileName.c_str());              //input file streaming object
    int l = 0;
    int countt = 0;
    double temp =0;
    string str, str1("#E"),str2("#H"), str3("#s"), str4("#w"),str5("#q"),str6("#t"),str7("#l"),str8("#k");

    while (inputFile)                                //returns -1 for EOF or EOL "\n"
    {
        l++;
        string s, f;
        if (!getline(inputFile, s)) break;          //reads line from file stream into string s - checks if valid, breaks otherwise

        str = s.substr(0,2);
        //strll = "#H";

        //cout<<str<<endl;
        if(str == str3) simtype = 1;          // 3D Cubic domain
        if(str == str4) simtype = 0;          // waveguide
        if(str == str5) simtype = 2;          // 2D Plane
        if(str == str6) simtype = 3;          // 2D Plane
        if(str == str7) simtype = 4;          // Waveguide HSCN
        if(str == str8) simtype = 5;          //  3D Cubic domain HSCN


        if ( str == str1 )
        {
            istringstream ss(s);

            while (ss)
            {
                string line;
                if (!getline(ss, line, ',') ) break;
                if (countt >=2)                       //interested in the nth column.
                {
                    e_filenames.push_back(line);
                    //cout<<line<<endl;
                }
                countt++;                           // have to read every element of string
            }
            countt = 0;
        } //"#H"

        if  (str == str2)
        {
            istringstream ss(s);

            while (ss)
            {
                string line;
                if (!getline(ss, line, ',') ) break;
                if (countt >=2)                       //interested in the nth column.
                {
                    h_filenames.push_back(line);
                    //cout<<line<<endl;
                }
                countt++;                           // have to read every element of string
            }
            countt = 0;
        }

        if (s[0] != '#')
        {
            istringstream ss(s);
            vector<double> v_temp;

            while (ss)                             // splits the string up and converts each component to a double.
            {
                string line;
                if (!getline(ss, line, ',') ) break;
                if (countt >=2)                       //interested in the nth column.
                {
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

    if (!inputFile.eof())
    {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return data;
}



vector<vector<double> > parseIncident_2File(string inputFileName)
{

    vector< vector <double> > data;

    ifstream inputFile(inputFileName.c_str());              // input file streaming object
    int l = 0;
    int countt = 0;
    double temp =0;

    while (inputFile)
    {
        vector< double > v_inc;
        l++ ;                          // lines count
        string line,s,f;
        if (!getline(inputFile, line, '\n'))break; // READS FIRST ROW INTO LINE

        int sz = 5, siz = 0, i= 0, k =0;                             // maximum length in first row is 2+2+1 = 5

        if (l>1) sz = line.length();

        for(int ii = 1; ii<=sz; ii = ii+0)
        {
            double d=0;
            string words ("");
            // separate into words
            do
            {
                words = words + line[i];
                i = i+1;
                if(i==sz) break;//cout<<" here "<<words<<endl;
            }
            while   (line[i] != ',');

            if(words == ",") break;  // takes care of '',''

            siz = words.length()+1 ;  //
            ii = ii + siz;
            i=ii-1;
            istringstream ss(words);
            ss>>d;


            v_inc.push_back(d);

            //        cout<<d<<endl;
            //        cout<<v_inc[k]; //int k =0;
            // if (l==3) { cout<<" CONVERTED INTEGER " <<d<< " " << v_inc[k]<<endl; cout<<ii<<" "<<sz<<endl; }
            //k = k+1;
        }

        //cout<<" row = "<< l <<" ...size of vector = "<<v_inc.size()<<endl;
        data.push_back(v_inc);
        // line.erase ( remove(line.begin(), line.end(), ','), line.end()  );
    }

    //cout<<"size "<<data[18][4]<<endl;
    //cin>>temp;
    //cout<<"size "<<data[3][data[2].size()-1]<<endl;

    return data;
}


int get_coordinate_iD_WG(double width, double height, double length,int npml, double dl,
                          double x, double y, double z) //only waveguide
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

int get_coordinate_iD_3D(double width, double height, double length,int n_PMl,
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

int get_coordinate_iD_2D(double width, double height,int npmlx,int npmly,
                          double dl, double x, double y, double z)
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

