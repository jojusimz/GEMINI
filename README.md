- [GEMINI: A TIME DOMAIN TLM ELECTROMAGNETIC FIELD SOLVER DEVELOPED IN C++](#gemini-a-time-domain-tlm-electromagnetic-field-solver-developed-in-c)
- [The TLM Method](#the-tlm-method)
- [How GEMINI Works](#how-gemini-works)
- [Features](#features)
- [C++ Source Code Features](#c-source-code-features)
- [Geometry/Material Features](#geometrymaterial-features)
- [Getting Setup](#getting-setup)
  * [Building GEMINI from source on Linux / Windows (MSYS)](#building-gemini-from-source-on-linux--windows-msys)
- [Publications Featuring GEMINI](#publications-featuring-gemini)

<!-- toc -->

# GEMINI: A TIME DOMAIN TLM ELECTROMAGNETIC FIELD SOLVER DEVELOPED IN C++

This repository hosts a C++ TLM electromagnetic field solver suitable for simulating E and H field interactions in prescribed 2D and 3D cubic geometries discretized with a uniform mesh. 

This software was developed to facilitate my PhD research. Through it I implemented novel formulations of PML algorithms and benchmarked these against existing published methods. I have now decided to release the source code as it may be useful to any researcher interested in simulating similar geometries or wanting to build on the existing work. 
The implementation and approach taken to mesh and solve are simple and can be easily understood my researchers familiar with the TLM/Finite difference methods. The novelty however is the stable PML algorithms implemented the details of which are published in my journal papers  (See reference list [1]-[5]).

Author: Jomiloju Odeyemi

**Contact me** : 

email: simi.odeyemi@gmail.com

https://www.linkedin.com/in/jomiloju-odeyemi/

https://www.researchgate.net/profile/Jomiloju_Odeyemi

# The TLM Method

# How GEMINI Works


![image](https://user-images.githubusercontent.com/60849864/80833434-6d41f700-8be6-11ea-800b-de9427535782.png)


# Features

   * Uniform rectangular mesh TLM solver for 2D (TLM shunt node) and 3D (SCN and HSCN) problems.
   
   * Boundary conditions:
   
        * Perfect Electrical Conductor
        
        * Magnetic walls
        
        * Matched Boundary Conditions
        
        * **Perfectly Matched Layer (Split field and Stretched coordinate implementations included)**
        
   * Material Space:
   
        * isotropic / anisotropic permeability and permittivity material space
        
        * lossy material space (electric and magnetic losses)
        
   * Outputs E and H fields with point , line and plane output options.
   
        * Output file format: .txt/ .csv files
        
   * Sequential/Parallel (OpenMP) implementation options
   
   * Source code is compilable on Linux and Windows using MSY. See setUP instructions section.
   
   * Excitation options:
   
        * point source
        
        * line source (2D)
        
        * Plane source : TE or TEM sources 
        
        * Gaussian derivative pulse
        
   * Provides option to run batch of simulations in an automated manner
   
   * **Stable Perfectly Matched Layer implementation for all parameter choices!!**
  
# C++ Source Code Features

   * Object Oriented C++ development
   
        * Modularity of code: Class libraries and their implementations can be used independently of the GEMINI package
        
        * Scalable: Addition of new features
        
   * Portable to Linux / Windows platforms
   
   * Parallel/Sequential mode
   
   * Key features:
   
        * **<Sim_handler>** : a wrapper class which handles the whole simulation process
        
        * **<Mesh_handler>** : a wrapper class for handling a matrix of TLM nodes
        
        * **data_parser()** : a function created to parse simulation data from an external file to a sim_handler class instance
        
            * Simulation/Model specification defined by user in a .csv file.
            
            * Supports an automated running of a batch of simulation.
            
            * Avoids hardcoding of simulation data.
            
        * Additon of new Geometry/Features easily achieved via class methods interface.
        
# Geometry/Material Features

The following material/geometrical features are the insertable in the computational domain:

  * Waveguide capacitive iris
  
  * Dipole Antenna
  
  * PEC cube 
  
  * Frequency Selective Surface : Square Aperture FSS and Jerusalem Cross FSS
  
  * Dielectric substrate.
  
The insertion of new features is straightforward and can be achieved by a combination of the class_methods:

   *  Mesh_handler.get_coordinate_iD (x,y,z)
   
   *  Mesh_handler.Set_material (material_parameter_identifier)
   
For more info on this feature please contact me as this may not be obvious to the novice user.

# Getting Setup
 
  ## Building GEMINI from source on Linux / Windows (MSYS)
  
  Compilation Instructions
  
  1. **CHECKLIST**
  
     * Msys terminal (Mandatory for Windows. Download via MingGW Installation Manager)
     
     * g++ compiler (Mandatory)
     
     * MATLAB or Python to run Analytic scripts on Results (Mandatory)
  
  2.  **START** BASH terminal (Linux) / MSYS shell (Windows)
    
  3.  **NAVIGATE** to project folder and **CLONE** GEMINI repository to
  
        ```cd folder/to/clone-into/```
      
        ```git clone https://github.com/jojusimz/GEMINI.git```
    
  4.  **NAVIGATE** to GEMINI source code directory and **CREATE** a results directory 
  
         ```cd src```
       
         ```mkdir GEMINI_results```
  
  5.  **BUILD** the source code to executable software 
    
        ``` make ```
      
  6.  **RUN** the executable
  
        ``` .\simulate_GEMINI ```
  
# Publications Featuring GEMINI
1.	_J. Odeyemi, A. Vukovic, T. Benson, P. Sewell, “An improved PML implementation in the transmission line method,” IET 10th International Conference on Computational Electromagnetics, June 2019._ 

2.	_J. Odeyemi, A. Vukovic, T. Benson, P. Sewell, “Stretched-Coordinate PML in 2D TLM simulations,” IET Science, Measurement and Technology, vol.14, no.3, pp. 272-277, 2020._

3.  _J. Odeyemi, A. Vukovic, T. Benson, P. Sewell, “A complex domain mapping of the SCN for an effective PML implementation in TLM,” IEEE Open Journal of Antenna and Propagation, to be published DOI: 10.1109/OJAP.2020.2986293._

4.  _J. Odeyemi, C. Smartt, A. Vukovic, T. Benson, P. Sewell, “ PML Effectiveness in the Transmission Line Modelling Method for Radiation and Scattering Applications,”  European Conference on Antennas and Propagation, Copenhagen, 2020 (to be published)._

5.	_J. Odeyemi, M. Panitz, A. Vukovic, T. Benson, P. Sewell, “An Effective Stretched Coordinate TLM-PML Suitable for Analyzing Planar Periodic Structures,” IEEE Microwave and Wireless Components Letters (submitted and under review)._

