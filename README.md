
**CONTENT**

- [GEMINI: A TIME DOMAIN TLM ELECTROMAGNETIC FIELD SIMULATOR DEVELOPED IN C++](#gemini-a-time-domain-tlm-electromagnetic-field-simulator-developed-in-c)
- [The TLM Method](#the-tlm-method)
- [How GEMINI Works](#how-gemini-works)
- [Features](#features)
- [C++ Source Code Features](#c-source-code-features)
- [Geometry/Material Features](#geometrymaterial-features)
- [Getting Setup](#getting-setup)
  * [Building GEMINI from source on Linux / Windows (MSYS)](#building-gemini-from-source-on-linux--windows-msys)
- [Publications Featuring GEMINI](#publications-featuring-gemini)
---------------------------------------------------------------------------------------------------------------------------------------

<!-- toc -->

# GEMINI: A TIME DOMAIN TLM ELECTROMAGNETIC FIELD SIMULATOR DEVELOPED IN C++
  ![image](https://user-images.githubusercontent.com/60849864/81127431-12860380-8f36-11ea-9360-1b0f7b29596e.png)

This repository hosts a C++ TLM electromagnetic field software suitable for simulating E and H field interactions in prescribed 2D and 3D cubic geometries discretized with a uniform mesh. 

At the start of my PhD there was no open source C++ TLM code available and so I developed GEMINI. Through it I have implemented novel formulations of PML algorithms and benchmarked their performance against existing published methods. Several numerical experiments have been carried out using GEMINI. I have now decided to release my source code as it may be useful to any researcher interested in simulating similar problems or wanting to contribute to the work. 

The implementation and approach taken to mesh and solve are simple and can be easily understood by researchers familiar with the TLM/Finite difference methods. The novelty however is the stable PML algorithms implemented the details of which are published in my journal papers  (See reference list [1]-[5]).

**Author**: Jomiloju Odeyemi

**Email**: Jomiloju.odeyemi@nottingham.ac.uk

**Weblink**: https://www.linkedin.com/in/jomiloju-odeyemi/

**Weblink**: https://www.researchgate.net/profile/Jomiloju_Odeyemi

# The TLM Method (At a glance)

The TLM method is a well established numerical method developed based on the discrete Huygens principle. It explores the analogy between wave propagation theory and transmission line theory to obtain field solutions to **hyperbolic partial differential equations**. Given the wave equation which is hyperbolic in its form:


![\frac{\partial ^2 E_z }{\partial x^2} + \frac{\partial^2 E_z}{ \partial y^2} = c^2 \frac{\partial ^2 E_z}{\partial t^2} ](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial%20%5E2%20E_z%20%7D%7B%5Cpartial%20x%5E2%7D%20%2B%20%5Cfrac%7B%5Cpartial%5E2%20E_z%7D%7B%20%5Cpartial%20y%5E2%7D%20%3D%20c%5E2%20%5Cfrac%7B%5Cpartial%20%5E2%20E_z%7D%7B%5Cpartial%20t%5E2%7D%20)


the TLM method obtains an equivalent circuit equation of a similar form by applying a finite difference discretization and circuit theory. The resulting TLM equation to (1) is given as


![\frac{\partial ^2 V_z }{L_x \partial x^2} + \frac{\partial^2 V_z}{ L_y \partial y^2} = C_z \frac{\partial ^2 V_z}{\partial t^2}](https://render.githubusercontent.com/render/math?math=%5Cfrac%7B%5Cpartial%20%5E2%20V_z%20%7D%7BL_x%20%5Cpartial%20x%5E2%7D%20%2B%20%5Cfrac%7B%5Cpartial%5E2%20V_z%7D%7B%20L_y%20%5Cpartial%20y%5E2%7D%20%3D%20C_z%20%5Cfrac%7B%5Cpartial%20%5E2%20V_z%7D%7B%5Cpartial%20t%5E2%7D)


According to circuit theory it is possible to obtain exact solutions to equation (2). In this manner the TLM method differs from other popular Finite difference based methods because it avoids a full discretization approach. Seeing that exact solutions to equation 2 it is therefore satisfactory to assume unconditional conditional stability is achieved.

The spatial domain is discretized into a matrix of TLM cells and an explicit time stepping algorithm ( referred to as the scatter and connect ) is applied at all points in space to obtain the field solutions.

# How GEMINI Works 

*At a glance*
![image](https://user-images.githubusercontent.com/60849864/80833434-6d41f700-8be6-11ea-800b-de9427535782.png)


# Features

   * Uniform rectangular mesh TLM solver for 2D (TLM shunt node) and 3D (SCN and HSCN) problems.
   
   * Boundary conditions:
   
        * Perfect Electrical Conductor
        
        * Magnetic walls
        
        * Matched Boundary Conditions
        
        * **Perfectly Matched Layer (Stretched coordinate PML implementation included)**
        
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
   
   * **Stable Perfectly Matched Layer implementation that can terminate arbitrary inhomogeneous media!!**
  
# C++ Source Code Features

   * Object Oriented C++ development
   
        * Modularity of code: Class libraries and their implementations can be used independently of the GEMINI package
        
        * Scalable: Addition of new features
        
   * Portable to Linux / Windows platforms
   
   * Parallel/Sequential mode (OpenMP)
   
   * Key features:
   
        * **<Sim_handler>** : a wrapper class which handles the whole simulation process
        
        * **<Mesh_handler>** : a wrapper class for handling a matrix of TLM nodes
        
        * **sim_file_parser()** : a function created to parse simulation data from an external file to a sim_handler class instance
        
            * Simulation/Model specification defined by user in a .csv file
            
            * Supports an automated running of simulation batches
            
            * Avoids hard typing of simulation data
            
        * Additon of new Geometry/Features easily achieved via class methods interface
        
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
  
        ``` make RUN_GEMINI ```
        
        
 *The terminal window*    
 
![image](https://user-images.githubusercontent.com/60849864/81107341-e99e4800-8f0e-11ea-81ab-bc9ee1486939.png)

 

---------------------------------------------------------------------------------------------------------------------------------------
  
# Publications Featuring GEMINI
1.	_J. Odeyemi, A. Vukovic, T. Benson, P. Sewell, “An improved PML implementation in the transmission line method,” IET 10th International Conference on Computational Electromagnetics, June 2019._ 

2.	_J. Odeyemi, A. Vukovic, T. Benson, P. Sewell, “Stretched-Coordinate PML in 2D TLM simulations,” IET Science, Measurement and Technology, vol.14, no.3, pp. 272-277, 2020._

3.  _J. Odeyemi, A. Vukovic, T. Benson, P. Sewell, “A complex domain mapping of the SCN for an effective PML implementation in TLM,” IEEE Open Journal of Antenna and Propagation, to be published DOI: 10.1109/OJAP.2020.2986293._

4.  _J. Odeyemi, C. Smartt, A. Vukovic, T. Benson, P. Sewell, “ PML Effectiveness in the Transmission Line Modelling Method for Radiation and Scattering Applications,”  European Conference on Antennas and Propagation, Copenhagen, 2020 (to be published)._

5.	_J. Odeyemi, M. Panitz, A. Vukovic, T. Benson, P. Sewell, “An Effective Stretched Coordinate TLM-PML Suitable for Analyzing Planar Periodic Structures,” IEEE Microwave and Wireless Components Letters (submitted and under review)._

