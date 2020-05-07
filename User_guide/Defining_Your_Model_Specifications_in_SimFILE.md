## USING THE SIM.FILE 

To the best of my ability I have endeavoured to make the GEMINI simulation process as seamless as possible. To some degree I have been able to achieve this in a number of simple ways. One of which is the use of a simulaiton file **SIM.FILE**

By allowing the users to define their simulation parameters in a SIM.FILE I have removed the need for typing simulation data which in itself can be a inhibitive procedure especially if you are running a batch of simulations with lengthy run times. Also the **SIM.FILE** can be viewed as soft copy of the data of each simulation which can be referred to by the user at a later date.

### HOW IT WORKS:

The simulation data saved in **SIM.FILE** is parsed into a struct container called MODEL. The parser function is defined in gemini_utility.cpp. The struct container MODEL is a member variable of the wrapper class sim_handler. It is possible to modify MODEL, however all changes must be consistent with changes made in **SIM.FILE**. I am happy to guide you on this process as I imagine many of you would require a simulation process better tailored to the specific problems you are studying. 


### DEFINING YOUR SIMULATION DATA IN SIM.FILE

The whole **SIM.FILE** is a delimited text file that uses a comma to separate values. Each line of the file represents a simulation data. Each record consists of one or more fields, separated by commas. **The **SIM.FILE** can be edited in any text editing software. My preference editor is Microsoft Excel. When editing in Excel, make sure to save it as a .csv file.** 

An example template of **SIM.FILE** is shown below

![image](https://user-images.githubusercontent.com/60849864/81116681-b7e0ad80-8f1d-11ea-9a1b-de9d18fb9f56.png)

1.  The very first thing to do is to determine the type of simulation you wish to run i.e. 2D or 3D EM simulation. 

       **Modify Row 1, Col 1 in the simulation file.**
       
       * #w or #s: is used to specify 3D geometries      
       
       * #p or #t, is used to specify 2D geometries
       
2.  Determine the labels to your results. 
      
      **Modify Row 2, Col 3,4,5....... in the simulation file.**

3.  Define the model data e.g. geometry dimensions, excitation frequency, type of excitation etc. **See the example SIM.FILEs in the simulated examples depository.**
