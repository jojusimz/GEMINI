## USING THE SIM.FILE 

To the best of my ability I have endeavoured to make the GEMINI simulation process as seamless as possible. To some degree I have been able to achieve this in a number of simple ways. One of which is the use of a simulaiton file **SIM.FILE**

Allowing the users to define their simulation parameters in a SIM.FILE removes the need for hard coding simulation data which can be inhibitive procedure especially if you are running a batch of simulations with lengthy run times. Also the **SIM.FILE** can be viewed as soft copy of the data of each simulation which can be referred to by the user at a later date.

### HOW IT WORKS:

The simulation data saved in **SIM.FILE** is parsed into a struct container called MODEL. The parser function is defined in gemini_utility.cpp. The struct container MODEL is a member variable of the wrapper class sim_handler. It is possible to modify MODEL, however all changes must be consistent with changes made in **SIM.FILE**. I am happy to guide you on this process as I imagine many of you would require a simulation process better tailored to the specific problems you are studying. 


### DEFINING YOUR SIMULATION DATA IN SIM.FILE

The **SIM.FILE** template is given below. 

![image](https://user-images.githubusercontent.com/60849864/81116681-b7e0ad80-8f1d-11ea-9a1b-de9d18fb9f56.png)

1.  The very first thing to do is to determine the type of simulation you wish to run i.e. 2D or 3D structures. 
       
2.  Determine the labels to your results. 

3.  Define the model data e.g. geometry dimensions, excitation frequency, type of excitation etc.
