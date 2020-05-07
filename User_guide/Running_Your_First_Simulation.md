### Step by step instructions for running your first GEMINI simulation

* Successfully compile GEMINI into an .exe

* Define your Model parameters by editing a copy of **SIM.FILE**. (Edit SIMfile.csv )

* Move the **SIM.FILE** into the same directory as the compiled code.

* Create a folder called *Results* in the same directory as the compiled file.

* Run the following command in the BASH shell

       make RUN_GEMINI
      
    * Your Terminal should display the following information. Choose option '1' or '0' 
    
![image](https://user-images.githubusercontent.com/60849864/81107341-e99e4800-8f0e-11ea-81ab-bc9ee1486939.png)  


* GEMINI will now:

    * Create the Geometry
       
    * Discretize the Geometry into TLM cells
       
     * Iterate TLM algorithm for the total number of time steps set out in the **SIM.FILE** 

* Upon completion the *Results_directory* will be populated with results from the simulation

* **WELL DONE!! YOU HAVE NOW SIMULATED YOUR FIRST ELECTROMAGNETIC PROBLEM USING GEMINI YOU CAN PROCEED TO RUNNING YOUR ANALYTICS SCRIPTS ON THE RESULT FILES**


