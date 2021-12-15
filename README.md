# LBcuda

CUDA Fortran code to solve Lattice Boltzmann equations of a multi-component flow interacting with rigid spherical particles
                                                                        
Fabio Bonaccorso         IIT-CLNS,Rome                     Italy        
Marco Lauricella         IAC-CNR, Rome                     Italy        
Andrea Montessori        IAC-CNR, Rome                     Italy        
Adriano Tiribocchi       IAC-CNR, Rome                     Italy        
                                                                                                                
with contributions from:                                                
                                                                        
Giorgio Amati         CINECA-CED, Rome                     Italy        
Massimo Bernaschi        IAC-CNR, Rome                     Italy        
Sauro Succi              IAC-CNR, Rome                     Italy        
                                                                        
The software development process has received funding from the          
European Research Council under the Horizon 2020 Programme              
Grant Agreement n. 739964 (COPMAT). Marco Lauricella acknowledges       
the project “3D-Phys" (PRIN 2017PHRM8X) for support from MIUR.          
                                                                        
This is an experimental code. The authors accept no responsibility      
for the performance of the code or for the correctness of the results.  
                                                                        
The code is licensed under the 3-Clause BSD License (BSD-3-Clause).     
                                                                        
Structure                                                               
                                                                        
LBcuda is supplied as a main UNIX directory with 3 subdirectories.      
All the source code files are contained in the main directory.          
The 'Tests' sub-directory contains different test cases that can        
help the user to edit new input files. Further the 'Tests' sub-directory
contains the input files to run the test cases reported in Table 2 and  
3 of the main article 'LBcuda: a high-performance CUDA port of LBsoft   
for simulation of colloidal systems'.                                   
                                                                        
Compiling LBcuda                                                        
                                                                        
The main directory stores a UNIX makefile that assembles the            
executable versions of the code both in GPU and multi-GPUs version. 
Note the makefile could be eventually modified for several common        
workstations and NVIDIA devices into the main directory,                
where the code is compiled and linked. Finally, the     
binary executable file can be run in the main directory or copy in,     
an external location which is intended to be the working directory from 
which jobs are submitted for execution and the data files manipulated.  
                                                                        
On Windows system we advice the user to compile LBcuda under the        
Windows native Ubuntu shell. Note that the NVIDIA HPC Software          
Development  Kit (SDK) is necessary to compile and run LBcuda on        
the GPU and multi-GPUs parallel mode.                                                        
                                                                        
Executing LBcuda                                                        
                                                                        
To run LBcuda, it is necessary first to ensure that the program is      
compiled (from the source sub-directory) and that the file input.dat    
is present in the main directory.                                       
All output data files will be returned to the main directory.           
Remember that the input file HAS TO BE NAMED 'input.dat' and put in     
the main directory. Few example input files are contained in            
the 'Tests' sub-directory which can be used as test cases.              
                                                                        
Example command to run LBcuda in serial mode:                           
                                                                        
mpirun -np 1 ./main.x                                                                
                                                                        
Example command to run LBcuda in parallel mode on 8 GPUs:               
                                                                        
mpirun -np 8 ./main_mpi.x                                               
                                                                        
