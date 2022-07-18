# Nonlinear isogeometric analysis of thin-walled structures 
# - nligaStruct

------------------------------------------------------------
This is the README file of the open-source library nligaStruct: 
1. An open isogeometric framework for linear and nonlinear 
analysis of thin-walled structures is developed. 
2. The Beizer extraction of NURBS and T-spline is embedded. 
3. Both Reissner-Mindlin and Kirchhoff-Love theories are 
studied for the linear analysis of plates and shells. 
4. Large deformation of Kirchhoff-Love shell considering 
compressible and incompressible hyperelastic materials is 
implemented and validated. 
5. The complete source codes and data are fully provided 
and free to access.

------------------------------------------------------------
First of all, add the "nligaStruct" folder to the path. 
Right-click on the "nligaStruct" folder and select "Add to 
Path" -> "Selected Folders and Subfolders". 

To have a quick look at nligaStruc, directly switch to the 
subfolder "outputs" and run the script 'demo'.


------------------------------------------------------------
**NOTE: If you want to build geometrical models with NURBS 
basis functions, please add NLIGA library or nurbs toolbox 
to your path. 
NLIGA: https://sourceforge.net/projects/nliga/
nurbs: https://octave.sourceforge.io/nurbs/index.html


------------------------------------------------------------
The detailed contents of the folders in NLIGA are listed below:
1. examples
This folder includes all numerical examples named with the 
prefix ex.
(1) Two benchmarks for static bending of RM and KL plates;
(2) One benchmark for free vibration of RM and KL plates;
(3) Three benchmarks for small deformation of RM and KL 
shells, widely known as shell obstacle course problems;
(4) One benchmark for free vibration of RM and KL shells;
(5) Five benchmarks for large deformation of KL shells 
considering both geometrical and material nonlinearities.
(6) Data for linear and nonlinear analysis of shell examples.

------------------------------------------------------------
2. functions
The necessary functions used for numerical examples are 
collected in this subfolder, e.g., evaluation of the basis 
functions of T-spline and Bezier functions, and their 
derivatives.

------------------------------------------------------------
3. igafiles
Include all *.iga files extracted from T-spline models 
for simulation. Some functions for NURBS-based modeling 
are also provided.

------------------------------------------------------------
4. outputs
The subfolder includes the converged displacement results 
saved as *.msh files at all load steps for large deformation 
of thin shells using KL assumption.


------------------------------------------------------------
*** If you have any comments or suggestions, please feel free 
to contact us! 
- Xiaoxiao Du, 
- Beihang University,
- duxxhf@gmail.com or 
- duxiaoxiao@buaa.edu.cn. 

Citation:
Xiaoxiao Du, Gang Zhao, Ran Zhang, Wei Wang, Jiaming Yang. Numerical implementation for isogeometric analysis of thin-walled structures based on a Bézier extraction framework: nligaStruct. Thin-Walled Structures.
