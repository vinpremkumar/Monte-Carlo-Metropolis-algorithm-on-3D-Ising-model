# Monte-Carlo-Meteropolis-algorithm-on-3D-Ising-model
Code was edited from : MathWorks Physics Team (2019). Ising Model and Metropolis Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/62194-ising-model-and-metropolis-algorithm), MATLAB Central File Exchange. Retrieved December 5, 2019.

Ising Model and Metropolis Algorithm by the  MathWorks Physics Team is an excellent code for simulating 2D Ising model using Monte-Carlo Meteropolis algorithm. However, the modelling of the 3D Ising model was left as an *Exercise (no. 9)*. In this repository, I have coded the 3D Ising model, as well as visualized it.

There are two ways to parallelize the program.
    1. Parallelize the Temperature iteration
    2. Parallelize the Monte-Carlo steps in the Metropolis algorithm.
Both parallelization methods are provided in this repository. After downloading the common files (which are not inside any specific folder), copy the files from the type of parallelization that suites your requirement.

After completing the above step, run **_runthis_containsvariablesfortuning.m_** file for most of the basic tuning of Monte-Carlo metropolis algorithm on 3D Ising model.
For more advanced tuning, edit **_MCmetropolis_on_3DIsing_x_parallelized.m_** file.

Additional files required for visualization:
**_PATCH_3Darray.m_** (https://www.mathworks.com/matlabcentral/fileexchange/28497-plot-a-3d-array-using-patch?focused=5224890&tab=function)
