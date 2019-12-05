%%--------------------------Monte-Carlo metropolis algorithm on 3D Ising model--------------------------
%% Exercise 9: https://www.mathworks.com/matlabcentral/fileexchange/62194-ising-model-and-metropolis-algorithm#functions_tab
%%-------------------------------------------------------------------------------------------------------
%% This file contains all the variables that can be tuned for achieving the result of 3D Ising model.
%% For more advanced tuning, please refer MCmetropolis_on_3DIsing.m file
%%-------------------------------------------------------------------------------------------------------

par_process         = 1;    % 1=yes, 0=no
visual_BHJgrid      = 1;    % 1=yes, 0=no % Enable in order to visualize every loop but set par_process to false.
                            % If parallel processing is enabled, visualization cannot occur every loop
numSpins_xDim       = 8;    % Number of spins in x dimension 
numSpins_yDim       = 8;    % Number of spins in y dimension
numSpins_zDim       = 8;    % Number of spins in z dimension
kT                  = [1];  % Temperature value in k (Boltzman constant) units (runs in parfor loop), make as an array
montecarlo_steps    = [100];% To vary the number of Monte-Carlo steps (runs in for loop), make as an array

MCmetropolis_on_3DIsing(par_process, visual_BHJgrid, numSpins_xDim, numSpins_yDim, numSpins_zDim, kT, montecarlo_steps);