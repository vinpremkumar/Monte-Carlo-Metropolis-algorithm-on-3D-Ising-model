%%--------------------------Monte-Carlo metropolis algorithm on 3D Ising model--------------------------
%% Exercise 9: https://www.mathworks.com/matlabcentral/fileexchange/62194-ising-model-and-metropolis-algorithm#functions_tab
%%-------------------------------------------------------------------------------------------------------

function MCmetropolis_on_3DIsing(par_process_var, visual_BHJgrid_var, numSpins_xDim_var, numSpins_yDim_var, numSpins_zDim_var, kT_var, montecarlo_steps_var)

clc;
close all;

%Parallel processing to speed it up? = no loop-by-loop visualization
delete(gcp('nocreate'));
par_process = par_process_var;    % edit this; 1=yes, 0=no
if(par_process == 1)
    parforArg = Inf;
else
    parforArg = 0;
end

%Visualizing BHJ grid (Metropolis algorithm implemented to 3D Ising model)
visual_BHJgrid = visual_BHJgrid_var; % edit this; 1=yes, 0=no % Enable in order to visualize every loop but set par_process to false.
if(par_process == 1)
    visual_BHJgrid = 0; %If parallel processing is enabled, visualization cannot occur every loop
end

%Initialization
%Edit this for making n x n x n BHJ matrix
numSpins_xDim = numSpins_xDim_var; %Number of spins in x dimension
numSpins_yDim = numSpins_yDim_var; %Number of spins in y dimension
numSpins_zDim = numSpins_zDim_var; %Number of spins in z dimension

probSpinUp = 0.5;
J = 1;

% Temperatures to sample
kTc = 2*J / log(1+sqrt(2)); % Curie temperature
% Enable if you require a graph similar to the example program at https://www.mathworks.com/matlabcentral/fileexchange/62194-ising-model-and-metropolis-algorithm#functions_tab
% numTemps = 2^4;           
% kT = linspace(0, 2*kTc, numTemps); %change this for Temperature
kT = kT_var; %Temperature; For more than 1, make an array
numTemps = length(kT); %number of temperature samples;

montecarlo_steps = montecarlo_steps_var; %change this to vary the number of Monte-Carlo steps
numMCsteps = length(montecarlo_steps);

% Preallocate to store results
Emean = zeros(numTemps,numMCsteps);
Mmean = zeros(numTemps,numMCsteps);
spin = zeros(numSpins_xDim,numSpins_yDim,numSpins_zDim,numTemps,numMCsteps);

% Replace 'for' with 'parfor' to run in parallel with Parallel Computing Toolbox.
parfor (tempIndex = 1 : numTemps, parforArg)
    for (tempIndex2 = 1 : numMCsteps)
        spin_temp = initSpins_3D(numSpins_xDim, numSpins_yDim, numSpins_zDim, probSpinUp);
        [spin_temp, Emean_output_after_metropolis, Mmean_output_after_metropolis] = metropolis_3D(spin_temp, kT(tempIndex), J, visual_BHJgrid, montecarlo_steps(tempIndex2));
        
        if(visual_BHJgrid == 1)
            Emean(tempIndex,tempIndex2) = Emean_output_after_metropolis;
            Mmean(tempIndex,tempIndex2) = Mmean_output_after_metropolis;
        else
            Emean(tempIndex,tempIndex2) = energyIsing_3D(spin_temp, J);
            Mmean(tempIndex,tempIndex2) = magnetizationIsing_3D(spin_temp);
        end
        spin(:,:,:,tempIndex,tempIndex2) = spin_temp;
    end
end

%Plot the Ising model of BHJ
close all;
for tempIndex = 1:numTemps
    for tempIndex2 = 1:numMCsteps
        figure;
        A = (spin(:,:,:,tempIndex,tempIndex2)+1)*128;
        A_xyz_dims = size(A);
        x_image_dim = 1:1:A_xyz_dims(1);
        y_image_dim = 1:1:A_xyz_dims(2);
        z_image_dim = 1:1:A_xyz_dims(3);
        cmap = jet(16);
        clim = [0,256];
        PATCH_3Darray(A,x_image_dim,y_image_dim,z_image_dim,cmap,clim,'col');
        xlabel(sprintf('T = %0.2f, M = %0.2f, E = %0.2f, Monte-Carlo steps = %d', kT(tempIndex), Mmean(tempIndex,tempIndex2), Emean(tempIndex,tempIndex2), montecarlo_steps(tempIndex2)));
        set(gca,'YTickLabel',[],'XTickLabel',[]);
        axis square;
        BlueAndRedColormap = [repmat([0,0,1], [1,1]); repmat([1,0,0], [1,1])];
        colormap(BlueAndRedColormap);
        colorbar;
        view(-40.3,37.2);
        drawnow;
    end
end

%Let's plot the energy vs. temperature and its moving mean and median using the functions movmean and movmedian.
figure;
plot(kT / kTc, Emean, '.');
hold on;
window = (2^-3)*numTemps - 1;
grid on;
if numTemps >= 2^4
    plot(kT / kTc, movmean(  Emean, window));
    plot(kT / kTc, movmedian(Emean, window));
    hold off;
    title('Mean Energy Per Spin vs Temperature');
    xlabel('kT / kTc');
    ylabel('<E>');
    legend('raw',...
        [num2str(window) ' pt. moving mean'],...
        [num2str(window) ' pt. moving median'],...
        'Location', 'NorthWest');
end
%Notice the inflection point at a particular temperature, known as the Curie temperature.
%This significant change is referred to as a phase change of the ferromagnet.


%Let's plot the magnetization vs. temperature.
figure;
plot(kT / kTc, Mmean, '.');
grid on;
title('Magnetization vs Temperature');
xlabel('kT / kTc');
ylabel('M');
%Notice that when started from a distribution of equally like up and down spins,
%the system is equally likely to magnetize in the up and down orientations.


%Let's look at the absolute value of the magnetization in the same way we looked at the energy.
figure;
plot(kT / kTc, abs(Mmean), '.');
hold on;
window = (2^-3)*numTemps - 1;
grid on;
if numTemps >= 2^4
    plot(kT / kTc, movmean(  abs(Mmean), window));
    plot(kT / kTc, movmedian(abs(Mmean), window));
    hold off;
    title('Magnetization vs Temperature');
    xlabel('kT / kTc');
    ylabel('|M|');
    legend('raw',...
        [num2str(window) ' pt. moving mean'],...
        [num2str(window) ' pt. moving median'],...
        'Location', 'NorthEast');
end
%Notice that, like the energy, the magnetization changes most significantly at the Curie temperature.