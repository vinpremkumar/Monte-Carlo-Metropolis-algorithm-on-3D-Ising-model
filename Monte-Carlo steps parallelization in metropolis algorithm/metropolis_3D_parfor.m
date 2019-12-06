function [spin_metro, Emean_output_after_metropolis, Mmean_output_after_metropolis] = metropolis_3D_parfor(par_process, spin, kT, J, visual_BHJgrid, montecarlo_steps)
%METROPOLIS The Metropolis algorithm.
%   spin = METROPOLIS(spin, kT, J) runs the Metropolis algorithm on a
%   configuration of spins with an coupling coefficient |J| at a
%   temperature |kT|. |spin| is a matrix of +/- 1's.

%   Copyright 2017 The MathWorks, Inc.

if(par_process == 1)
    parforArg = Inf;
else
    parforArg = 0;
end

numIters = montecarlo_steps;
Emean_output_after_metropolis = zeros(1,numIters);
Mmean_output_after_metropolis = zeros(1,numIters);
spin_metro                    = zeros([size(spin),numIters]);

% for (iter = 1 : numIters)
parfor(iter = 1: numIters, parforArg)
    spin_metro_temp = spin;
    Emean_output_after_metropolis_temp = 0;
    Mmean_output_after_metropolis_temp = 0;
    
    % Pick a random spin
    linearIndex = randi(numel(spin_metro_temp));
    [x, y, z]  = ind2sub(size(spin_metro_temp), linearIndex);
    
    % Find its nearest neighbors
    above = mod(x - 1 - 1, size(spin_metro_temp,1)) + 1;
    below = mod(x + 1 - 1, size(spin_metro_temp,1)) + 1;
    left  = mod(y - 1 - 1, size(spin_metro_temp,2)) + 1;
    right = mod(y + 1 - 1, size(spin_metro_temp,2)) + 1;
    front = mod(z - 1 - 1, size(spin_metro_temp,3)) + 1;
    back  = mod(z + 1 - 1, size(spin_metro_temp,3)) + 1;
    
    neighbors = [   spin_metro_temp(above,y,z); 
                    spin_metro_temp(x,left,z); 
                    spin_metro_temp(x,right,z); 
                    spin_metro_temp(below,y,z); 
                    spin_metro_temp(x,front,z); 
                    spin_metro_temp(x, back,z);];
    
    % Calculate energy change if this spin is flipped
    dE = 2 * J * spin_metro_temp(x, y, z) * sum(neighbors);
    
    % Boltzmann probability of flipping
    prob = exp(-dE / kT);
    
    % Spin flip condition
    if dE <= 0 || rand() <= prob
        spin_metro_temp(x, y, z) = - spin_metro_temp(x, y, z);
        
        if(visual_BHJgrid == 1)
            Emean_output_after_metropolis_temp = energyIsing_3D(spin_metro_temp, J);
            Mmean_output_after_metropolis_temp = magnetizationIsing_3D(spin_metro_temp);
            
            A = (spin_metro_temp+1)*128;
            A_xyz_dims = size(A);
            x_image_dim = 1:1:A_xyz_dims(1);
            y_image_dim = 1:1:A_xyz_dims(2);
            z_image_dim = 1:1:A_xyz_dims(3);
            cmap = jet(16);
            clim = [0,256];
            PATCH_3Darray(A,x_image_dim,y_image_dim,z_image_dim,cmap,clim,'col');
            xlabel(sprintf('T = %0.2f, M = %0.2f, E = %0.2f', kT, Mmean_output_after_metropolis_temp, Emean_output_after_metropolis_temp));
            set(gca,'YTickLabel',[],'XTickLabel',[]);
            axis square; 
            BlueAndRedColormap = [repmat([0,0,1], [1,1]); repmat([1,0,0], [1,1])];
            colormap(BlueAndRedColormap);
            colorbar;
            drawnow;
        end
    end
    spin_metro(:,:,:,iter) = spin_metro_temp;
    Emean_output_after_metropolis(1,iter) = Emean_output_after_metropolis_temp;
    Mmean_output_after_metropolis(1,iter) = Mmean_output_after_metropolis_temp;
end
spin_metro = spin_metro(:,:,:,end);
Emean_output_after_metropolis = Emean_output_after_metropolis(1,end);
Mmean_output_after_metropolis = Mmean_output_after_metropolis(1,end);
end
