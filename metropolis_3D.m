function [spin_metro, Emean_output_after_metropolis, Mmean_output_after_metropolis] = metropolis(spin_metro, kT, J, visual_BHJgrid, montecarlo_steps)
%METROPOLIS The Metropolis algorithm.
%   spin = METROPOLIS(spin, kT, J) runs the Metropolis algorithm on a
%   configuration of spins with an coupling coefficient |J| at a
%   temperature |kT|. |spin| is a matrix of +/- 1's.

%   Copyright 2017 The MathWorks, Inc.

numIters = montecarlo_steps;
Emean_output_after_metropolis = 0;
Mmean_output_after_metropolis = 0;

for iter = 1 : numIters
    % Pick a random spin
    linearIndex = randi(numel(spin_metro));
    [x, y, z]  = ind2sub(size(spin_metro), linearIndex);
    
    % Find its nearest neighbors
    above = mod(x - 1 - 1, size(spin_metro,1)) + 1;
    below = mod(x + 1 - 1, size(spin_metro,1)) + 1;
    left  = mod(y - 1 - 1, size(spin_metro,2)) + 1;
    right = mod(y + 1 - 1, size(spin_metro,2)) + 1;
    front = mod(z - 1 - 1, size(spin_metro,3)) + 1;
    back  = mod(z + 1 - 1, size(spin_metro,3)) + 1;
    
    neighbors = [   spin_metro(above,y,z); 
                    spin_metro(x,left,z); 
                    spin_metro(x,right,z); 
                    spin_metro(below,y,z); 
                    spin_metro(x,front,z); 
                    spin_metro(x, back,z);];
    
    % Calculate energy change if this spin is flipped
    dE = 2 * J * spin_metro(x, y, z) * sum(neighbors);
    
    % Boltzmann probability of flipping
    prob = exp(-dE / kT);
    
    % Spin flip condition
    if dE <= 0 || rand() <= prob
        spin_metro(x, y, z) = - spin_metro(x, y, z);
        
        if(visual_BHJgrid == 1)
            Emean_output_after_metropolis = energyIsing_3D(spin_metro, J);
            Mmean_output_after_metropolis = magnetizationIsing_3D(spin_metro);
            
            A = (spin_metro+1)*128;
            A_xyz_dims = size(A);
            x_image_dim = 1:1:A_xyz_dims(1);
            y_image_dim = 1:1:A_xyz_dims(2);
            z_image_dim = 1:1:A_xyz_dims(3);
            cmap = jet(16);
            clim = [0,256];
            PATCH_3Darray(A,x_image_dim,y_image_dim,z_image_dim,cmap,clim,'col');
            xlabel(sprintf('T = %0.2f, M = %0.2f, E = %0.2f', kT, Mmean_output_after_metropolis, Emean_output_after_metropolis));
            set(gca,'YTickLabel',[],'XTickLabel',[]);
            axis square; 
            BlueAndRedColormap = [repmat([0,0,1], [1,1]); repmat([1,0,0], [1,1])];
            colormap(BlueAndRedColormap);
            colorbar;
            drawnow;
        end
    end
end
end
