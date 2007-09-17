%This function is to be used with scatter_simulation and absorption_cutoff.
%It takes the current positions of the photons, the amount absorbed of each
%photons, the absorption grid, and the size of the grid elements to update
%the absorption grid.

function [grid] = absorption_grid(PhotonsX,PhotonsY,PhotonsZ,wa,grid,delta_r,...
    delta_z)

s = numel(wa); % Number of photons

r = sqrt(PhotonsX.^2 + PhotonsY.^2); % radial position of photons

for i = 1:s
    Z = floor(PhotonsZ(i)/delta_z)+1;
    R = floor(r(i)/delta_r) + 1;
    if( Z <= 512 & Z >= 0)
        if(R <= 512)
          grid(Z,R) = grid(Z,R) + wa(i);
        else
          grid(Z,513) = grid(Z,513) + wa(i);
        end
    end
end
        
