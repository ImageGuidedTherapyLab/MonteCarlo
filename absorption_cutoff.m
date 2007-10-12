%This function is meant to be used in conjunction with scatter_simulation.
%It takes the current positions of the photons, the current weights of the
%photons, the current absorption count, the weight-absorption ratio, weight
%cutoff and weight check number, the absorption grid and the size of the
%grid elements to check if the weights of the photons are small enough to
%be considered gone, and put into the grid the amount of each photons that
%was absorbed in the appropriate position. It returns the new photon
%weights, absorption count and absorption grid.

function [PhotonsWN, AN,grid] = absorption_cutoff(PhotonsX,PhotonsY,PhotonsZ,...
    PhotonsW, A, wa_mult, m, e,grid,delta_r,delta_z,dim1,dim2)
      
    % Calculate amount of photon that is absorbed
    wa = PhotonsW*wa_mult;
    [grid] = absorption_grid(PhotonsX,PhotonsY,PhotonsZ,wa,grid,delta_r,delta_z,dim1,dim2);
    AN = A + sum(wa);
         
    % Calculate weight of photon that scatters
    PhotonsW = PhotonsW - wa;
         
    % Check for weight cutoff
    WeightIndicies = find(PhotonsW < e); % Indicies of photons whose weight is 
                                         % below cutoff
    RN = rand(numel(WeightIndicies),1); % Random Number matrix corresponding 
                                        % to photons whose weight is below cutoff
    RNIndicies = find(RN < 1/m); % Indicies where the random number was less 
                                 % than 1/m
    
    PhotonsW(WeightIndicies(RNIndicies)) = PhotonsW(WeightIndicies(RNIndicies))*m; 
    
    RNIndicies = find(RN >= 1/m); % Indicies where the random number was greater 
                                  % than 1/m
    
    PhotonsW(WeightIndicies(RNIndicies)) = 0; 
    
    PhotonsWN = PhotonsW;
