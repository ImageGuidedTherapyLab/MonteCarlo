%This function is meant to used in conjunction with scatter_simulation. It
%takes the total attentuation coefficient, the current positions and the 
%directional cosines of the photons and returns the new positions of the
%photons.

function [PhotonsXN,PhotonsYN,PhotonsZN] = move_photons(PhotonsX,PhotonsY,...
    PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,mu_t)
    
    [m,n] = size(PhotonsZ);
    delta_s = -log(1-rand(m,n))/(mu_t);
    
    % New position of photons

    PhotonsXN = PhotonsX + delta_s.*PhotonsCX;
    PhotonsYN = PhotonsY + delta_s.*PhotonsCY;
    PhotonsZN = PhotonsZ + delta_s.*PhotonsCZ;