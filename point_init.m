%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.


%  point_init: assumes all photons are initiated with x=0.
function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = point_init(num_photons,beam_loc, diff_len,R_0,n1,n2,dist_r,dist_z)


disp('  initializing point source');
PhotonsX = zeros(num_photons,1);
PhotonsY = zeros(num_photons,1);
PhotonsZ = zeros(num_photons,1);

PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);

PhotonsCZ = cos((rand(num_photons,1))*pi);

PhotonsCY = rand(num_photons,1).*sqrt(1-PhotonsCZ.^2);

PhotonsCX = sqrt(1-PhotonsCZ.^2 - PhotonsCY.^2);

R_spec = (n1 - n2)^2/(n1+n2)^2; % specular reflection constant

