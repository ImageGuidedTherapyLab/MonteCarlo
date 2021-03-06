%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%  flat_init : is a beam with uniform irradiance over the circle centered
%              at the origin with radius R_0

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = flat_init(num_photons,beam_loc, diff_len,R_0,n1,n2,dist_r,dist_z)


disp('  initializing flat beam');

PhotonsR     = R_0*rand(num_photons,1).^dist_r;
PhotonsTheta = 2*pi*rand(num_photons,1);
PhotonsX = PhotonsR.*cos(PhotonsTheta);
PhotonsY = PhotonsR.*sin(PhotonsTheta);
PhotonsZ = diff_len*rand(num_photons,1).^dist_z + (beam_loc - .5*diff_len)...
      *ones(num_photons,1);


PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);

R_spec = (n1 - n2)^2/(n1+n2)^2; % specular reflection constant

