%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%  gauss_init: is a beam of irradiance that goes E(z) =
%              E_0*exp(-2*r^2/R_0^2) where R_0 is the 1/exp(2) radius
function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = gauss_init(num_photons,beam_loc, diff_len,R_0,n1,n2,dist_r,dist_z)


disp('  initializing gauss beam');
PhotonsX = R_0/sqrt(2)*sqrt(-log(rand(num_photons,1)));
PhotonsY = zeros(num_photons,1);
PhotonsZ = zeros(num_photons,1);

PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);
R_spec = (n1 - n2)^2/(n1+n2)^2; % specular reflection constant

