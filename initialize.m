%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%This program initializes values of position in the z direction and the ...
%directional cosines simulating a diffusing laser tip.

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ] = diff_init(num_photons,beam_loc, diff_len)

PhotonsX = zeros(num_photons,1);
PhotonsY = zeros(num_photons,1);
PhotonsZ = diff_len*rand(num_photons,1).^(1/4) + (beam_loc - .5*diff_len)...
      *ones(num_photons,1);


PhotonsCZ = cos((rand(num_photons,1)).^7*pi);

PhotonsCY = rand(num_photons,1).*sqrt(1-PhotonsCZ.^2);

PhotonsCX = sqrt(1-PhotonsCZ.^2 - PhotonsCY.^2);



%  flat_init : is a beam with uniform irradiance over the circle centered
%              at the origin with radius R_0

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ] = flat_init(num_photons,beam_loc, diff_len)


PhotonsX = R_0*sqrt(rand(num_photons,1));
PhotonsY = zeros(num_photons,1);
PhotonsZ = zeros(num_photons,1);

PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);

%  gauss_init: is a beam of irradiance that goes E(z) =
%              E_0*exp(-2*r^2/R_0^2) where R_0 is the 1/exp(2) radius
function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ] = gauss_init(num_photons,beam_loc, diff_len)


PhotonsX = R_0/sqrt(2)*sqrt(-log(rand(num_photons,1)));
PhotonsY = zeros(num_photons,1);
PhotonsZ = zeros(num_photons,1);

PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);

%  point_init: assumes all photons are initiated with x=0.
function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ] = point_init(num_photons,beam_loc, diff_len)


PhotonsX = zeros(num_photons,1);
PhotonsY = zeros(num_photons,1);
PhotonsZ = zeros(num_photons,1);

PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);

