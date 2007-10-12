%This program initializes values of position in the z direction and the ...
%directional cosines simulating a diffusing laser tip.

function [PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ] = diff_init(num_photons,....
      beam_loc, diff_len)

PhotonsZ = diff_len*rand(num_photons,1).^(1/4) + (beam_loc - .5*diff_len)...
      *ones(num_photons,1);


PhotonsCZ = cos((rand(num_photons,1)).^7*pi);

PhotonsCY = rand(num_photons,1).*sqrt(1-PhotonsCZ.^2);

PhotonsCX = sqrt(1-PhotonsCZ.^2 - PhotonsCY.^2);


%PhotonsCX = zeros(num_photons,1);
%PhotonsCY = zeros(num_photons,1);
%PhotonsCZ = ones(num_photons,1);
