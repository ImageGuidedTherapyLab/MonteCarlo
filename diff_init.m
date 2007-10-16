%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%This program initializes values of position in the z direction and the ...
%directional cosines simulating a diffusing laser tip.

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = diff_init(num_photons,beam_loc, diff_len,R_0,n1,n2)

disp('initializing interstitial laser');
PhotonsX = zeros(num_photons,1);
PhotonsY = zeros(num_photons,1);
PhotonsZ = diff_len*rand(num_photons,1).^(1/4) + (beam_loc - .5*diff_len)...
      *ones(num_photons,1);


PhotonsCZ = cos((rand(num_photons,1)).^7*pi);

PhotonsCY = rand(num_photons,1).*sqrt(1-PhotonsCZ.^2);

PhotonsCX = sqrt(1-PhotonsCZ.^2 - PhotonsCY.^2);
R_spec = 0;       % specular reflection constant



