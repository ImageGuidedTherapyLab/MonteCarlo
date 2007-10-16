%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%This program initializes values of position in the z direction and the ...
%directional cosines simulating a diffusing laser tip.

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = diff_init(num_photons,beam_loc, diff_len,R_0,n1,n2,dist_r,dist_z)

disp('initializing interstitial laser');
PhotonsR     = R_0*rand(num_photons,1).^dist_r;
PhotonsTheta = 2*pi*rand(num_photons,1);
PhotonsX = PhotonsR.*cos(PhotonsTheta);
PhotonsY = PhotonsR.*sin(PhotonsTheta);
PhotonsZ = diff_len*rand(num_photons,1).^dist_z + (beam_loc - .5*diff_len)...
      *ones(num_photons,1);


PhotonsCX = zeros(num_photons,1);
PhotonsCY = zeros(num_photons,1);
PhotonsCZ = ones(num_photons,1);

R_spec = 0;       % specular reflection constant



