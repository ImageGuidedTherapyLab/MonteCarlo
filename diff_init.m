%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%This program initializes values of position in the z direction and the ...
%directional cosines simulating a diffusing laser tip.

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = diff_init(num_photons,beam_loc, diff_len,R_0,n1,n2,dist_r,dist_z)

global DELTA_R DELTA_Z 

disp('  initializing interstitial laser');


PhotonsTheta = 2*pi*rand(num_photons,1);
% uniformly distribute energy along radius i.e.
% see Welch section 4.5.1 equation (4.58)
% Assume independent random variable for the radial and axial photon
% distribution. the photon radial position is uniformly distributed over
% the area i.e.
% 
%   int(p(r),r,0,R_0) = int( 2*pi*r /pi/R_0^2* ,r,0,R_0)
% 
%   the probability that a photon with within [0,r] is equal to the 
%   probability that a uniform random variable is between [0,\xsi] (\xsi < 1)
% 
%   \eta = int(1,\eta,0,\eta)= int( 2*pi*r /pi/R_0^2* ,r,0,r) = r^2/R_0^2
%   r = R_0 * sqrt(\eta) 
% 
PhotonsRad   = R_0 *rand(num_photons,1).^.5;

PhotonsX = PhotonsRad.*cos(PhotonsTheta);
PhotonsY = PhotonsRad.*sin(PhotonsTheta);
PhotonsZ = diff_len*rand(num_photons,1);

PhotonsCX = PhotonsX;
PhotonsCY = PhotonsY;
PhotonsCZ = PhotonsZ;

magnitude = (PhotonsCX.^2 + PhotonsCY.^2  + PhotonsCZ.^2).^.5;
PhotonsCX = PhotonsCX.* magnitude.^-1 ;    
PhotonsCY = PhotonsCY.* magnitude.^-1 ;
PhotonsCZ = PhotonsCZ.* magnitude.^-1 ;

PhotonsZ = PhotonsZ + (beam_loc - .5*diff_len)*ones(num_photons,1);

R_spec = 0;       % specular reflection constant



