%Initialize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%This program initializes values of position in the z direction and the ...
%directional cosines simulating a diffusing laser tip.

function [PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = diff_init(num_photons,beam_loc, diff_len,R_0,n1,n2,dist_r,dist_z)

global DELTA_R DELTA_Z 

disp('  initializing interstitial laser');

% create a list of all possible positions with the boundaries of the limacon
Naxl=floor(diff_len/DELTA_Z);
Nrad=floor(diff_len/DELTA_R);
options = optimset('fsolve');
options = optimset(options,'MaxFunEvals',1000);
xguess = [1  1/sqrt(2) 0;1 1 0; 0 0 1] \ [R_0;diff_len;pi/4];
[x,f,exit] = fsolve(@(x) limacon(x,R_0,diff_len),xguess,options);
a = x(1);
b = x(2);
thetamax = x(3);
%Nrad = floor(log(1-diff_len/(DELTA_R*(dist_r-1))*(1-dist_r))/log(dist_r) - 1)
%radlist = diff_len*ones(Nrad+1,1);
%for i = 0:Nrad
%  for j = i+1:Nrad+1
%    radlist(j) = radlist(j) - dist_r^i* (DELTA_R*(dist_r-1));
%  end
%end
radlist = diff_len*[1/Nrad/2:1/Nrad:1];
axllist = diff_len*([   -1  :1/Naxl:1].^dist_z);


poslist=[];
for i = 1:size(radlist,2)
   for j = 1: size(axllist,2)
       r= radlist(i) ;
       z= axllist(j) ;
       rad = sqrt(z^2 + r^ 2);
       theta = atan2(r,z);
       limaconrad = a + b * cos(theta);
       if(rad <= limaconrad)
          poslist=[poslist;r,z];
       end
   end
end

RandPos = ceil(size(poslist,1).*rand(num_photons,1).^dist_r);

PhotonsTheta = 2*pi*rand(num_photons,1);

PhotonsX = poslist(RandPos,1).*cos(PhotonsTheta);
PhotonsY = poslist(RandPos,1).*sin(PhotonsTheta);
PhotonsZ = poslist(RandPos,2);

%PhotonsX = zeros(num_photons,1);
%PhotonsY = CardioidY;
%%PhotonsY = R_0*rand(num_photons,1).^dist_r;
%PhotonsZ = -CardioidX;

PhotonsCX = PhotonsX;
PhotonsCY = PhotonsY;
PhotonsCZ = PhotonsZ;


magnitude = (PhotonsCX.^2 + PhotonsCY.^2  + PhotonsCZ.^2).^.5;
PhotonsCX = PhotonsCX.* magnitude.^-1 ;    
PhotonsCY = PhotonsCY.* magnitude.^-1 ;
PhotonsCZ = PhotonsCZ.* magnitude.^-1 ;

%PhotonsCX = ones(num_photons,1);
%PhotonsCY = ones(num_photons,1);
%PhotonsCZ = zeros(num_photons,1);

%PhotonsZ = PhotonsZ + (beam_loc - .5*diff_len)*ones(num_photons,1);
PhotonsZ = PhotonsZ + beam_loc*ones(num_photons,1);



R_spec = 0;       % specular reflection constant



