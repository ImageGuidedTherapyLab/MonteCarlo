close all
clear all

% [usage]
% 
%    set BeamType  and spacing parameters and then type:
% 
%            multi_run for matlab command prompt
% 
% 
% 
% 
% set global variables
global BEAM_LOC   DIFF_LEN   DELTA_R   DELTA_Z 

DIFF_LEN = .5;    % length of diffusing tip [cm]
DELTA_R = .0117;   %  DELTA_R  = spacing along the r-axis [cm]
DELTA_Z = .0117;   %  DELTA_Z  = spacing along the z-axis [cm]
dimz = 257; % grid dimensions
dimr = 128; % grid dimensions

% beam radius
R_0 = .1 ;  % [cm]

% set function pointer representing type of beam
%                    diff_init : interstitial beam
%                    flat_init : Flat beam of radium R_0
%                    gauss_init: Gaussian beam with 1/exp(2) radius R_0
%                    point_init: Point source at the origin
BeamType = 4;
if(BeamType == 1)
   functer = @gauss_init ;
   namebase = 'Gauss_%04d';
   BEAM_LOC = 0.0;  % location in z-axis of center of diffusing tip [cm]
elseif(BeamType == 2)
   functer = @flat_init  ;
   namebase = 'Flat_%04d';
   BEAM_LOC = 0.0;  % location in z-axis of center of diffusing tip [cm]
elseif(BeamType == 3)
   functer = @point_init ;
   namebase = 'Point_%04d';
   BEAM_LOC = 0.0;  % location in z-axis of center of diffusing tip [cm]
elseif(BeamType == 4)
   functer = @diff_init  ;
   namebase = 'Diff_%04d';
   BEAM_LOC = 1.5;  % location in z-axis of center of diffusing tip [cm]
end


% create power grid
P = ones(dimz,dimr);
%c = floor(BEAM_LOC/DELTA_Z)+1; %grid location of the center of the diff tip
%l = floor(.5*DIFF_LEN/DELTA_Z); %number of grid points from center of diff tip to end
%w = floor(.4/DELTA_R);   %number of grid points out from laser for cooling
%P(c-l:c+l,1:1+w) = 0;
P = ones(dimz,dimr); % comment/uncomment for uniform power


% setup for parameter study 
mu_a = [.046,.44,.704,1.008,1.312,1.616,1.92,2.224,2.528,2.832,3.136,3.44,3.744,4.048,4.352,4.656,4.96,5.0];
mu_s = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25];
g = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,.99];
dist_r = [.1,.3,.5,.7,.9];
dist_z = [.1,.3,.5,.7,.9];
mu_a = [.046];           % comment/uncomment for single run 
mu_s = [14.7];           % comment/uncomment for single run
g    = [.7];             % comment/uncomment for single run
%dist_r = [.5];           % comment/uncomment for single run
%dist_z = [.5];           % comment/uncomment for single run

% set number of times to run
nrun = 10 ;


% Create distance grid for isotropic comparison
dist = zeros(dimz,dimr-1);
for iii = 1:dimz
    for jjj = 1:dimr-1
      % should be units of meters
      dist(iii,jjj) = .01* sqrt(  (            (jjj-.5)*DELTA_R ) ^ 2 ...
                                + ( BEAM_LOC - (iii-.5)*DELTA_Z ) ^ 2     ) ;
    end
end

ntotal = size(g,2)*size(mu_a,2)*size(mu_s,2)*size(dist_r,2)*size(dist_z,2);
icount = 0;  
for ig = 1:size(g,2)
  for ia = 1:size(mu_a,2)
    for is = 1:size(mu_s,2)
       for ir = 1:size(dist_r,2)
          for iz = 1:size(dist_z,2)
             disp(sprintf('\non grid %d of %d total grids',icount,ntotal))
             disp(sprintf('g = %f, mu_a = %f, mu_s = %f',...
                                                    g(ig),mu_a(ia),mu_s(is) ) )
             disp(sprintf('dist_r = %f, dist_r = %f', dist_r(ir),dist_z(iz) ) )
             % calculate isotropic source term for comparison
             mu_tr  = (mu_a(ia) *100)+ (mu_s(is)*100) * (1-g(ig)); % 1/meters
             mu_eff = sqrt(3*(mu_a(ia)*100)*mu_tr); % 1/meters
             %isotropic should be units of  [W/m^3]
             isotropic = 3/4*(mu_a(ia)*100)*mu_tr*pi*...
                             (P(:,1:dimr-1).*exp(-mu_eff*dist)).* dist.^-1;
             %initialize monte carlo grid
             heating4 = zeros(dimz,dimr-1);
             for i = 1:nrun
               disp(sprintf('  photon set %d of %d',i,nrun))
               [R,A,T,grid,fluence,HEATING] = ...
               scatter_simulation(10000,R_0,g(ig),mu_a(ia),mu_s(is),P,...
                                  functer,dist_r(ir),dist_r(iz));
               heating4 = heating4 + HEATING;
             end
             %monte carlo heating should be units of  [W/m^3]
             heating4 = heating4./nrun;  
             %output results
             name = sprintf(namebase,icount); icount = icount + 1 ;
             ascii_write(heating4,isotropic,name,DELTA_R,DELTA_Z,...
                         g(ig),mu_a(ia),mu_s(is),dist_r(ir),dist_z(iz));
         end
       end
     end
  end
end


