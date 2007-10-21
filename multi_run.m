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

DIFF_LEN = .005;    % length of diffusing tip [cm]
DELTA_R = .0117;   %  DELTA_R  = spacing along the r-axis [cm]
DELTA_Z = .0117;   %  DELTA_Z  = spacing along the z-axis [cm]
dimz = 256; % grid dimensions
dimr = 129; % grid dimensions

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
   BEAM_LOC = 0.8;  % location in z-axis of center of diffusing tip [cm]
end


% create power grid
P = ones(dimz,dimr);
%c = floor(BEAM_LOC/DELTA_Z)+1; %grid location of the center of the diff tip
%l = floor(.5*DIFF_LEN/DELTA_Z); %number of grid points from center of diff tip to end
%w = floor(.4/DELTA_R);   %number of grid points out from laser for cooling
%P(c-l:c+l,1:1+w) = 0;
P = ones(dimz,dimr); % comment/uncomment for uniform power


% setup for parameter study 
% from Welch ch 8 Appendix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  in-vitro 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          lamdba     mu_t      mu_a      mu_s    mu_s(1-g)     g     mu_eff
%           (nm)    (cm^-1)   (cm^-1)   (cm^-1)   (cm^-1)             (cm^-1)
%          --------------------------------------------------------------------
%###Brain###   
% Calf      633        -       .19         -         -          -     2.5-3.4
%           1064       -       .36         -         -          -       2.5
%           1320       -       .84         -         -          -       4.0
% Cat       488        -         -         -         -          -       10.9
%           514.5      -         -         -         -          -      13.3
%           630        -         -         -         -          -     5.3-8.9
% Human     488        -         -         -         -          -    14.0-25.0
%  (Adult)  514        -         -         -         -          -    14.0-16.7
%           660        -         -         -         -          -     7.0-12.5
%           1060       -         -         -         -          -     2.3- 3.4
%           630        -     0.3-1.0       -      30.0-40.0     -       8.3
%  (Neonate)488        -         -         -         -          -     5.9-7.9
%           514        -         -         -         -          -     5.8-9.0
%           660        -         -         -         -          -     2.5-3.3
%           1060       -         -         -         -          -     1.1-1.4
%  (Adult white matter)
%          632.8       -      2.2+-.2   532+-41   91+-5     0.82+-.01  23.8+-1.4
%          1064        -      3.2+-.4   469+-34 60.4+-2.55  0.87+-.007 24.1+-1.7
%  (Adult grey matter)
%          632.8       -      2.7+-.2   354+-37 20.6+-2     0.94+-.004 13.3+-0.9
%          1064        -      5.0+-.5   134+-34 11.8+-0.9   0.90+-.007 15.7+-1.3
% Pig
%          633     0.26-0.64     -       52-57       -        0.945     4.3-14.2
%
%
%
dist_r = [.5];
dist_z = [10];
mu_a = [.046];
mu_s = [3;14.7;21];
g = [.7;.85;.99];
%dist_r = [.1;.3;.5;.7;.9];
%dist_z = [.1;.3;.5;.7;.9];
%mu_a = [.046];           % comment/uncomment for single run 
%mu_s = [14.7];           % comment/uncomment for single run
%g    = [.7];             % comment/uncomment for single run
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

ntotal = size(g,1)*size(mu_a,1)*size(mu_s,1)*size(dist_r,1)*size(dist_z,1);
icount = 0;  
for ig = 1:size(g,1)
  for ia = 1:size(mu_a,1)
    for is = 1:size(mu_s,1)
       for ir = 1:size(dist_r,1)
          for iz = 1:size(dist_z,1)
             disp(sprintf('\non grid %d of %d total grids',icount,ntotal))
             disp(sprintf('g = %f, mu_a = %f, mu_s = %f',...
                                                    g(ig),mu_a(ia),mu_s(is) ) )
             disp(sprintf('dist_r = %f, dist_z = %f', dist_r(ir),dist_z(iz) ) )
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
                                  functer,dist_r(ir),dist_z(iz));
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


