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

DIFF_LEN = 0.5;    % length of diffusing tip [cm]
DELTA_R = .0117;   %  DELTA_R  = spacing along the r-axis [cm]
DELTA_Z = .0117;   %  DELTA_Z  = spacing along the z-axis [cm]
dimz = 256; % grid dimensions
dimr = 129; % grid dimensions

% beam radius
R_0 = 0.1;  % [cm]

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          lamdba     mu_t      mu_a      mu_s    mu_s(1-g)     g     mu_eff
%           (nm)    (cm^-1)   (cm^-1)   (cm^-1)   (cm^-1)             (cm^-1)
%          --------------------------------------------------------------------
%###H20###   
%           633        -       .005       <.003      -          -       -    
%           810        -       .05        <.003      -          -       -
%           1064       -       .50        <.003      -          -       -
% setup for parameter study 
% from Welch ch 8 Appendix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN-VITRO 
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
% Pig      633     0.26-0.64     -       52-57       -        0.945     4.3-14.2
%###Prostate###   
% Canine    355        -       9.2         -       47.1         -       39.4
%           532        -       2.4         -       23.8         -       13.7
%           1064       -      0.04         -       13.0         -       1.25
% Human     633       9.6     0.7+-2       -      8.6+-.5      0.0    4.3+-.5
%           612        -         -         -         -          -       6.3
%           594        -         -         -         -          -       9.1
%           543        -         -         -         -          -      18.2
%           1064       -      1.47+-.24  47+-13     6.43      0.862      -   
%###Blood###   
% Hct=.47  partially oxygenated
%           760      2840      15.5     2820       7.9        0.9972    33.0
%###Tumors###   
% Brain,intracranial,Human
%           488        -         -         -         -          -     7.1-20.0
%           513        -         -         -         -          -     7.1-20.0
%           635        -         -         -         -          -     5.9- 3.9
%           1060       -         -         -         -          -     3.3- 1.9
% Brain,Human
%           630        -         -         -         -          -     3.8- 8.3
% Prostate,rat,r3327-AT
%           630       271      0.49       270     8.1-5.4    .97-.98  3.6- 2.9
% Sarcoma,rat
%           630        -         -         -         -          -       2.3
%         514.5        -         -         -         -          -       4.8
% fibrosarcoma,rat
%           630        -         -         -         -          -      4.4-9.8
%###Myocardium (specialized cardiac muscle cells)###   
% Dog native (<=45degC) coagulated at 45 -> 74degC
%           632.8      -     4.2+-.2       -         -       .71+-.02     -    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IN-VIVO 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%          lamdba     mu_t      mu_a      mu_s    mu_s(1-g)     g     mu_eff
%           (nm)    (cm^-1)   (cm^-1)   (cm^-1)   (cm^-1)             (cm^-1)
%          --------------------------------------------------------------------
%###Brain###   
% Cat    405-410       -         -         -         -          -      44.1
%          545         -         -         -         -          -      34.4
%          577         -         -         -         -          -      25.9
%          631         -         -         -         -          -     5.0-9.8
% Human    630         -         -         -         -          -     4.8-10.0
% Pig      630         -         -         -         -          -     3.7-4.5
% Brain Tumors
%          630         -         -         -         -          -     2.2-6.6
%###Prostate (rat)###   
% R3327-At 630         -       0.97        -         13.2       -       -
%          789         -       0.6         -         6.2        -       -
% R3327-H  630         -       1.53        -         11.4       -       -
%          789         -       0.96        -         7.86       -       -
%###Tumors###   
% Human retinoblastoma (athymic mice)
%          488/514     -         -         -         -          -       6.25
%          630         -         -         -         -          -       3.03
%          668         -         -         -         -          -       2.8
%          1064        -         -         -         -          -       1.3
% Mammary carcinoma (C3H/HEJ mice)
%          488/514     -         -         -         -          -       9.1
%          630         -         -         -         -          -       5.0
%          668         -         -         -         -          -       4.3
%          1064        -         -         -         -          -       2.7
% B16 melanotic melanoma (C57/B16)
%          630         -         -         -         -          -      20.0
%          668         -         -         -         -          -      20.0
%          1064        -         -         -         -          -       5.0
%
%
dist_r = [0.5];  % distribute uniformly in circular radius
dist_z = [1.0];   % uniformly distribute over length of laser tip
% mu_a arranged so that the fastest runs are done first...
mu_a = [ 15.5 ;... %blood
        1.23;... % human prostate
        .36;...  %calf brain  in-vitro
        .04 ]; %canine prostate in-vitro (also close to water)
mu_s = [3;...      % agar (as a guess ~ water *100)
        47.0;...   % human prostate in-vitro
        435.0;...  % brain adult white matter 
        2820];     % blood
g    = [.71;...    % myocardium
        .862;...   % human prostate in-vitro
        .97];      % prostate rat tumor
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
                               + ( BEAM_LOC  - (iii-.5)*DELTA_Z ) ^ 2     ) ;
    end
end

ntotal = size(g,1)*size(mu_a,1)*size(mu_s,1)*size(dist_r,1)*size(dist_z,1);
icount = 0;  
for is = 1:size(mu_s,1)
  for ia = 1:size(mu_a,1)
    for ig = 1:size(g,1)
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
           isotropic = 3/4*(mu_a(ia)*100)*mu_tr/pi*...
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


