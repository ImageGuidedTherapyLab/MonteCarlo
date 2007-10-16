%Andrea Hawkins
%BME 385J
%Program #7

%This program is a Monte Carlo simulation of light attenuation in tissue 
%incorporating absorption and non-isotropic scattering assuming the 
%Heyney-Greenstein distribution for the deflection angle.We are assuming that 
%the indicies of reflection are not matched at tissue boundary, n_1 = 1, 
%n_2 = 1.37. In this simulation, variable weight photons are used and the 
%photons can be initialized in multiple ways: flat beam, gaussian beam, 
%point. The program returns the percentage of photons that are reflected, 
%R, transmitted, T, and absorbed, A. It also returns a matrix representing 
%a grid (3mm x 1.5 mm) of the tissue on the r-z variables whose entries 
%represent the amount of photons absorbed in that volume. Each entry 
%represents a DELTA_R = .15 mm and DELTA_Z = .05 mm. The last column of 
%this matrix, however, represents all photons that went beyond the 3 mm 
%represented by the rest of the matrix. It also returns a matrix giving 
%the fluence in each of the elements representing the same region.The 
%expected value for R is .225 and for T is .015.

function [R,A,T,grid,fluence,heating] = scatter_simulation(num_photons,R_0,g,mu_a,mu_s,P,init,dist_r,dist_z)

% mu_a =  absorption coefficient [cm^-1]
% mu_s =  scattering coefficient [cm^-1]
% g    =  expected vale for cosine of angle scattering
% P    =  total Power [W] (a 2d grid of power values)
% dist_r,dist_z = coefficients to determine initial photon distribution for
%                  interstitial scatter initialization


tic; % function that keeps track f cputime
%initialization

% DELTA_R  = distance along the r-axis each element in grid represents [cm]
% DELTA_Z  = distance along the z-axis each element in grid represents [cm]
% BEAM_LOC = location in z-axis of center of diffusing tip [cm]
% DIFF_LEN = length of diffusing tip [cm]

 global DELTA_R DELTA_Z BEAM_LOC DIFF_LEN 

% num_photons        number of photons to be simulated
% R_0                parameter for beam type [cm]
% delta_s            distance photon moves between events
% phi                azimuthal angle
% theta              deflection angle
% PhotonsW           current weight of current photon
% wa                 portion of photon weight that is absorbed at current 
%                    event

 A = 0;            % counter for photons absorbed
 R = 0;            % counter for photons reflected
 T = 0;            % counter for photons transmitted
 mu_t = mu_a + mu_s; % attenuation coefficient [cm^-1]
 wa_multiplier = mu_a/mu_t; % fraction mu_a/mu_t
 e = .005;         % weight cutoff for photon
 m = 20;           % constant used to determine photon cutoff 
 d = 10;           % depth of tissue [cm]
 n1 = 1;           % refractive index of air
 n2 = 1.37;        % refractive index of tissue
 n = n1/n2;        % ratio to calculate Snell's Law
 theta_c = asin(n); % critical angle;
 dim1 = size(P,1);  % 1st dimension of grid
 dim2 = size(P,2);  % 2nd dimension of grid
 grid = zeros(dim1,dim2); % matrix representing the absorbed photons


% PhotonsX           1st coordinate of photon position
% PhotonsY           2nd coordinate of photon position
% PhotonsZ           3rd coordinate of photon position
% PhotonsCX          angle in 1st direction of photon travel
% PhotonsCY          angle in 2nd direction of photon travel
% PhotonsCZ          angle in 3rd direction of photon travel
% R_spec             specular reflection constant


% init is used as a function pointer to:
%                    diff_init : interstitial beam
%                    flat_init : Flat beam of radium R_0
%                    gauss_init: Gaussian beam with 1/exp(2) radius R_0
%                    point_init: Point source at the origin
                    
[PhotonsX,PhotonsY,PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,R_spec] = init(num_photons,BEAM_LOC,DIFF_LEN,R_0,n1,n2,dist_r,dist_z);

PhotonsW = ones(num_photons,1);

AliveIndicies = [1:num_photons];  % Indicies of photons with weight 
                                  % above cutoff

while(numel(AliveIndicies) > 0)
    
    
    % Move Photons
    [PhotonsX,PhotonsY,PhotonsZ] = move_photons(PhotonsX,PhotonsY,...
        PhotonsZ,PhotonsCX,PhotonsCY,PhotonsCZ,mu_t);
    
    % Check for internal reflection
    [PhotonsZ, PhotonsCZ, PhotonsW,R, T] = internal_reflection(...
        PhotonsZ, PhotonsCZ, PhotonsW,R,T,d,theta_c,n);
    
    
    % Calculate absorption and check for weight cutoff
    [PhotonsW, A, grid] = absorption_cutoff(PhotonsX,PhotonsY,...
        PhotonsZ,PhotonsW, A, wa_multiplier, m, e,grid,DELTA_R,DELTA_Z,...
        dim1,dim2);
    
    % Throw out Dead Photons
    AliveIndicies = find(PhotonsW > 0);
    
    PhotonsX = PhotonsX(AliveIndicies);
    PhotonsY = PhotonsY(AliveIndicies);
    PhotonsZ = PhotonsZ(AliveIndicies);
 
    PhotonsCX = PhotonsCX(AliveIndicies);
    PhotonsCY = PhotonsCY(AliveIndicies);
    PhotonsCZ = PhotonsCZ(AliveIndicies);

    PhotonsW = PhotonsW(AliveIndicies);
    
    % Update Angles
    [PhotonsCX, PhotonsCY, PhotonsCZ] = dir_cosine_update(PhotonsCX,...
        PhotonsCY, PhotonsCZ, g);
   
    
    end

R = R_spec + R/num_photons;
T = T/num_photons;
A = A/num_photons;


%Calculate fluence and shrink grid and fluence for output
[p,n] = size(grid);
Volume = zeros(p,n); % matrix holding the volume that each element
                     % in the grid matrix represents
Vol_mult = pi*DELTA_R^2*DELTA_Z;

for i=0:n-1
    Volume(:,i+1) = (2*i+1)*Vol_mult;
end

fluence = (grid.*P)./(num_photons.*Volume*mu_a); % matrix holding the
                                               % fluence values for
                                               %locations in the tissue
                                               %[W/cm^2]
                                               
fluence = fluence * 100^2;                 %Done to convert to [W/m^2]      
heating = mu_a*100 * fluence;              %Done such that in [W/m^3]
                                               

off_side = sum(grid(:,dim2)); % Number of photons that were absorbed 
                            % beyond the side boundary
                          
heating = heating(1:dim1,1:dim2-1);

toc % function that keeps track of cputime
