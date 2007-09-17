%This function is meant to be used in conjunction with scatt-simulation. It
%takes the critical angle, the ration of the refractive indicies, the depth
%of the tissue, the current amount of reflectance and transmission, the
%current weight of the photons, and the current z position and z directional 
%cosine of the photonsto figure out which photons went outside of the tissue 
%in order to calculate internal reflectance from the top and bottom of the 
%tissue. It returns the new z position and directional cosine for those
%photons that were internally reflected as well as the new reflectance and
%transmission counts.

function [PhotonsZN, PhotonsCZN, PhotonsWN,RN, TN] = internal_reflection(PhotonsZ,...
    PhotonsCZ, PhotonsW,R,T,d,theta_c,n)

   PhotonsZN = PhotonsZ;
   PhotonsCZN = PhotonsCZ;
   PhotonsWN = PhotonsW;
   
    ReflectionIndicies = find(PhotonsZ < 0); % Indicies of Photons that are out 
                                             % of tissue at top
    TransmittedIndicies = find(PhotonsZ > d); % Indicies of Photons that are out
                                              % of tissue at bottom
    Rtheta_t = acos(abs(PhotonsCZ(ReflectionIndicies))); % Incident angle of 
                                                         % Photons out of top 
    Ttheta_t = acos(abs(PhotonsCZ(TransmittedIndicies))); % Incident angle of
                                                          %Photons out of
                                                          %bottom
    
    
    % Deal with Photons out of top of Tissue
    if(numel(ReflectionIndicies) > 0)
        % Indicies of photons whose incident angle was below the critical angle
        
        ThetaRIndicies = find(Rtheta_t < theta_c);
    
        if(numel(ThetaRIndicies) >0)
       
            theta_a = asin(1/n.*sin(Rtheta_t(ThetaRIndicies))); % Angle predicted
                                                                %by Snell's
                                                                %Law 
            theta_p = Rtheta_t(ThetaRIndicies) + theta_a;
            theta_m = Rtheta_t(ThetaRIndicies) - theta_a;

            R_F =.5*((sin(theta_m).^2./sin(theta_p).^2)+...
                (tan(theta_m).^2./tan(theta_p).^2)); % Fresnal Reflection ratio
                    
            RN = R + sum(PhotonsW(ReflectionIndicies(ThetaRIndicies)).*(1-R_F));           
            PhotonsWN(ReflectionIndicies(ThetaRIndicies)) = PhotonsW(...
                ReflectionIndicies(ThetaRIndicies)).*R_F; 
        else
            RN = R;
        end
    
    PhotonsZN(ReflectionIndicies) = -PhotonsZ(ReflectionIndicies);
    PhotonsCZN(ReflectionIndicies) = -PhotonsCZ(ReflectionIndicies);   
    else
        RN=R;
    end
    
     % Deal with Photons out of bottom of Tissue
    if(numel(TransmittedIndicies) > 0)
       % Indicies of photons whose incident angle was below critical angle
      ThetaTIndicies = find(Ttheta_t < theta_c);
    
        if(numel(ThetaTIndicies) >0)
            
            theta_a = asin(1/n.*sin(Ttheta_t(ThetaTIndicies))); % Angle predicted
                                                                % by
                                                                % Snell's Law
            theta_p = Ttheta_t(ThetaTIndicies) + theta_a;
            theta_m = Ttheta_t(ThetaTIndicies) - theta_a;

            R_F =.5*((sin(theta_m).^2./sin(theta_p).^2)+...
                (tan(theta_m).^2./tan(theta_p).^2)); % Fresnal Reflection ratio
                    
            TN = T + sum(PhotonsW(TransmittedIndicies(ThetaTIndicies)).*(1-R_F));           
            PhotonsWN(TransmittedIndicies(ThetaTIndicies)) = PhotonsW(...
                TransmittedIndicies(ThetaTIndicies)).*R_F; 
        else
            TN = T;
        end
    
    PhotonsZ(TransmittedIndicies) = -PhotonsZ(TransmittedIndicies)+ones(size(...
        PhotonsZ(TransmittedIndicies)))*2*d;
    PhotonsCZ(TransmittedIndicies) = - PhotonsCZ(TransmittedIndicies);          
    else
        TN = T;
    end