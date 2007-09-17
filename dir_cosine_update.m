%This function is meant to be used in conjunction with the
%scatter_simulation. It takes the current directional cosines and the
%expected value of the scattering angle, and returns the new directional
%cosines. 

function [PhotonsCXP, PhotonsCYP, PhotonsCZP] = dir_cosine_update(...
    PhotonsCX, PhotonsCY, PhotonsCZ, g)
    

    [m,n] = size(PhotonsCX);
    % Update angles
    phi = rand(m,n)*2*pi;
    if(g == 0)
        theta = acos(2*rand(m,n) - 1);
    else
        theta = acos((1/(2*g)).*(1+g^2-((1-g^2)./(1-g+2*g*rand(m,n))).^2));
    end
    
    PhotonsCXP = PhotonsCX;
    PhotonsCYP = PhotonsCY;
    PhotonsCZP = PhotonsCZ;
    
    % Take care of Indicies where cz is close to 1
    AngleIndicies = find(PhotonsCZ > .99999);
    Mult = sin(theta(AngleIndicies));
    PhotonsCXP(AngleIndicies) = Mult.*cos(phi(AngleIndicies));
    PhotonsCYP(AngleIndicies) = Mult.*sin(phi(AngleIndicies));
    PhotonsCZP(AngleIndicies) = PhotonsCZ(AngleIndicies).*cos(theta(...
        AngleIndicies))./abs(PhotonsCZ(AngleIndicies));
    
    
    % Take care of Indicies where cz is not close to 1
    
    AngleIndicies = find(PhotonsCZ <= .99999);
    
    Multiplier = sin(theta(AngleIndicies))./sqrt(1-PhotonsCZ(AngleIndicies).^2);
    PhotonsCXP(AngleIndicies) = Multiplier.*(PhotonsCX(AngleIndicies).*...
        PhotonsCZ(AngleIndicies).*cos(phi(AngleIndicies)) - PhotonsCY(AngleIndicies).*...
        sin(phi(AngleIndicies)))+PhotonsCX(AngleIndicies).*cos(theta(AngleIndicies));
   
    PhotonsCYP(AngleIndicies) = Multiplier.*(PhotonsCY(AngleIndicies).*...
        PhotonsCZ(AngleIndicies).*cos(phi(AngleIndicies)) + PhotonsCX(AngleIndicies).*...
        sin(phi(AngleIndicies))) + PhotonsCY(AngleIndicies).*cos(theta(AngleIndicies));

    PhotonsCZP(AngleIndicies) = -sin(theta(AngleIndicies)).*cos(phi(AngleIndicies)).*...
        sqrt(1-PhotonsCZ(AngleIndicies).*PhotonsCZ(AngleIndicies))+PhotonsCZ(...
        AngleIndicies).*cos(theta(AngleIndicies));

     
