%Initizlize is a program meant to be used with scatter_simulation. It takes
%the numbers of photons needing to be simulated and returns the appropriate
%starting value for the 1st coordinate of position.

%beam_type is meant to be an integer 1,2 or 3
%  1: is a beam with uniform irradiance over the circle centered at the
%         origin with radius R_0
%  2: is a beam of irradiance that goes E(z) =
%             E_0*exp(-2*r^2/R_0^2) where R_0 is the 1/exp(2) radius
%  3: assumes all photons are initiated with x=0.

function [PhotonsX] = initialize(num_photons,beam_type,R_0)


if(beam_type == 1)
    PhotonsX = R_0*sqrt(rand(num_photons,1));
end
    
if(beam_type == 2)
    PhotonsX = R_0/sqrt(2)*sqrt(-log(rand(num_photons,1)));
end

if(beam_type == 3)
    PhotonsX = zeros(num_photons,1);
end