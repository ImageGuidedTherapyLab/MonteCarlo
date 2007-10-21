%This function is meant to take a matrix and write it into an ascii 
%text file with the first line giving the 4-d dimension sizes of 
%the matrix


function ascii_write(A,B, filebase,delta_r,delta_z,...
                         g,mu_a,mu_s,dist_r,dist_z);

global BEAM_LOC 

%  A is the matrix being written
%  name is a base string of the file name to be written
%  delta_r  = spacing along the r-axis [units can be anything]
%  delta_z  = spacing along the z-axis [units can be anything]
%  g,mu_a,mu_s,dist_r,dist_r used to tag file w/ parameters used to make it.

fid = fopen( sprintf('%s.asc',filebase)  ,'w');

[m,n,p,q] = size(A);


% .asc files are written row major first index varies the fastest

fprintf(fid, '%d %d %d %d %f %f %f\n', m,n,p,q,...
              delta_z*.01,delta_r*.01,BEAM_LOC*.01);
% axial direction varies over rows
% radial direction varies over columns
for j = 1:n
    for i = 1:m
        fprintf(fid, '%f ', A(i,j));
    end
end
% remember params used for this file
fprintf(fid,'\n# g = %f, mu_a = %f, mu_s = %f, dist_r = %f, dist_z = %f \n',...
                                                   g,mu_a,mu_s,dist_r,dist_z);

fclose(fid);

% write avs field
fid = fopen( sprintf('%s.fld',filebase) ,'w');
% write out avs file header
fprintf(fid,'# AVS field file \n');
% remember params used for this file
fprintf(fid,'# g = %f, mu_a = %f, mu_s = %f, dist_r = %f, dist_z = %f \n',...
                                                   g,mu_a,mu_s,dist_r,dist_z);
fprintf(fid,'ndim=2 # number of dimensions in the field \n');
fprintf(fid,'dim1= %d  # dimension of r-axis \n',size(A,2));
fprintf(fid,'dim2= %d  # dimension of z-axis \n',size(A,1));
fprintf(fid,'nspace=2        # number of physical coordinates per point\n');
fprintf(fid,'veclen=2        # number of components at each point\n');
fprintf(fid,'label =monte_carlo  # data labels \n');
fprintf(fid,'label =isotropic    # data labels \n');
fprintf(fid,'unit = Watt/m_cubed  # units \n');
fprintf(fid,'unit = Watt/m_cubed  # units \n');
fprintf(fid,'min_val = %f %f # array bounds \n',min(min(A)),min(min(B)));
fprintf(fid,'max_val = %f %f # array bounds \n',max(max(A)),max(max(B)));
fprintf(fid,'data=float     # native format of linux\n');
fprintf(fid,'field=rectilinear  # field type(uniform,rectilinear,irregular)\n');
fprintf(fid,'\f\f');
% buffer the data
buffer = zeros(2*size(A,1)*size(A,2),1);
icnt = 1;
% TRANSPOSING THE MATRIX!!!!!!!!!!!!!!!!!!!!!!
% .fld files are written column major second index varies the fastest
for i = 1:size(A,1)
   for j = 1:size(A,2)
       buffer(icnt) = A(i,j); icnt = icnt + 1;
       buffer(icnt) = B(i,j); icnt = icnt + 1;
   end
end
fwrite(fid,buffer,'single');
% create rectilinear coordinate buffer 
% write on the coordinates of the variable that varies the fastest first
bounds = zeros(size(A,1)+size(A,2),1);
for i=0:size(A,2)-1
    bounds(i+1) = i*delta_r * .01; % meters
end
for i=0:size(A,1)-1
    bounds(size(A,2)+i+1) = i*delta_z * .01; % meters
end
fwrite(fid,bounds,'single');

fclose(fid);
