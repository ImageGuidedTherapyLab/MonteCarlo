%This function is meant to take a matrix and write it into an ascii 
%text file with the first line giving the 4-d dimension sizes of 
%the matrix


function ascii_write(A,B, filebase,delta_r,delta_z)

%  A is the matrix being written
%  name is a base string of the file name to be written
%  delta_r  = spacing along the r-axis [units can be anything]
%  delta_z  = spacing along the z-axis [units can be anything]

fid = fopen( sprintf('%s.asc',filebase)  ,'w');

[m,n,p,q] = size(A);


fprintf(fid, '%d %d %d %d \n', m,n,p,q);

if(p == 1 && q ==1)
  for j = 1:n
      for i = 1:m
          fprintf(fid, '%f ', A(i,j));
      end
      fprintf(fid, '\n');
  end
end

fclose(fid);

fid = fopen( sprintf('%s.fld',filebase) ,'w')
% create rectilinear coordinate buffer
bounds = zeros(size(A,1)+size(A,2),1);
for i=0:size(A,1)-1
    bounds(i+1) = i*delta_r;
end
for i=0:size(A,2)-1
    bounds(size(A,1)+i+1) = i*delta_z;
end
% write out avs file header
fprintf(fid,'# AVS field file \n');
fprintf(fid,'ndim=2 # number of dimensions in the field \n');
fprintf(fid,'dim1= %d  # dimension of axis 1 \n',size(A,1));
fprintf(fid,'dim2= %d  # dimension of axis 2 \n',size(A,2));
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
for i = 1:size(A,1)
   for j = 1:size(A,2)
       buffer(icnt) = A(i,j); icnt = icnt + 1;
       buffer(icnt) = B(i,j); icnt = icnt + 1;
   end
end
fwrite(fid,buffer,'single');
fwrite(fid,bounds,'single');

fclose(fid);
