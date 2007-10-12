%This function is meant to take a matrix and pixel dimensions and 
% write it into an avs field file 



function createavsfield(A,delta_r,delta_z)


fid = fopen('mc.fld','w')
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
fprintf(fid,'veclen=1        # number of components at each point\n');
fprintf(fid,'data=float     # native format of linux\n');
fprintf(fid,'field=rectilinear  # field type(uniform,rectilinear,irregular)\n');
fprintf(fid,'\f\f');
fwrite(fid,A','single');
fwrite(fid,bounds,'single');

fclose(fid);
