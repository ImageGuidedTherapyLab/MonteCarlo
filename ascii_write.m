%This function is meant to take a matrix and write it into an ascii 
%text file with the first line giving the 4-d dimension sizes of 
%the matrix


function ascii_write(A, name)

%A is the matrix being written
%name is a string that is the name of the file to be written

fid = fopen(name,'w');

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