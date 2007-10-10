function sensitivity4(k)


heating4 = zeros(512,512);

   for i = 1:k
        i
    
        [R,A,T,grid,fluence,HEATING] = scatter_simulation(100000,2,.38);
        heating4 = heating4 + HEATING;
    end
        heating4 = heating4./k;  

    ascii_write(heating4,'Original.asc')


a = [.44,.704,1.008,1.312,1.616,1.92,2.224,2.528,2.832,3.136,3.44,3.744,4.048,4.352,4.656,4.96,5.0];


for j = 1:17
   
    heating4 = zeros(512,512);
    for i = 1:k
        i
    
        [R,A,T,grid,fluence,HEATING] = scatter_simulation_mua(100000,2,.38,a(j));
        heating4 = heating4 + HEATING;
    end
        heating4 = heating4./k;  
        name = sprintf('Gauss_mua%d.asc',j);
        ascii_write(heating4,name);
end

for j = 1:17

    heating4 = zeros(512,512);
    for i = 1:k
        i
        [R,A,T,grid,fluence,HEATING] = scatter_simulation_mua(100000,1,.38,a(j));
        heating4 = heating4 + HEATING;
    end
        heating4 = heating4./k;
        name = sprintf('Flat_mua%d.asc',j);
        ascii_write(heating4,name);
end


g = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,.99];

for j = 1:11
   
    heating4 = zeros(512,512);
    for i = 1:k
        i
   
        [R,A,T,grid,fluence,HEATING] = scatter_simulation_g(100000,2,.38,g(j));
        heating4 = heating4 + HEATING;
    end
        heating4 = heating4./k;  
        name = sprintf('Gauss_g%d.asc',j);
        ascii_write(heating4,name);
end


for j = 1:11

    heating4 = zeros(512,512);
    for i = 1:k
        i

        [R,A,T,grid,fluence,HEATING] = scatter_simulation_g(100000,1,.38,g(j));
        heating4 = heating4 + HEATING;
    end
        heating4 = heating4./k;
        name = sprintf('Flat_g%d.asc',j);
        ascii_write(heating4,name);
end



s = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25];
      

for j = 1:25
   
    heating4 = zeros(512,512);
    for i = 1:k
        i
    
        [R,A,T,grid,fluence,HEATING] = scatter_simulation_mus(100000,2,.38,s(j));
        heating4 = heating4 + HEATING;
    end
        heating4 = heating4./k;  
        name = sprintf('Gauss_mus%d.asc',j);
        ascii_write(heating4,name);
end


for j = 1:25

    heating4 = zeros(512,512);
    for i = 1:k
        i
        [R,A,T,grid,fluence,HEATING] = scatter_simulation_mus(100000,1,.38,s(j));
       heating4 = heating4 + HEATING;
    end
       heating4 = heating4./k;
       name = sprintf('Flat_mus%d.asc',j);
       ascii_write(heating4,name);
end

