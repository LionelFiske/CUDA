clear
fid = fopen('StokesU.out', 'r');

N=50;
%U=zeros(N-1,N);
%     for i=(0:N*(N-1)-1)  
%         
%         U(floor(i/(N))+1,mod(i,N)+1)=fread(fid,1,'*double');
%          
%     end
    for i=(1:N*(N-1))  
        
        U(i)=fread(fid,1,'*double');
         
    end
    
    U=(reshape(U,N,N-1)');



fclose(fid)


fid = fopen('StokesV.out', 'r');

% V=zeros(N-1,N);
%     for i=(0:N*(N-1) -1  ) 
%         V(floor(i/(N))+1,mod(i,N)+1)=fread(fid,1,'*double');
%     end
%     
     for i=(1:N*(N-1))  
        
        V(i)=fread(fid,1,'*double');
         
    end
    
    V=(reshape(V,N-1,N));


    
   
fclose(fid)

% 
fid = fopen('StokesP.out', 'r');
 
%P=zeros(N-1,N-1);
%    for i=(0:(N-1)*(N-1) -1 )
%        P(floor(i/(N-1))+1,mod(i,N-1)+1)=fread(fid,1,'*double');
%    end

%     
     for i=(1:(N-1)*(N-1))   
        
        P(i)=fread(fid,1,'*double');
         
    end
    
    P=(reshape(P,N-1,N-1)');




fclose(fid)


figure(2) 

subplot(3,1,1)
imagesc(U)
title('u')
colorbar  


subplot(3,1,2)
imagesc(V)
title('v')
colorbar  


subplot(3,1,3)
imagesc(P)
title('p ')
colorbar  


