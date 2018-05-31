function [u,v,p] = sorstokes(N,mu,P,omega,tol)
% sorstokes(N,mu,omega,P,tol)
% Use SOR to solve Stokes flow in a channel
% N = number of samples
% mu = viscosity
% P = pressure drop
% omega = relaxation parameter 0<omega<2
% tol = max residual tolerance

global M;
M = N;

figure(1)
fx = 0;
fy = 0;
x = linspace(0,1,N);
y = x;
[X,Y] = meshgrid(x,y);
X = X'; Y = Y';
dx = x(2)-x(1);
dy = dx;
u = zeros(N,N-1);
v = zeros(N-1,N);
p = zeros(N+1,N-1);
uresid = zeros(size(u));
vresid = zeros(size(v));
presid = zeros(size(p));

uresid(1,1) = 1000;
iter = 0;
while max(abs([uresid(:);vresid(:);presid(:)])) > tol
    % solve u equation

    i = 1;
 
    j = 1;
    uresid(i,j) = dy/dx*(-u(i,j)+u(i+1,j))+dx/dy*(-3*u(i,j)+u(i,j+1)) ...
        -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
    u(i,j) = u(i,j)+omega*uresid(i,j);

    for j=2:N-2
        uresid(i,j) = dy/dx*(-u(i,j)+u(i+1,j))+dx/dy*(u(i,j-1)-2*u(i,j)+u(i,j+1)) ...
            -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
        u(i,j) = u(i,j)+omega*uresid(i,j);
    end

    j = N-1;
    uresid(i,j) = dy/dx*(-u(i,j)+u(i+1,j))+dx/dy*(-3*u(i,j)+u(i,j-1)) ...
        -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
    u(i,j) = u(i,j)+omega*uresid(i,j);

    for i=2:N-2
        
        j = 1;
        uresid(i,j) = dy/dx*(u(i-1,j)-2*u(i,j)+u(i+1,j))+dx/dy*(-3*u(i,j)+u(i,j+1)) ...
            -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
        u(i,j) = u(i,j)+omega*uresid(i,j);

        for j=2:N-2
            uresid(i,j) = dy/dx*(u(i-1,j)-2*u(i,j)+u(i+1,j))+dx/dy*(u(i,j-1)-2*u(i,j)+u(i,j+1)) ...
                -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
            u(i,j) = u(i,j)+omega*uresid(i,j);            
        end
        
        j = N-1;
        uresid(i,j) = dy/dx*(u(i-1,j)-2*u(i,j)+u(i+1,j))+dx/dy*(-3*u(i,j)+u(i,j-1)) ...
            -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
        u(i,j) = u(i,j)+omega*uresid(i,j);
   end
    
    i = N-1;
 
    j = 1;
    uresid(i,j) = dy/dx*(u(i-1,j)-2*u(i,j)+u(i+1,j))+dx/dy*(-3*u(i,j)+u(i,j+1)) ...
        -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
    u(i,j) = u(i,j)+omega*uresid(i,j);

    for j=2:N-2
        uresid(i,j) = dy/dx*(u(i-1,j)-2*u(i,j)+u(i+1,j))+dx/dy*(u(i,j-1)-2*u(i,j)+u(i,j+1)) ...
            -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
        u(i,j) = u(i,j)+omega*uresid(i,j);
    end

    j = N-1;
    uresid(i,j) = dy/dx*(u(i-1,j)-2*u(i,j)+u(i+1,j))+dx/dy*(u(i,j-1)-3*u(i,j)) ...
        -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
    u(i,j) = u(i,j)+omega*uresid(i,j);

    i = N;
 
    j = 1;
    uresid(i,j) = dy/dx*(u(i-1,j)-u(i,j))+dx/dy*(-3*u(i,j)+u(i,j+1)) ...
        -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
    u(i,j) = u(i,j)+omega*uresid(i,j);

    for j=2:N-2
        uresid(i,j) = dy/dx*(u(i-1,j)-u(i,j))+dx/dy*(u(i,j-1)-2*u(i,j)+u(i,j+1)) ...
            -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
        u(i,j) = u(i,j)+omega*uresid(i,j);
    end

    j = N-1;
    uresid(i,j) = dy/dx*(u(i-1,j)-u(i,j))+dx/dy*(u(i,j-1)-3*u(i,j)) ...
        -dy/mu*(p(i+1,j)-p(i,j))+fx*dx*dy/mu;
    u(i,j) = u(i,j)+omega*uresid(i,j);
    
    % solve v equation
    
    i = 1;
 
    j = 1;
    vresid(i,j) = v(i,j);
    v(i,j) = 0;

    for j=2:N-2
        vresid(i,j) = dy/dx*(-v(i,j)+v(i+1,j))+dx/dy*(v(i,j-1)-2*v(i,j)+v(i,j+1)) ...
            -dx/mu*(p(i+1,j)-p(i+1,j-1))+fy*dx*dy/mu;
        v(i,j) = v(i,j)+omega*vresid(i,j);
    end

    j = N-1;
    vresid(i,j) = dy/dx*(-v(i,j)+v(i+1,j))+dx/dy*(v(i,j-1)-2*v(i,j)+v(i,j+1)) ...
        -dx/mu*(p(i+1,j)-p(i+1,j-1))+fy*dx*dy/mu;
    v(i,j) = v(i,j)+omega*vresid(i,j);

    j = N;

    v(i,j) = 0;  

    for i=2:N-2
        
        j = 1;
        vresid(i,j) = v(i,j);
        v(i,j) = 0;
        
        for j=2:N-2
            vresid(i,j) = dy/dx*(v(i-1,j)-2*v(i,j)+v(i+1,j))+dx/dy*(v(i,j-1)-2*v(i,j)+v(i,j+1)) ...
                -dx/mu*(p(i+1,j)-p(i+1,j-1))+fy*dx*dy/mu;
            v(i,j) = v(i,j)+omega*vresid(i,j);
        end
        
        j = N-1;
        vresid(i,j) = dy/dx*(v(i-1,j)-2*v(i,j)+v(i+1,j))+dx/dy*(v(i,j-1)-2*v(i,j)+v(i,j+1)) ...
            -dx/mu*(p(i+1,j)-p(i+1,j-1))+fy*dx*dy/mu;
        v(i,j) = v(i,j)+omega*vresid(i,j);
        
        j = N;
        
        v(i,j) = 0;        

    end
    
    i = N-1;
 
    j = 1;
    vresid(i,j) = 0;
    v(i,j) = 0;

    for j=2:N-2
        vresid(i,j) = dy/dx*(v(i-1,j)-v(i,j))+dx/dy*(v(i,j-1)-2*v(i,j)+v(i,j+1)) ...
            -dx/mu*(p(i+1,j)-p(i+1,j-1))+fy*dx*dy/mu;
        v(i,j) = v(i,j)+omega*vresid(i,j);
    end

    j = N-1;
    vresid(i,j) = dy/dx*(v(i-1,j)-v(i,j))+dx/dy*(v(i,j-1)-2*v(i,j)+v(i,j+1)) ...
        -dx/mu*(p(i+1,j)-p(i+1,j-1))+fy*dx*dy/mu;
    v(i,j) = v(i,j)+omega*vresid(i,j);

    j = N;

    v(i,j) = 0;        

    % solve p equation
    
    i = 1;
 
    j = 1;
    presid(i,j) = 2*P-p(i+1,j)-p(i,j);
    p(i,j) = 2*P-p(i+1,j);

    for j=2:N-2
        presid(i,j) = 2*P-p(i+1,j)-p(i,j);
        p(i,j) = 2*P-p(i+1,j);
    end

    j = N-1;
    presid(i,j) = 2*P-p(i+1,j)-p(i,j);
    p(i,j) = 2*P-p(i+1,j);

    for i=2:N-2
        
        j = 1;
        presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
        p(i,j) = p(i,j)+omega*presid(i,j);
        
        for j=2:N-2
            presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
            p(i,j) = p(i,j)+omega*presid(i,j);
        end
        
        j = N-1;
        presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
        p(i,j) = p(i,j)+omega*presid(i,j);
        
    end
    
    i = N-1;
 
    j = 1;
    presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
    p(i,j) = p(i,j)+omega*presid(i,j);

    for j=2:N-2
        presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
        p(i,j) = p(i,j)+omega*presid(i,j);
    end

    j = N-1;
    presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
    p(i,j) = p(i,j)+omega*presid(i,j);

    i = N;
 
    j = 1;
    presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
    p(i,j) = p(i,j)+omega*presid(i,j);

    for j=2:N-2
        presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
        p(i,j) = p(i,j)+omega*presid(i,j);
    end

    j = N-1;
    presid(i,j) = -(u(i,j)-u(i-1,j))-dx/dy*(v(i-1,j+1)-v(i-1,j));
    p(i,j) = p(i,j)+omega*presid(i,j);
    
    i = N+1;
    for j=1:N-1
        p(i,j) = -p(i-1,j);
    end
    
     if mod(iter,100) == 0
        max(abs([uresid(:);vresid(:);presid(:)]))        
        subplot(1,3,1)
        surf(X(:,1:N-1),Y(:,1:N-1)+dy/2,u);
        set(gca,'DataAspectRatio',[1,1,1e-1]);
        title('U');
        xlabel('x');
        ylabel('y');
        subplot(1,3,2)
        surf(X(1:N-1,:)+dx/2,Y(1:N-1,:),v);
        title('V');
        xlabel('x');
        ylabel('y');
        axis([0,1,0,1,-0.5,0.5]);
        subplot(1,3,3)
        surf([X(1,1:N-1)-dx/2;X(:,1:N-1)+dx/2],[Y(1,1:N-1)+dy/2;Y(:,1:N-1)+dy/2],p);
        title('P');
        xlabel('x');
        ylabel('y');
        axis equal
        drawnow
     end
    iter = iter+1;
end

iter

end