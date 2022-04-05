%% Flux Calculation
% Calculates the Flux using the SMART function
% INPUTS:
% U: Velocity values on the vertical lines of the grid
% V: Velocity values on the horizontal lines of the grid
% N: Grid Size
% h: distance between each grid point, L/N
% nu: viscocity
% 
% OUTPUTS:
% F: x-momentum flux on cell centers
% G: y-momentum flux on cell centers
% Hx: x-momentum flux on the grid nodes
% Hy: y-momentum flux on the grid nodes
function[F,G,Hx,Hy]=Flux_calc(U,V,N,h,nu)
F=zeros(N,N); G=zeros(N,N);
Hx=zeros(N+1,N+1); Hy=zeros(N+1,N+1);
%% F
for i=1:N
    for j=1:N
        q=0.5*(U(i+1,j+2)+U(i+1,j+1));
        if q>0
            phi_half=SMART(U(i+1,j),U(i+1,j+1),U(i+1,j+2));
        else
            phi_half=SMART(U(i+1,j+3),U(i+1,j+2),U(i+1,j+1));
        end
        F(i,j)=(q*phi_half)-(nu*(U(i+1,j+2)-U(i+1,j+1))/h);
    end
end

%% G
for i=1:N
    for j=1:N
        q=0.5*(V(i+2,j+1)+V(i+1,j+1));
        if q>0
            phi_half=SMART(V(i,j+1),V(i+1,j+1),V(i+2,j+1));
        else
            phi_half=SMART(V(i+3,j+1),V(i+2,j+1),V(i+1,j+1));
        end
        G(i,j)=(q*phi_half)-(nu*(V(i+2,j+1)-V(i+1,j+1))/h);
    end
end

%% Hx
% Filling the interior points
for i=2:N
    for j=2:N
        q=0.5*(V(i+1,j)+V(i+1,j+1));
        if q>0
            phi_half=SMART(U(i-1,j+1),U(i,j+1),U(i+1,j+1));
        else
            phi_half=SMART(U(i+2,j+1),U(i+1,j+1),U(i,j+1));
        end
        Hx(i,j)=(q*phi_half)-(nu*(U(i+1,j+1)-U(i,j+1))/h);
    end
end

% Filling in the top and bottom
for i=1:N+1
    Hx(1,i)=-(nu*(U(2,i+1)-U(1,i+1))/h);
    Hx(N+1,i)=-(nu*(U(N+2,i+1)-U(N+1,i+1))/h);
end

% Filling the boundary
for i=2:N
    Hx(i,1)=-(nu*(U(i+1,2)-U(i,2))/h);
    Hx(i,N+1)=-(nu*(U(i+1,N+2)-U(i,N+2))/h);
end

%% Hy
% Filling the interior points
for i=2:N
    for j=2:N
        q=0.5*(U(i,j+1)+U(i+1,j+1));
        if q>0
            phi_half=SMART(V(i+1,j-1),V(i+1,j),V(i+1,j+1));
        else
            phi_half=SMART(V(i+1,j+2),V(i+1,j+1),V(i+1,j));
        end
        Hy(i,j)=(q*phi_half)-(nu*(V(i+1,j+1)-V(i+1,j))/h);
    end
end

% Filling in the top and bottom
for i=1:N+1
    Hy(1,i)=-(nu*(V(2,i+1)-V(2,i))/h);
    Hy(N+1,i)=-(nu*(V(N+2,i+1)-V(N+2,i))/h);
end

% Filling in the boundary
for i=2:N
    Hy(i,1)=-(nu*(V(i+1,2)-V(i+1,1))/h);
    Hy(i,N+1)=-(nu*(V(i+1,N+2)-V(i+1,N+1))/h);
end