% This function calculates the Velocity and the Residual for Lid driven
% cavity problem for specified Reynolds number and grid size using the
% Gauss Seidel iterative method.
%
% INPUTS:
% N: Grid size
% L: length of the edge od the cavity
% U_wall: wall Velocity
% Re: Reynolds number
% Beta: Stability factor
% w: omega for succesive over relaxation in Gauss Seidel
%
% OUTPUTS:
% U: Velocity values on the vertical lines of the grid
% V: Velocity values on the Horizontal lines of the grid
% R: Residuals
% iter: number of iterations

function [U,V,R,iter]=Solver_GS(N,L,U_wall,Re,Beta,w)
    h=L/N;
    nu=U_wall*L/Re;
    dt=Beta*min((h^2)/(4*nu),(4*nu)/(U_wall^2));

    P=zeros(N,N);  U=zeros(N+2,N+3);   V=zeros(N+3,N+2);   R=zeros();
    iter=1; time=0; n=1; Residual=1;
    
    while Residual>1e-5
        %% Velocity initialization
        [U,V]=velocity(N,U,V,U_wall);
        
        %% Flux Calculation
        [F,G,Hx,Hy]=Flux_calc(U,V,N,h,nu);

        %% Updating Velocity
        [C,D]=size(U);
        for i=1:C-2
            for j=1:D-4
                Fx=(F(i,j+1)-F(i,j))/h;
                Hxy=(Hx(i+1,j+1)-Hx(i,j+1))/h;
                U(i+1,j+2)=U(i+1,j+2)-dt*(Fx+Hxy);

                Gy=(G(j+1,i)-G(j,i))/h;
                Hyx=(Hy(j+1,i+1)-Hy(j+1,i))/h;
                V(j+2,i+1)=V(j+2,i+1)-dt*(Gy+Hyx);
            end
        end
        
        %% Pressure Poisson
        for e=1:30
        for i=1:N
            for j=1:N
                ux=(U(i+1,j+2)-U(i+1,j+1));
                vy=(V(i+2,j+1)-V(i+1,j+1));
                if i==1
                    if j==1
                        P(i,j)=(1-w)*(P(i,j))+(w/2)*(P(i,j+1)+P(i+1,j)-(h/dt)*(ux+vy));
                    elseif j>1&&j<N
                        P(i,j)=(1-w)*(P(i,j))+(w/3)*(P(i,j-1)+P(i,j+1)+P(i+1,j)-(h/dt)*(ux+vy));
                    elseif j==N
                        P(i,j)=(1-w)*(P(i,j))+(w/2)*(P(i,j-1)+P(i+1,j)-(h/dt)*(ux+vy));
                    end
                elseif i>1 && i<N
                    if j==1
                        P(i,j)=(1-w)*(P(i,j))+(w/3)*(P(i,j+1)+P(i+1,j)+P(i-1,j)-(h/dt)*(ux+vy));
                    elseif j>1&&j<N
                        P(i,j)=(1-w)*(P(i,j))+(w/4)*(P(i,j+1)+P(i+1,j)+P(i-1,j)+P(i,j-1)-(h/dt)*(ux+vy));
                    elseif j==N
                        P(i,j)=(1-w)*(P(i,j))+(w/3)*(P(i,j-1)+P(i+1,j)+P(i-1,j)-(h/dt)*(ux+vy));
                    end
                elseif i==N
                    if j==1
                        P(i,j)=(1-w)*(P(i,j))+(w/2)*(P(i,j+1)+P(i-1,j)-(h/dt)*(ux+vy));
                    elseif j>1&&j<N
                        P(i,j)=(1-w)*(P(i,j))+(w/3)*(P(i,j-1)+P(i,j+1)+P(i-1,j)-(h/dt)*(ux+vy));
                    elseif j==N
                        P(i,j)=(1-w)*(P(i,j))+(w/2)*(P(i,j-1)+P(i-1,j)-(h/dt)*(ux+vy));
                    end
                end
            end
        end
        end

        %% Correcting Velocity field
        [C,D]=size(U);
        for i=1:C-2
            for j=1:D-4
                U(i+1,j+2)=U(i+1,j+2)-(dt/h)*(P(i,j+1)-P(i,j));

                V(j+2,i+1)=V(j+2,i+1)-(dt/h)*(P(j+1,i)-P(j,i));
            end
        end

        %% Residuals
        Rx=zeros();
        Ry=zeros();
        k=1;
        for i=1:N
            for j=1:N-1
                Rx(k)=abs(h*((F(i,j+1)+P(i,j+1))-(F(i,j)+P(i,j)))+h*(Hx(i+1,j+1)-Hx(i,j+1)));
                Ry(k)=abs(h*((G(j+1,i)+P(j+1,i))-(G(j,i)+P(j,i)))+h*(Hy(j+1,i+1)-Hy(j+1,i)));
                k=k+1;
            end
        end
        Residual=sum(Rx)+sum(Ry); 

        iter=iter+1;
        time=time+dt;
        if rem(iter,100)==0
            R(n)=Residual;
            disp(['Iterations: ',num2str(iter),' Time: ',num2str(time),' Residual: ',num2str(Residual)])
            n=n+1;
        end
    end
end