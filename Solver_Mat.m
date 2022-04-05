function [U,V,R,iter]=Solver_Mat(N,L,U_wall,Re,Beta)
    h=L/N;
    nu=U_wall*L/Re;
    dt=Beta*min((h^2)/(4*nu),(4*nu)/(U_wall^2));
    
    %% Velocity initialization
    U=zeros(N+2,N+3);
    V=zeros(N+3,N+2);
    [U,V]=velocity(N,U,V,U_wall);
    R=zeros();
    iter=1;
    time=0;
    n=1;
    Residual=1;
    while Residual>1e-5

        %% Flux Calculation
        [F,G,Hx,Hy]=Flux_calc(U,V,N,h,nu);

        %% Updating Velocity
        %U
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
        [U,V]=velocity(N,U,V,U_wall);
        
        %% Pressure Poisson
        rhs=zeros(N^2,1);
        k=1;
        for i=1:N
            for j=1:N
                ux=(U(i+1,j+2)-U(i+1,j+1));
                vy=(V(i+2,j+1)-V(i+1,j+1));
                rhs(k)=(1/dt)*(ux+vy);
                k=k+1;
            end
        end
        rhs(1)=0;
        A=laplace(N)/(h^2);
        B=A\rhs;
        P=(reshape(B,N,N))';

        %% Correcting Velocity field
        [C,D]=size(U);
        for i=1:C-2
            for j=1:D-4
                U(i+1,j+2)=U(i+1,j+2)-(dt/h)*(P(i,j+1)-P(i,j));

                V(j+2,i+1)=V(j+2,i+1)-(dt/h)*(P(j+1,i)-P(j,i));
            end
        end
        [U,V]=velocity(N,U,V,U_wall);

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
% clc; clear; close all;
% 
% N=32; L=1; U_wall=1; B=0.5;
% Re=100; h=L/N;
% nu=U_wall*L/Re;
% dt=B*min((h^2)/(4*nu),(4*nu)/(U_wall^2));
% 
% % function [U,V,R,n]=Solver()
% P=zeros(N,N);
% %% Velocity initialization
% U=zeros(N+2,N+3);
% V=zeros(N+3,N+2);
% [U,V]=velocity(N,U,V,U_wall);
% n=1;
% t=0;
% k=1;
% Residual=1;
% while Residual>1e-5
% 
%     %% Flux Calculation
%     [F,G,Hx,Hy]=Flux_calc(U,V,N,h,nu);
% 
%     %% Updating Velocity
%     %U
%     [C,D]=size(U);
%     for i=1:C-2
%         for j=1:D-4
%             Fx=(F(i,j+1)-F(i,j))/h;
%             Hxy=(Hx(i+1,j+1)-Hx(i,j+1))/h;
%             U(i+1,j+2)=U(i+1,j+2)-dt*(Fx+Hxy);
% 
%             Gy=(G(j+1,i)-G(j,i))/h;
%             Hyx=(Hy(j+1,i+1)-Hy(j+1,i))/h;
%             V(j+2,i+1)=V(j+2,i+1)-dt*(Gy+Hyx);
%         end
%     end
%     
%     [U,V]=velocity(N,U,V,U_wall);
%     %% Pressure Poisson
%     % RHS
%     rhs=zeros(N^2,1);
%     k=1;
%     for i=1:N
%         for j=1:N
%             ux=(U(i+1,j+2)-U(i+1,j+1));
%             vy=(V(i+2,j+1)-V(i+1,j+1));
% %             rhs(k)=(1/dt)*(ux+vy);
% %             k=k+1;
%             if i==1
%                 if j==1
%                     P(i,j)=(1/2)*(P(i,j+1)+P(i+1,j)-(h/dt)*(ux+vy));
%                 elseif j>1&&j<N
%                     P(i,j)=(1/3)*(P(i,j-1)+P(i,j+1)+P(i+1,j)-(h/dt)*(ux+vy));
%                 elseif j==N
%                     P(i,j)=(1/2)*(P(i,j-1)+P(i+1,j)-(h/dt)*(ux+vy));
%                 end
%             elseif i>1 && i<N
%                 if j==1
%                     P(i,j)=(1/3)*(P(i,j+1)+P(i+1,j)+P(i-1,j)-(h/dt)*(ux+vy));
%                 elseif j>1&&j<N
%                     P(i,j)=(1/4)*(P(i,j+1)+P(i+1,j)+P(i-1,j)+P(i,j-1)-(h/dt)*(ux+vy));
%                 elseif j==N
%                     P(i,j)=(1/3)*(P(i,j-1)+P(i+1,j)+P(i-1,j)-(h/dt)*(ux+vy));
%                 end
%             elseif i==N
%                 if j==1
%                     P(i,j)=(1/2)*(P(i,j+1)+P(i-1,j)-(h/dt)*(ux+vy));
%                 elseif j>1&&j<N
%                     P(i,j)=(1/3)*(P(i,j-1)+P(i,j+1)+P(i-1,j)-(h/dt)*(ux+vy));
%                 elseif j==N
%                     P(i,j)=(1/2)*(P(i,j-1)+P(i-1,j)-(h/dt)*(ux+vy));
%                 end
%             end
%         end
%     end
% %     rhs(1)=0;
% %     A=laplace(N)/(h^2);
% %     B=A\rhs;
% %     P=(reshape(B,N,N))';
% 
%     %% Correcting Velocity field
%     [C,D]=size(U);
%     for i=1:C-2
%         for j=1:D-4
%             U(i+1,j+2)=U(i+1,j+2)-(dt/h)*(P(i,j+1)-P(i,j));
% 
%             V(j+2,i+1)=V(j+2,i+1)-(dt/h)*(P(j+1,i)-P(j,i));
%         end
%     end
% 
%     %% Residuals
%     Rx=zeros();
%     Ry=zeros();
%     for i=1:N
%         for j=1:N-1
%             Rx(k)=abs(h*((F(i,j+1)+P(i,j+1))-(F(i,j)+P(i,j)))+h*(Hx(i+1,j+1)-Hx(i,j+1)));
%             Ry(k)=abs(h*((G(j+1,i)+P(j+1,i))-(G(j,i)+P(j,i)))+h*(Hy(j+1,i+1)-Hy(j+1,i)));
%             k=k+1;
%         end
%     end
%     Residual=sum(Rx)+sum(Ry); 
%     R(n)=Residual;
%     k=k+1;
%     [U,V]=velocity(N,U,V,U_wall);
%     n=n+1;
%     t=t+dt;
%     if rem(n,100)==0
%         disp(['Iterations: ',num2str(n),' Time: ',num2str(t),' Residual: ',num2str(Residual)])
%     end
% %     figure(1)
% %     contourf(U)
% %     colorbar
% end
% % figure()
% % contourf(flip(U(2:end-1,2:end-1)),20,'LineStyle','None')
% % colorbar
% % figure()
% % contourf(V)
% % colorbar