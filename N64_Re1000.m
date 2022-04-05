%% PROJECT 3
% This File runs the solver for the specified Reynolds number (Re=1000) 
% and Grid size (N=64) for the provided dimensions and Uwall. 

clc; clear; close all;
format long; format compact;

N=64; L=1; U_wall=1;
Re=5000; 
w=1.5;                % SOR Factor
B=3;              % Stability Factor

for i=1:length(Re)
% Using the Solver to compute U and V values as well as Residuals and
% number of iterations
[U,V,R,n]=Solver_GS(N,L,U_wall,Re(i),B,w);

%% Convergence plot
figure(1)
semilogy(100:100:n,R,'LineWidth',1.3)
grid on
xlabel('No. of iterations')
ylabel('Residual |R_L|')
title(['Plot of iterations vs Residual for Re = ',num2str(Re),' and grid N = ',num2str(N)])

%% Streamfunction
h=L/N;
psi=zeros(N+1,N+1);
for ii=1:N+1
    for j=1:N+1
        if ii==1&&j==1
            psi(ii,j)=0;
        elseif ii==1&&j~=1
            psi(ii,j)=psi(ii,j-1)-V(ii+1,j)*h;
        else
            psi(ii,j)=psi(ii-1,j)+U(ii,j+1)*h;
        end
    end
end
figure()
lvls=[-1e-10 -1e-7 -1e-5 -1e-4 -0.01 -0.03 -0.05 -0.07 -0.09 -0.1 -0.11 -0.115 -0.1175 1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
[C,h]=contourf(psi,lvls,'LineWidth',1);
colormap(white)
title(['Contour plot for Re = ',num2str(Re),' and Grid size N = ',num2str(N)])

%% Velocity Graph

ghia_v=[0 0.27485 0.29012 0.30353 0.32627 0.37095 0.33075 0.32235 0.02526 -0.31966 -0.42665 -0.51550 -0.39188 -0.33714 -0.27669 -0.21388 0];
ghia_u=[0 -0.18109 -0.20196 -0.22220 -0.29730 -0.38289 -0.27805 -0.10648 -0.06080 0.05702 0.18719 0.33304 0.46604 0.51117 0.57492 0.65928 1];

yg=linspace(0,1,129);
guu=[yg(1) yg(8) yg(9) yg(10) yg(14) yg(23) yg(37) yg(59) yg(65) yg(80) yg(95) yg(110) yg(123) yg(124) yg(125) yg(126) yg(129)];
guv=[yg(1) yg(9) yg(10) yg(11) yg(13) yg(21) yg(30) yg(31) yg(65) yg(104) yg(111) yg(117) yg(122) yg(123) yg(124) yg(125) yg(129)];
y=linspace(0,1,N);
figure()
plot(y,U((2:end-1),ceil(N/2)+1),'LineWidth',1.3)
hold on
scatter(guu,ghia_u,'o','filled')
grid on
xlabel('y','FontSize',15)
ylabel('U (x-velocity)')
title(['Plot of U vs y for Re = ',num2str(Re),' for grid N = ',num2str(N)])
legend('Solver','Ghia et al')
figure()
plot(y,(V(ceil(N/2)+1,(2:end-1))),'LineWidth',1.3)
hold on
scatter(guv,ghia_v,'o','filled')
grid on
xlabel('x','FontSize',15)
ylabel('V (y-velocity)')
title(['Plot of V vs x for Re = ',num2str(Re),' for grid N = ',num2str(N)])
legend('Solver','Ghia et al')
end