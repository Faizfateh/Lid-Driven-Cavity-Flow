%% PROJECT 3
% This File runs the solver for the specified Reynolds number (Re=3200 & 5000) 
% and Grid size (N=128) for the provided dimensions and Uwall. 

clc; clear; close all;
format long; format compact;

N=128; L=1; U_wall=1;
Re=[3200,5000]; 
w=1.6;                % SOR Factor
B=5;              % Stability Factor

for i=1:length(Re)
% Using the Solver to compute U and V values as well as Residuals and
% number of iterations
[U,V,R,n]=Solver_RBGS(N,L,U_wall,Re(i),B,w);

%% Convergence plot
figure(1)
semilogy(100:100:n,R,'LineWidth',1.3)
grid on
xlabel('No. of iterations')
ylabel('Residual |R_L|')
title(['Plot of iterations vs Residual for Re = ',num2str(Re(i)),' and grid N = ',num2str(N)])

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
title(['Contour plot for Re = ',num2str(Re(i)),' and Grid size N = ',num2str(N)])

%% Velocity Graph
if Re(i)==3200
    ghia_u=[0 -0.32407 -0.35344 -0.37827 -0.41933 -0.34323 -0.024427 -0.04272...
        0.07156 0.19791 0.34682 0.46101 0.46547 0.48296 0.53236 1];
    ghia_v=[0 0.39560 0.40917 0.41906 0.42768 0.37119 0.29030 0.28188 0.00999 -0.31184...
        -0.37401 -0.44307 -0.54053 -0.52357 -0.47425 -0.39017 0];
elseif Re(i)
    ghia_u = [0.0, -0.41165, -0.42901, -0.43643, -0.40435, -0.33050, -0.22855,...
        -0.07404, 0.03039, 0.08183, 0.20087, 0.33556, 0.46036, 0.45992, 0.46120, 0.48223,1];
    ghia_v =[0.0, 0.42447, 0.43329, 0.43648, 0.42951, 0.35368, 0.28066, 0.27280, 0.00945...
         , -0.30018, -0.36214, -0.41442, -0.52876, -0.55408, -0.55069, -0.49774, 0];
end
yg=linspace(0,1,129);
if Re(i)==3200
    guu=[yg(1) yg(8) yg(9) yg(10) yg(14) yg(23) yg(37) yg(65) yg(80) yg(95) yg(110) yg(123) yg(124) yg(125) yg(126) yg(129)];
else
    guu=[yg(1) yg(8) yg(9) yg(10) yg(14) yg(23) yg(37) yg(59) yg(65) yg(80) yg(95) yg(110) yg(123) yg(124) yg(125) yg(126) yg(129)];
end
guv=[yg(1) yg(9) yg(10) yg(11) yg(13) yg(21) yg(30) yg(31) yg(65) yg(104) yg(111) yg(117) yg(122) yg(123) yg(124) yg(125) yg(129)];
y=linspace(0,1,N);
figure()
plot(y,U((2:end-1),ceil(N/2)+1),'LineWidth',1.3)
hold on
scatter(guu,ghia_u,'o','filled')
grid on
xlabel('y','FontSize',15)
ylabel('U (x-velocity)')
title(['Plot of U vs y for Re = ',num2str(Re(i)),' for grid N = ',num2str(N)])
legend('Solver','Ghia et al')
figure()
plot(y,(V(ceil(N/2)+1,(2:end-1))),'LineWidth',1.3)
hold on
scatter(guv,ghia_v,'o','filled')
grid on
xlabel('x','FontSize',15)
ylabel('V (y-velocity)')
title(['Plot of V vs x for Re = ',num2str(Re(i)),' for grid N = ',num2str(N)])
legend('Solver','Ghia et al')
end