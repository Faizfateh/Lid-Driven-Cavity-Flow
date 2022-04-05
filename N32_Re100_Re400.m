%% PROJECT 3
% This File runs the solver for the specified Reynolds number (Re=100 & 400) 
% and Grid size (N=32) for the provided dimensions and Uwall. 

clc; clear; close all;
format long; format compact;

N=32; L=1; U_wall=1;
Re=[100,400]; 
w=1;                % SOR Factor
B=0.5;              % Stability Factor

for i=1:length(Re)
% Using the Solver to compute U and V values as well as Residuals and
% number of iterations
[U,V,R,n]=Solver_GS(N,L,U_wall,Re(i),B,w);

%% Convergence plot
figure()
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
if Re(i)==100
    ghia_v=[0 0.09233 0.10091 0.10890 0.12317 0.16077 0.17507 0.17527 0.05454 -0.24533 -0.22445 -0.16914 -0.10313 -0.08864 -0.07391 -0.05906 0];
    ghia_u=[0 -0.03717 -0.04192 -0.04775 -0.06434 -0.10150 -0.15662 -0.21090 -0.20581 -0.13641 0.00332 0.23151 0.68717 0.73722 0.78871 0.84123 1];
elseif Re(i)==400
    ghia_v=[0 0.18360 0.19713 0.20920 0.22965 0.28124 0.30203 0.30174 0.05186 -0.38598 -0.44993 -0.23827 -0.22847 -0.19254 -0.15663 -0.12146 0];
    ghia_u=[0 -0.08186 -0.09266 -0.10338 -0.14612 -0.24299 -0.32726 -0.17119 -0.11477 0.02135 0.16256 0.29093 0.55892 0.61756 0.68439 0.75837 1];
end
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
title(['Plot of U vs y for Re = ',num2str(Re(i)),' for grid N = ',num2str(N)])
legend('Solver','Ghia et al','Location','northwest')
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