%% Velocity Calculations
% This function enforces the boundary condition to calculate the ghost
% cells of u and v velocity.
% INPUT: 
% N: Grid Size
% U: Velocity values on the vertical lines of the grid
% V: Velocity values on the horizontal lines of the grid
% U_wall: wall Velocity
%
% OUTOUTS:
% U_vel: Velocity values on the vertical lines of the grid with ghost cells
%        calculations
% V_vel: Velocity values on the horizontal lines of the grid with ghost cells
%        calculations
function [U_vel,V_vel]=velocity(N,U,V,U_wall)
    %% U
    U_vel=U;
    % Filling boundary
    for i=2:N+1
        U_vel(i,1)=U(i,3);
        U_vel(i,N+3)=U(i,N+1);
    end
    % Filling bottom and top
    for j=1:N+3
        U_vel(1,j)=-U(2,j);
        U_vel(N+2,j)=2*U_wall-U(N+1,j); %
    end

    %% V
    V_vel=V;
    % Filling boundary
    for i=1:N+2
        V_vel(i,1)=-V(i,2);
        V_vel(i,N+2)=-V(i,N+1);
    end
    % Filling bottom and top
    for i=1:N+2
        V_vel(1,i)=V(3,i);
        V_vel(N+3,i)=V(N+1,i);
    end
end
