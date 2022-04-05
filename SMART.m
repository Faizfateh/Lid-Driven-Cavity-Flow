%% SMART Function
% Calculates the transported quantity phi_1/2 used for flux calculation
%INPUTS:
% U1 ,U2, U3: upwinded or downwinded velocity values depending upon the
%             sign of q
% OUTPUT: 
% phi_half: Transported quantity
function phi_half=SMART(U1,U2,U3)
    phi_hat=(U2-U1)/(U3-U1);
    if U3==U1
        phi_hat_half=U2;
    else
        if phi_hat<=0 || phi_hat>=1
            phi_hat_half=phi_hat;
        elseif phi_hat>0 && phi_hat<=(1/6)
            phi_hat_half=3*phi_hat;
        elseif phi_hat>(1/6) && phi_hat<=(5/6)
            phi_hat_half=(3/8)*(2*phi_hat+1);
        elseif phi_hat>(5/6) && phi_hat<1
            phi_hat_half=1;
        end
    end
    phi_half=phi_hat_half*(U3-U1)+U1;
end