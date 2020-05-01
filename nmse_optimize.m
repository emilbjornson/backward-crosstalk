% This Matlab script is the function "nmse_optimize" used to generate
% some figures and finds the optimal input reference power minimizing
% the maximum of analytical NMSEs in the paper:
%
% Peter Händel, Özlem Tugfe Demir, Emil Björnson, and Daniel Rönnow,
% "Impact of Backward Crosstalk in 2 × 2 MIMO Transmitters on NMSE and Spectral Efficiency,"
% IEEE Transactions on Communications, To appear.
%
% Download article: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9072380
%
% This is version 1.0 (Last edited: 2020-04-28)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

function sigmaStar=nmse_optimize(xi,betaa,kappa1,kappa2,...
    sigmaw2,rho1,rho2,gamma1,gamma2)

t11=gamma1^2+2*gamma1^2*gamma2*betaa*real(kappa2'*xi)...
    +gamma1^2*gamma2^2*abs(kappa2)^2*betaa^2;

t22=gamma2^2*betaa^2+2*gamma1*gamma2^2*betaa*real(kappa1*xi)...
    +gamma1^2*gamma2^2*abs(kappa1)^2;

% Prepare array for the candidate solutions for Case 1,2, and 3.
candidates=[];
% Case 1
kokler=roots([12*rho1^2*t11^3,  ...
    4*gamma1^2*gamma2*t11*rho1*...
    (gamma2*betaa^2*abs(kappa2)^2+betaa*real(kappa2*xi')),...
    0, -sigmaw2]);
kokler=kokler(real(kokler)>0);
[~,indexx]=min(abs(imag(kokler)));
candidates=[candidates real(kokler(indexx))];

% Case 2
kokler=roots([12*rho2^2*t22^3,  ...
    4*gamma1*gamma2^2*t22*rho2*...
    (gamma1*abs(kappa1)^2+betaa*real(kappa1*xi)),...
    0, -sigmaw2]);
kokler=kokler(real(kokler)>0);
[~,indexx]=min(abs(imag(kokler)));
candidates=[candidates real(kokler(indexx))];

% Case 3

kokler=roots([ 6*(rho1^2*t11^3*gamma2^2*betaa^2-rho2^2*t22^3*gamma1^2)...
    /(gamma1^2*gamma2^2*betaa^2), ...
    4*gamma2*t11*rho1*(gamma2*betaa^2*abs(kappa2)^2+betaa*real(kappa2*xi'))-4*gamma1*t22*rho2/(betaa^2)*...
    (gamma1*abs(kappa1)^2+betaa*real(kappa1*xi)), ...
    betaa^2*gamma2^2*abs(kappa2)^2-gamma1^2*abs(kappa1)^2/(betaa^2),...
    (sigmaw2/(gamma1^2)-sigmaw2/(gamma2^2*betaa^2))]);
kokler=kokler(real(kokler)>0);
kokler=kokler(abs(imag(kokler))./abs(real(kokler))<10^(-3));
candidates=[candidates real(kokler.')];
minn=inf;
for sigmax2=candidates
    
    nmse1=6*rho1^2*t11^3/(gamma1^2)*sigmax2^2 ...
        +4*gamma2*t11*rho1*(gamma2*betaa^2*abs(kappa2)^2+betaa*real(kappa2*xi'))*sigmax2...
        +betaa^2*gamma2^2*abs(kappa2)^2+sigmaw2/(gamma1^2*sigmax2);
    
    nmse2=6*rho2^2*t22^3/(gamma2^2*betaa^2)*sigmax2^2 ...
        +4*gamma1*t22*rho2/(betaa^2)*...
        (gamma1*abs(kappa1)^2+betaa*real(kappa1*xi))*sigmax2...
        +gamma1^2*abs(kappa1)^2/(betaa^2)+sigmaw2/(gamma2^2*betaa^2*sigmax2);
    if max(nmse1,nmse2)<minn
        sigmaStar=sigmax2;
        minn=max(nmse1,nmse2);
    end
end
