% This Matlab script generates Figure 5 in the paper:

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

% Amplification gain of the power amplifiers
gamma=10.^(30/20);
gamma1=gamma;
gamma2=gamma;


figure


% Correlation coefficient between x_1 and x_2
xi=0.7;

% Ratio of second input signal's amplitude to the first one
betaa=1.3;

% Crosstalk parameters \kappa_1 and \kappa_2

kappa1=10.^(-48/20);
kappa2=10.^(-52/20);

% Thermal noise variance
sigmaw2=db2pow(-40);

% Compression Parameters
rho1=-0.023;
rho2=-0.027;

% Construct the matrices in the paper
G=diag([rho1,rho2]);
L=diag([gamma1,gamma2]);
K=[0, gamma1*kappa2; gamma2*kappa1, 0];
Q=[gamma1, gamma1*gamma2*kappa2; gamma1*gamma2*kappa1, gamma2];

% Number of random input signal realizations for the simulation results
% (non-linear system solution)
NbrOfRealizations=10000;

% Prepare the array for the NMSE values obtained by non-linear system
% solution
MSE_sim=zeros(2,length(db2pow(-60:1:-20)));

% Prepare the values of u_1 and u_2 obtained by simulation and approximate
% analytical result in the paper
U_sim=zeros(2,NbrOfRealizations,length(db2pow(-60:1:-20)));
U_approx=zeros(2,NbrOfRealizations,length(db2pow(-60:1:-20)));

ppp=0;
options = optimset('Display','off');

% Go through different values of the input reference power P_x
for sigmax2=db2pow(-60:1:-20)
    ppp=ppp+1;
    
    % Input covariance matrix
    Cx=sigmax2*[1, betaa*xi; betaa*xi', betaa^2];
    
    % Generate white Gaussian unit-variance random variables
    W=sqrt(0.5)*(randn(2, NbrOfRealizations)+1i*randn(2, NbrOfRealizations));
    
    % Generate random input with covariance matrix Cx
    X=sqrtm(Cx)*W;
    
    % The analytical values for u_1 and u_2 in the paper
    U_approx(:,:,ppp)=Q*X;
    
    % Go through all the random input realizations and solve the non-linear
    % system to find the simulation values for u_1, u_2, r_1, r_2, y_1,
    % y_2, y_{01}, and y_{02}
    for nnn=1:NbrOfRealizations
        x=X(:,nnn);
        f = @(u) backward_crosstalk(u,x,L,K,G);
        [out,fval]=fsolve(f,Q*x,options);
        
        U_sim(:,nnn,ppp)=out;
        r=out+G*[abs(out(1))^2*out(1); abs(out(2))^2*out(2)];
        y=r+sqrt(sigmaw2/2)*(randn(2,1)+1i*randn(2,1));
        y0=L*x;
        
        % Update MSE estimate
        MSE_sim(:,ppp)=MSE_sim(:,ppp)+abs(y-y0).^2;
    end
    % Estimate NMSE for the non-linear system solution
    MSE_sim(:,ppp)=MSE_sim(:,ppp)/NbrOfRealizations./[sigmax2*gamma1^2; ...
        sigmax2*gamma2^2*betaa^2];
end

% Find the analytical NMSE values and optimal input reference power for the
% minimum of the maximum NMSEs of two branches
sigmax2=db2pow(-60:1:-20);
[nmse1, nmse2]=nmse_evaluate(xi,betaa,kappa1,kappa2,...
    sigmaw2,rho1,rho2,gamma1,gamma2,sigmax2);
sigmaStar_nmse=nmse_optimize(xi,betaa,kappa1,kappa2,...
    sigmaw2,rho1,rho2,gamma1,gamma2);
[nmse11,nmse22]=nmse_evaluate(xi,betaa,kappa1,kappa2,...
    sigmaw2,rho1,rho2,gamma1,gamma2,sigmaStar_nmse);

plot(pow2db(sigmax2)+30,10*log10(nmse1))
hold on
plot(pow2db(sigmax2)+30,10*log10(nmse2),'m--')
plot(pow2db(sigmaStar_nmse)+30, 10*log10(max(nmse11,nmse22)),'go','MarkerSize',20)
plot(pow2db(sigmax2)+30, 10*log10(MSE_sim(1,:)))
plot(pow2db(sigmax2)+30, 10*log10(MSE_sim(2,:)))

