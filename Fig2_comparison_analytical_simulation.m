% This Matlab script generates Figure 2 in the paper:

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
% To track the exitflag of the "fsolve"
exitt=[];

% Initialize the NMSE values for cross covariances between u_1 and u_2
% obtained by simulation (non-linear system solution) and the analytical
% result
nmse_sim=[];

% Correlation coefficient between x_1 and x_2
xi=0;

% Ratio of second input signal's amplitude to the first one
betaa=1;

% Go through different values of the input reference power P_x

for sigmax2=db2pow([-50 -40 -30])
    
    % Crosstalk parameters
    kappa=10.^(-50/20);
    kappa1=kappa;
    kappa2=kappa;
    
    % Thermal noise variance
    sigmaw2=db2pow(-40);
    
    % Compression Parameters
    rho1=-0.025;
    rho2=-0.025;
    
    % Construct the matrices in the paper
    G=diag([rho1,rho2]);
    L=diag([gamma1,gamma2]);
    K=[0, gamma1*kappa2; gamma2*kappa1, 0];
    Q=[gamma1, gamma1*gamma2*kappa2; gamma1*gamma2*kappa1, gamma2];
    
    % Number of random input signal realizations for the CDF curve
    NbrOfRealizations=10000;
    
    % Prepare the values of u_1 and u_2 obtained by simulation
    U_sim=zeros(2,NbrOfRealizations);
    
    
    options = optimset('Display','off');
    
    % Input covariance matrix
    Cx=sigmax2*[1, betaa*xi; betaa*xi', betaa^2];
    
    % Generate white Gaussian unit-variance random variables
    W=sqrt(0.5)*(randn(2, NbrOfRealizations)+1i*randn(2, NbrOfRealizations));
    
    % Generate random input with covariance matrix Cx
    X=sqrtm(Cx)*W;
    
    % The analytical values for u_1 and u_2 in the paper
    U_approx=Q*X;
    
    % Go through all the random input realizations and solve the non-linear
    % system to find the simulation values for u_1 and u_2
    for nnn=1:NbrOfRealizations
        x=X(:,nnn);
        f = @(u) backward_crosstalk(u,x,L,K,G);
        [out,fval,exitflag]=fsolve(f,Q*x,options);
        exitt=[exitt exitflag];
        U_sim(:,nnn)=out;
    end
    
    % Plot the CDF curves
    ecdf([real(vec(U_approx(1,:)));imag(vec(U_approx(1,:)));...
        real(vec(U_approx(2,:)));imag(vec(U_approx(2,:))) ])
    hold on
    ecdf([real(vec(U_sim(1,:)));imag(vec(U_sim(1,:)));...
        real(vec(U_sim(2,:)));imag(vec(U_sim(2,:))) ])
    
    % Estimate the cross covariance between u_1 and u_2
    cross_var=vec(U_sim(2,:))'*vec(U_sim(1,:))/NbrOfRealizations;
    
    % The analytical covariance matrix U for u_1 and u_2 in the paper
    UUU=Q*Cx*Q';
    
    % Calculate the NMSE between cross covariance estimate and the analytical
    % one
    nmse_sim=[nmse_sim abs(cross_var-UUU(1,2))^2/abs(cross_var)^2];
end