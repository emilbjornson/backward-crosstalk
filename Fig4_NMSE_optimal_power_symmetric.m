% This Matlab script generates Figure 4 in the paper:

% Peter H�ndel, �zlem Tugfe Demir, Emil Bj�rnson, and Daniel R�nnow,
% "Impact of Backward Crosstalk in 2 � 2 MIMO Transmitters on NMSE and Spectral Efficiency,"
% IEEE Transactions on Communications, To appear.
%
% Download article: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9072380
%
% This is version 1.0 (Last edited: 2020-04-28)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

% Prepare arrays for the optimal input reference power and the
% corresponding optimal min-max NMSE
OptInputPower=[];
MinMaxNMSE=[];

% Go through different values of the amplification gain
for gamma=10.^((2:2:36)/20)
    gamma1=gamma;
    gamma2=gamma;
    
    % Correlation coefficient between x_1 and x_2
    xi=0;
    
    % Ratio of second input signal's amplitude to the first one
    betaa=1;
    
    % Crosstalk parameters \kappa, change it to -50 dB for the other
    % curves in the figure
    kappa=10^(-60/20);
    kappa1=kappa;
    kappa2=kappa;
    
    % Thermal noise variance
    sigmaw2=db2pow(-40);
    
    % Compression Parameters
    rho1=-0.025;
    rho2=-0.025;
    
    
    % Find the optimal input reference power and the corresponding NMSE values
    sigmaStar_nmse=nmse_optimize(xi,betaa,kappa1,kappa2,...
        sigmaw2,rho1,rho2,gamma1,gamma2);
    [nmse11,nmse22]=nmse_evaluate(xi,betaa,kappa1,kappa2,...
        sigmaw2,rho1,rho2,gamma1,gamma2,sigmaStar_nmse);
    
    OptInputPower=[OptInputPower pow2db(sigmaStar_nmse)+30];
    MinMaxNMSE=[MinMaxNMSE 10*log10(max(nmse11,nmse22))];
end
yyaxis left
plot(2:2:36,OptInputPower)
hold on
yyaxis right
plot(2:2:36,MinMaxNMSE,'m--')

