% This Matlab script is the function "se_evaluate" used to generate
% some figures and evaluare the SE in the paper:
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
function SE=se_evaluate(sigmaw2,rho1,rho2,h1,h2,sigman2,c1,c2)

htilde1=2*h1*rho1;
htilde2=2*h2*rho2;
sigma2=2*sigmaw2*(abs(h1)^2+abs(h2)^2)+2*sigman2;
NUM=2*abs(h1*c1+h2*c2+htilde1*abs(c1)^2*c1+htilde2*abs(c2)^2*c2)^2;
DEN=abs(htilde1*abs(c1)^2*c1+htilde2*abs(c2)^2*c2)^2+sigma2;
SE=log2(1+NUM/DEN);

