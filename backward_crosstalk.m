% This Matlab script is the function "backward_crosstalk" used to generate
% some figures in the paper:
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

function F = backward_crosstalk(u,x,L,K,G)

F(1) = L(1,1)*x(1)+K(1,2)*(u(2)+G(2,2)*abs(u(2))^2*u(2))-u(1);
F(2) = L(2,2)*x(2)+K(2,1)*(u(1)+G(1,1)*abs(u(1))^2*u(1))-u(2);
