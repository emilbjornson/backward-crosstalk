% This Matlab script is the function "se_optimize" used to generate
% some figures and finds the optimal precoder weights in Theorem 2
% maximizing the SE in the paper:
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

function [c1star,c2star, maxx]=se_optimize(...
    sigmaw2,rho1,rho2,h1,h2,sigman2)

% \tilde{h}_1, \tilde{h}_2, and \sigma^2 in the paper
htilde1=2*h1*rho1;
htilde2=2*h2*rho2;
sigma2=2*sigmaw2*(abs(h1)^2+abs(h2)^2)+2*sigman2;

% Prepare arrays for the candidate solutions for different cases
candidates1=[];
candidates2=[];

% Case1-B
candidates1=[candidates1 0 0];
candidates2=[candidates2 sqrt(-1/(2*rho2))];
kokler=roots([2*abs(htilde2)^2, 0,  -6*rho2*sigma2, -sigma2]);
kokler=kokler(real(kokler)>0);
[~,indexx]=min(abs(imag(kokler)));
candidates2=[candidates2 sqrt(real(kokler(indexx)))];

% Case1-C
candidates1=[candidates1 sqrt(-1/(2*rho1))];
candidates2=[candidates2 0 ];

candidates1=[candidates1 sqrt(-1/(2*rho1))*exp(1j*angle(h1'*h2))];
candidates2=[candidates2 sqrt(abs(rho1/rho2))*sqrt(-1/(2*rho1))];

candidates1=[candidates1 sqrt(-1/(2*rho1))*exp(1j*angle(h1'*h2)+1j*pi)];
candidates2=[candidates2 sqrt(abs(rho1/rho2))*sqrt(-1/(2*rho1))];

kokler=roots([2*abs(htilde1)^2, 0,  -6*rho1*sigma2, -sigma2]);
kokler=kokler(real(kokler)>0);
[~,indexx]=min(abs(imag(kokler)));

candidates1=[candidates1 sqrt(real(kokler(indexx)))];
candidates2=[candidates2 0 ];

% Case 1-D
AAA=abs(htilde1)+abs(htilde2)*abs(rho1)/abs(rho2)*sqrt(abs(rho1/rho2));
kokler=roots([2*AAA^2, 0,  -6*rho1*sigma2, -sigma2]);
kokler=kokler(real(kokler)>0);
[~,indexx]=min(abs(imag(kokler)));


candidates1=[candidates1 sqrt(real(kokler(indexx)))*exp(1j*angle(h1'*h2))];
candidates2=[candidates2 sqrt(abs(rho1/rho2))*sqrt(real(kokler(indexx)))];


% Case 2
AAA=abs(htilde1)-abs(htilde2)*abs(rho1)/abs(rho2)*sqrt(abs(rho1/rho2));
kokler=roots([2*AAA^2, 0,  -6*rho1*sigma2, -sigma2]);
kokler=kokler(real(kokler)>0);
[~,indexx]=min(abs(imag(kokler)));

candidates1=[candidates1 sqrt(real(kokler(indexx)))*exp(1j*angle(h1'*h2)+1j*pi)];
candidates2=[candidates2 sqrt(abs(rho1/rho2))*sqrt(real(kokler(indexx)))];

% Find the optimal solution among the candidates
maxx=0;
for ttt=1:length(candidates1)
    c1=candidates1(ttt);
    c2=candidates2(ttt);
    NUM=2*abs(h1*c1+h2*c2+htilde1*abs(c1)^2*c1+htilde2*abs(c2)^2*c2)^2;
    DEN=abs(htilde1*abs(c1)^2*c1+htilde2*abs(c2)^2*c2)^2+sigma2;
    SE=log2(1+NUM/DEN);
    if SE>maxx
        maxx=SE;
        c1star=c1;
        c2star=c2;
    end
end
