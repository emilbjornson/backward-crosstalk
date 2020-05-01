% This Matlab script generates Figure 11 in the paper:

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

rng(1702)

% Prepare arrays for the optimal SE, SE with conventional-MRT and
% distortion-aware MRT, respectively
optimumm=[];
Mrt=[];
daMrt=[];

% Go through different crosstalk power values

rangee2=-70:10:-40;
for kappa=10.^(rangee2/20)
    
    
    
    
    opt=0;
    mr=0;
    damr=0;
    
    % Go through different random channel realizations
    
    for triall=1:10000
        
        % Amplification gain of the amplifiers
        gamma=10.^(30/20);
        gamma1=gamma;
        gamma2=gamma;
        
        
        % Receiver noise variance
        sigman2=db2pow(0);
        
        % Correlation coefficient between x_1 and x_2
        xi=0;
        
        % Ratio of second input signal's amplitude to the first one
        betaa=1;
        
        % Crosstalk parameters \kappa_1 and \kappa_2
        kappa1=kappa;
        kappa2=kappa;
        
        % Thermal noise variance
        sigmaw2=db2pow(-40);
        
        % Random channel realizations
        h1=sqrt(0.5)*(randn+1j*randn);
        h2=sqrt(0.5)*(randn+1j*randn);
        
        % Compression Parameters
        rho1=-0.025;
        rho2=-0.025;
        
        % tilde{h}_1 and \tilde{h}_2 in the paper
        htilde1=2*rho1*h1;
        htilde2=2*rho2*h2;
        
        % Q matrix in the paper
        Q=[gamma1 gamma1*gamma2*kappa2; gamma1*gamma2*kappa1 gamma2];
        
        
        
        % Find the optimal precoder weights in Theorem 2 (\tilde{c}_1 and \tilde{c}_2)
        % and the corresponding optimal SE
        [c1star,c2star, SEstar]=se_optimize(sigmaw2,rho1,rho2,h1,h2,sigman2);
        
        % Conventional MRT
        ccc=[h1'; h2'];
        
        ddd=Q*ccc;
        cbar=ddd/abs(h1);
        k0=h1*cbar(1)+h2*cbar(2);
        k1=htilde1*abs(cbar(1))^2*cbar(1)...
            +htilde2*abs(cbar(2))^2*cbar(2);
        sigma2=2*sigmaw2*(abs(h1)^2+abs(h2)^2)+2*sigman2;
        
        % Find the roots of (75)
        kokler=roots([2*abs(k1)^2*real(k0*k1'),...
            2*abs(k1)^2*abs(k0)^2,...
            -3*abs(k1)^2*sigma2, ...
            -4*real(k0*k1')*sigma2, ...
            -abs(k0)^2*sigma2]);
        kokler=kokler(real(kokler)>0);
        
        % Find the limit for the candidate solutions
        % to guarantee that the Bussgang gains are non-negative,
        % i.e., input power is not too high
        limitt=min(1/abs(2*abs(cbar(1))^2*rho1),1/abs(2*abs(cbar(2))^2*rho2));
        
        kokler=kokler(real(kokler)<limitt);
        kokler=real(kokler(abs(imag(kokler))./abs(real(kokler))<0.0001));
        
        % Go through all the candidate solutions (roots of (75)) to find the
        % optimal SE
        
        maxx=0;
        
        for candi=kokler.'
            
            funcc=2*candi*(abs(k1)^2*candi^2+2*real(k0*k1')*candi...
                +abs(k0)^2)/(abs(k1)^2*candi^3+sigma2);
            if funcc>maxx
                sigmaStarMRT=candi;
                SEstarMRT=log2(1+funcc);
                maxx=funcc;
            end
        end
        
        % Distortion-aware MRT
        SE_DA_MRT=[];
        xeksen2=[];
        for eta=0.0001:0.001:50
            
            A1=2*rho1*abs(h1)*sqrt(eta);
            B1=-1;
            C1=abs(h1)*sqrt(eta);
            
            A2=2*rho2*abs(h2)*sqrt(eta);
            B2=-1;
            C2=abs(h2)*sqrt(eta);
            
            c1=(-B1-sqrt(B1^2-4*A1*C1))/(2*A1)*exp(1j*angle(h1'));
            
            c2=(-B2-sqrt(B2^2-4*A2*C2))/(2*A2)*exp(1j*angle(h2'));
            ddd=inv(Q)*[c1; c2];
            SE=se_evalaute(sigmaw2,rho1,rho2,h1,h2,sigman2,c1,c2);
            SE_DA_MRT=[SE_DA_MRT SE];
            xeksen2=[xeksen2 abs(ddd(1))^2];
        end
        
        
        opt=opt+SEstar;
        mr=mr+SEstarMRT;
        damr=damr+max(SE_DA_MRT);
    end
    % Average SE values
    optimumm=[optimumm opt/10000];
    Mrt=[Mrt mr/10000];
    daMrt=[daMrt damr/10000];
end

% Plot optimumm, Mrt, and daMrt to obtain Fig. 11.