function [HEst,HPreEst,HPreEst3,HPreEst4,YTildeSage,zetaNewton,xiNewton,aRISNewton,thetaNewton,phiNewton,tauNewton,beta] =...
    allD_est_newton_sage(zetaNewton,xiNewton,aRISNewton,thetaNewton,phiNewton,YTilde,zcRootOri,W,zcLen,...
    zcCpLen,cpLen,lTilde,P,Q,dx,dz,Mr,numTraining,Tp,maxNewtonIter,IterMax,minNewtonStepSize,upSampRate,roCoeff)
% ***************************************
%  Refined Channel Estimation Using Newton’s Method
%  author - Xuemeng Zhou
%  input - numTraining: the observation times of H
%  output - HEst: the time-domain channel response of the estimated channel from UE to RIS at present
%           HPreEst: the time-domain channel response of the projected channel from UE to RIS after 10 OFDM symbols
%           HPreEst3: the time-domain channel response of the projected channel from UE to RIS after 20 OFDM symbols
%           HPreEst4: the time-domain channel response of the projected channel from UE to RIS after 40 OFDM symbols
%  copyright - CSRL@Fudan,2021
%  **************************************

%% 牛顿迭代优化
for ii = 1 : IterMax
    %% 二维频偏+二维角度
    [zetaNewton,xiNewton,WXi] = zetaXi_est_newton(W,YTilde,Mr,numTraining,zcLen,zcCpLen,zetaNewton,xiNewton,aRISNewton,maxNewtonIter,minNewtonStepSize);
    [thetaNewton,phiNewton,aRISNewton] = angle_est_newton(zetaNewton,WXi,W,zcLen,YTilde,P,Q,dx,dz,thetaNewton,phiNewton,maxNewtonIter,minNewtonStepSize);

    xiIterNewton(ii) = xiNewton;
    thetaIterNewton(ii) = thetaNewton;
    phiIterNewton(ii) = phiNewton;
    if ii >= 2
        if xiIterNewton(ii)-xiIterNewton(ii-1) < 1e-10 && phiIterNewton(ii)-phiIterNewton(ii-1) < 1e-10 && thetaIterNewton(ii)-thetaIterNewton(ii-1) < 1e-10
            break
        end
    end
end
%% 频域恢复信道
tauNewton = (xiNewton-zetaNewton)*lTilde;
dZeta = exp(-1j*2*pi*zetaNewton*(-zcLen/2:zcLen/2-1).');
betaTilde = (aRISNewton'*WXi'*YTilde*dZeta)/(zcLen*norm(W*aRISNewton)^2);
beta = betaTilde/exp(1j*pi*tauNewton^2/lTilde);

rcfCoeffTmp = raised_cosine(upSampRate,roCoeff,tauNewton,Tp);
rcfCoeffH = rcfCoeffTmp.';
HEst = aRISNewton*beta*rcfCoeffH; 
HPreEst = aRISNewton*beta*rcfCoeffH*exp(1j*2*pi*xiNewton*10*zcCpLen); 
HPreEst3 = aRISNewton*beta*rcfCoeffH*exp(1j*2*pi*xiNewton*20*zcCpLen); 
HPreEst4 = aRISNewton*beta*rcfCoeffH*exp(1j*2*pi*xiNewton*40*zcCpLen); 

%% SAGE 接收端信号
YTildeSage = betaTilde*WXi*aRISNewton*exp(1j*2*pi*zetaNewton*(-zcLen/2:zcLen/2-1));
end

