function [HDEst,HDPreEst,HDPreEst3,HDPreEst4,YDTildeSage,YDTildeSageTest,zetaNewton,xiNewton,gammaNewton,tauNewton,alpha] =...
    Hd_est_newton_sage(zetaNewton,xiNewton,gammaNewton,YTilde,zcRootOri,zcLen,...
    zcCpLen,cpLen,lTilde,Mr,numTraining,numTrainingD,Tp,maxNewtonIter,IterMax,minNewtonStepSize,upSampRate,roCoeff,dg)
% ***************************************
%  Refined Channel Estimation Using Newton’s Method
%  author - Xuemeng Zhou
%  input - numTraining: the observation times of H
%          numTrainingD: the observation times of Hd
%  output - HDEst: the time-domain channel response of the estimated channel from UE to BS at present
%           HDPreEst: the time-domain channel response of the projected channel from UE to BS after 10 OFDM symbols
%           HDPreEst3: the time-domain channel response of the projected channel from UE to BS after 20 OFDM symbols
%           HDPreEst4: the time-domain channel response of the projected channel from UE to BS after 40 OFDM symbols
%  copyright - CSRL@Fudan,2021
%  **************************************

%% 牛顿迭代优化
for ii = 1 : IterMax
    %% 三维（频偏+角度）
    [zetaNewton,xiNewton,gammaNewton,cNewton,BXi] = zetaXiGamma_est_newton(YTilde,Mr,numTrainingD,zcLen,zcCpLen,zetaNewton,xiNewton,gammaNewton,maxNewtonIter,minNewtonStepSize,dg);
    zetaIterNewton(ii) = zetaNewton;
    xiIterNewton(ii) = xiNewton;
    gammaIterNewton(ii) = gammaNewton;
    if ii >= 2
        if xiIterNewton(ii)-xiIterNewton(ii-1) < 1e-10 && zetaIterNewton(ii)-zetaIterNewton(ii-1) < 1e-10 && gammaIterNewton(ii)-gammaIterNewton(ii-1) < 1e-10
            break
        end
    end
end

%% 频域恢复信道
tauNewton = (xiNewton-zetaNewton)*lTilde;
dZeta = exp(-1j*2*pi*zetaNewton*(-zcLen/2:zcLen/2-1).');
alphaTilde = (cNewton'*BXi'*YTilde*dZeta)/(zcLen*numTrainingD*Mr);
alpha = alphaTilde/exp(1j*pi*tauNewton^2/lTilde);

rcfCoeffTmp = raised_cosine(upSampRate,roCoeff,tauNewton,Tp);
rcfCoeffH = rcfCoeffTmp.';
HDEst = cNewton*alpha*rcfCoeffH*exp(1j*2*pi*xiNewton*numTrainingD*zcCpLen);
HDPreEst = HDEst*exp(1j*2*pi*xiNewton*10*zcCpLen); 
HDPreEst3 = HDEst*exp(1j*2*pi*xiNewton*20*zcCpLen);
HDPreEst4 = HDEst*exp(1j*2*pi*xiNewton*40*zcCpLen);

%% SAGE 当前径的接收端信号
YDTildeSage = alphaTilde*BXi*cNewton*exp(1j*2*pi*zetaNewton*(-zcLen/2:zcLen/2-1));

%% 反射径估计准备
BXiTest = kron(exp(1j*2*pi*xiNewton*((0:numTraining-1)+numTrainingD).'*zcCpLen),eye(Mr));
YDTildeSageTest = alphaTilde*BXiTest*cNewton*exp(1j*2*pi*zetaNewton*(-zcLen/2:zcLen/2-1));

end

