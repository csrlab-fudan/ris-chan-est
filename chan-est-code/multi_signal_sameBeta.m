function [HfSumNew,HfPreSumNew,HfPreNewSumNew,HfPreNew2SumNew,zcRoot1,aRIS,tauListOut,betaListOut,thetaListOut,phiListOut,cfoListOut] = multi_signal_sameBeta(U,upSampRate,roCoeff,Tp,zcCpLen,cpLen,lTilde,zcRootOri,dx,dz,P,Q,M)
% ***************************************
%  generate channel from UE to RIS
%  author - Xuemeng Zhou
%  input - U: the multipath number of H
%          upSampRate=1: sampling rate
%          roCoeff: rolloff coefficient 
%          Tp: the half-length of the pulse shaper
%          zcCpLen: the length of cp+lTilde
%          M: RIS element number
%          dx,dz：normalized inter-element distance for RIS
%  output - HfSumNew: the frequency-domain channel response of the estimated channel from UE to RIS
%           HfPreSumNew: the frequency-domain channel response of the projected channel from UE to RIS after 10 OFDM symbols
%           HfPreNewSumNew: the frequency-domain channel response of the projected channel from UE to RIS after 20 OFDM symbols
%           HfPreNew2SumNew: the frequency-domain channel response of the projected channel from UE to RIS after 40 OFDM symbols
%           cfoListOut: the frequency offset set of the reflection link
%  copyright - CSRL@Fudan,2021
%  **************************************

%% 反射径参数设置
tauList = [0.5 1.1 1.8 2.6 5.5 6.4];
betaList = [0.85*exp(1j*2*pi*rand),0.8*exp(1j*2*pi*rand),0.7*exp(1j*2*pi*rand),0.65*exp(1j*2*pi*rand),0.6*exp(1j*2*pi*rand),0.5*exp(1j*2*pi*rand)];
thetaList = deg2rad([90 40 65 150 100 120 90]); 
phiList = deg2rad([60 30 50 60 45 40 55]); 
cfoList = [3e-6 3e-4 1e-5 -6e-6 1e-4 4e-6 -2e-5];

%% UE2RIS信号
zcRoot1 = zeros(U,lTilde);
aRIS = zeros(M,U);

for ii = 1:U
    ARIS = generate_ar(dx, dz, P, Q, thetaList(ii), phiList(ii));
    aRIS(:,ii) = reshape(ARIS.',[],1); 
    rcfCoeff = raised_cosine(upSampRate,roCoeff,tauList(ii),Tp); 
    zcRootOriRC1 = conv(rcfCoeff.',zcRootOri);
    zcRootCp1 = zcRootOriRC1(end/2+(-zcCpLen/2+1:zcCpLen/2));
    zcRoot1(ii,:) = zcRootCp1(cpLen+1:end);   
end

%% 多径信道恢复（全频）：
HSum = 0;HPreSum = 0;HPreNewSum = 0;HPreNew2Sum = 0;
for ii = 1:U
    rcfCoeff = raised_cosine(upSampRate,roCoeff,tauList(ii),Tp);
    rcfCoeffH = rcfCoeff.'; 
   
    HTmp1 = aRIS(:,ii)* betaList(ii)*rcfCoeffH;
    HPreTmp2 = aRIS(:,ii)* betaList(ii)*rcfCoeffH*exp(1j*2*pi*cfoList(ii)*10*zcCpLen);
    HPreNewTmp3 = aRIS(:,ii)* betaList(ii)*rcfCoeffH*exp(1j*2*pi*cfoList(ii)*20*zcCpLen); 
    HPreNew2Tmp4 = aRIS(:,ii)* betaList(ii)*rcfCoeffH*exp(1j*2*pi*cfoList(ii)*40*zcCpLen); 

    HSum = HSum + HTmp1; 
    HPreSum = HPreSum + HPreTmp2;
    HPreNewSum = HPreNewSum + HPreNewTmp3;
    HPreNew2Sum = HPreNew2Sum + HPreNew2Tmp4;
end
HfSumNew = fft(HSum,lTilde,2);
HfPreSumNew = fft(HPreSum,lTilde,2);
HfPreNewSumNew = fft(HPreNewSum,lTilde,2);
HfPreNew2SumNew = fft(HPreNew2Sum,lTilde,2);

%%
tauListOut = tauList(1:U);
betaListOut = betaList(1:U);
thetaListOut = thetaList(1:U);
phiListOut = phiList(1:U);
cfoListOut = cfoList(1:U);

end
