function [HDfSum,HDfPreSum,HDfPreNewSum,HDfPreNew2Sum,zcRootD,c,tauDListOut,alphaDListOut,gammaDListOut,xiDListOut] = multi_signal_Hd(Ud,upSampRate,roCoeff,Tp,zcCpLen,cpLen,lTilde,zcRootOri,Mr,numTrainingD,dg)
% ***************************************
%  generate channel from UE to BS
%  author - Xuemeng Zhou
%  input - Ud: the multipath number of Hd
%          upSampRate=1: sampling rate
%          roCoeff: rolloff coefficient 
%          Tp: the half-length of the pulse shaper
%          zcCpLen: the length of cp+lTilde
%          Mr: the antenna number of BS
%          dg：normalized inter-element distance for BS
%  output - HDfSum: the frequency-domain channel response of the estimated channel from UE to BS
%           HDfPreSum: the frequency-domain channel response of the projected channel from UE to BS after 10 OFDM symbols
%           HDfPreNewSum: the frequency-domain channel response of the projected channel from UE to BS after 20 OFDM symbols
%           HDfPreNew2Sum: the frequency-domain channel response of the projected channel from UE to BS after 40 OFDM symbols
%           gammaDListOut: the angle set of the direct link
%  copyright - CSRL@Fudan,2021
%  **************************************

%% 直射径参数设置
gammaDList = deg2rad([80 110]);
tauDList = [0.1 3.5];
alphaDList = [exp(1j*2*pi*rand) 0.8*exp(1j*2*pi*rand)]; 
xiDList = [6e-5 -9e-5];
c = zeros(Mr,Ud);
zcRootD = zeros(Ud,lTilde);

%% UE2BS信号
% rxSig1 = zeros(Mr,lTilde,Ud);
% rxSig2 = zeros(Mr,lTilde,Ud);
% rxSig3 = zeros(Mr,lTilde,Ud);
% rxSig4 = zeros(Mr,lTilde,Ud);

for ii = 1:Ud
    c(:,ii) = exp(-1j*pi*(0:Mr-1).'*cos(gammaDList(ii)));
    rcfCoeffD = raised_cosine(upSampRate,roCoeff,tauDList(ii),Tp); 
    zcRootOriRCD = conv(rcfCoeffD.',zcRootOri);
    zcRootCpD = zcRootOriRCD(end/2+(-zcCpLen/2+1:zcCpLen/2));
    zcRootD(ii,:) = zcRootCpD(cpLen+1:end);

%     rxSigTmp1 = c(:,ii)* (alphaDList(ii)*zcRootD(ii,:).*exp(1j*2*pi*xiDList(ii)/upSampRate*(numTrainingD*zcCpLen+(-lTilde/2:lTilde/2-1))));
%     rxSigTmp2 = rxSigTmp1*exp(1j*2*pi*xiDList(ii)/upSampRate*(10*zcCpLen));
%     rxSigTmp3 = rxSigTmp1*exp(1j*2*pi*xiDList(ii)/upSampRate*(20*zcCpLen));
%     rxSigTmp4 = rxSigTmp1*exp(1j*2*pi*xiDList(ii)/upSampRate*(40*zcCpLen)); 
% 
%     rxSig1(:,:,ii) = rxSigTmp1; 
%     rxSig2(:,:,ii) = rxSigTmp2;
%     rxSig3(:,:,ii) = rxSigTmp3;
%     rxSig4(:,:,ii) = rxSigTmp4;
end

%% 信道估计及预测（全频）
HDSum = 0;HDPreSum = 0;HDPreNewSum = 0;HDPreNew2Sum = 0;
for ii = 1:Ud
    c(:,ii) = exp(-1j*2*dg*pi*(0:Mr-1).'*cos(gammaDList(ii)));
    
    rcfCoeff = raised_cosine(upSampRate,roCoeff,tauDList(ii),Tp);
    rcfCoeffH = rcfCoeff.'; 

    HDTmp1 = c(:,ii)* alphaDList(ii)*rcfCoeffH*exp(1j*2*pi*xiDList(ii)*numTrainingD*zcCpLen); 
    HDPreTmp2 = HDTmp1*exp(1j*2*pi*xiDList(ii)*10*zcCpLen); 
    HDPreNewTmp3 = HDTmp1*exp(1j*2*pi*xiDList(ii)*20*zcCpLen); 
    HDPreNew2Tmp4 = HDTmp1*exp(1j*2*pi*xiDList(ii)*40*zcCpLen); 

    HDSum = HDSum + HDTmp1; 
    HDPreSum = HDPreSum + HDPreTmp2;
    HDPreNewSum = HDPreNewSum + HDPreNewTmp3;
    HDPreNew2Sum = HDPreNew2Sum + HDPreNew2Tmp4;
end
HDfSum = fft(HDSum,lTilde,2);
HDfPreSum = fft(HDPreSum,lTilde,2);
HDfPreNewSum = fft(HDPreNewSum,lTilde,2);
HDfPreNew2Sum = fft(HDPreNew2Sum,lTilde,2);
%%
tauDListOut = tauDList(1:Ud);
alphaDListOut = alphaDList(1:Ud);
gammaDListOut = gammaDList(1:Ud);
xiDListOut = xiDList(1:Ud);

end
