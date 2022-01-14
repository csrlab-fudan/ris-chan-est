function [YTildeSage,zetaFft,xiFft,gammaFft,BXi,cFft] = Hd_est_fft_sage(RGammaFft,Mr,numTraining,zcLen,zcCpLen,YTilde,fftLenXi,fftLenZeta,fftLenGamma,dg)
% ***************************************
%  Coarse Solution to (88) Using FFTs
%  author - Xuemeng Zhou
%  input - fftLenXi: the length of the FFT transformation of \xi\bar
%          fftLenZeta: the FFT length of ζ\bar
%          fftLenGamma: the FFT length of θ\bar
%  output - zetaFft: the coarse estimate of ζ\bar
%           xiFft: the coarse estimate of \xi\bar
%           gammaFft: the coarse estimate of θ\bar
%  copyright - CSRL@Fudan,2021
%  **************************************

%% 离线角度+尽量小计算量的FFTs(KMr小于M时)
RZetaFft = fft(YTilde,fftLenZeta,2);
R = zeros(fftLenGamma,fftLenZeta,numTraining);
for nn = 1 : numTraining
    R(:,:,nn) = RGammaFft*RZetaFft((nn-1)*Mr+(1:Mr),:);
end
RGammaZetaXiFft = fft(R,fftLenXi,3);
BFft = abs(RGammaZetaXiFft); 
[~,indx] = max(BFft(:)); 
[gammaIndx, zetaIndx, xiIndx] = ind2sub([fftLenGamma,fftLenZeta,fftLenXi],indx);

%%
fzeta = (zetaIndx-1)/fftLenZeta;  
if fzeta > 1/2 
    fzeta = fzeta - 1;
end    
zetaFft = fzeta;
fxi = (xiIndx-1)/fftLenXi;  
if fxi > 1/2
    fxi = fxi - 1;
end
xiFft =  fxi/zcCpLen;
fgamma = (gammaIndx-1)/fftLenGamma;  
if fgamma > 1/2
    fgamma = fgamma - 1;
end
gammaFft = acos(-1/dg*fgamma);

%%
BXi = kron(exp(1j*2*pi*xiFft*(0:numTraining-1).'*zcCpLen),eye(Mr));
cFft = exp(-1j*2*dg*pi*(0:Mr-1).'*cos(gammaFft)); 
dZeta = exp(-1j*2*pi*zetaFft*(-zcLen/2:zcLen/2-1).');
alphaTilde = (cFft'*BXi'*YTilde*dZeta)/(zcLen*numTraining*Mr);

%% SAGE 当前径的接收端信号
YTildeSage = alphaTilde*BXi*cFft*exp(1j*2*pi*zetaFft*(-zcLen/2:zcLen/2-1));
end