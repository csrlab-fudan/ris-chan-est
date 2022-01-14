function [YTildeSage,zetaFft,xiFft,thetaFft,phiFft,WXi,aRISFft] = allD_est_fft_sage(TFft,AangleFft,W,M,Mr,numTraining,zcLen,zcCpLen,YTilde,fftLenXi,fftLenZeta,fftLen1,fftLen2,dx,dz,P,Q)
% ***************************************
%  Coarse Solution to (39) Using FFTs
%  author - Xuemeng Zhou
%  input - fftLenXi: the length of the FFT transformation of \xi
%          fftLenZeta: the FFT length of ζ
%          fftLen1 fftLen2: the FFT length of θand \varphi
%  output - zetaFft: the coarse estimate of ζ
%           xiFft: the coarse estimate of \xi
%           thetaFft: the coarse estimate of θ
%           phiFft: the coarse estimate of \varphi
%  copyright - CSRL@Fudan,2021
%  **************************************

%% 利用离线信息化简角度计算 
RZetaFft = fft(YTilde,fftLenZeta,2);
R = zeros(fftLen1*fftLen2,fftLenZeta,numTraining);
for nn = 1 : numTraining
    R(:,:,nn) = TFft(:,(nn-1)*Mr+(1:Mr))*RZetaFft((nn-1)*Mr+(1:Mr),:);
end
RZetaXiFft = fft(R,fftLenXi,3);
BangleFftTmp = reshape(RZetaXiFft,fftLen1,fftLen2,fftLenZeta,fftLenXi);
BangleFft = abs(BangleFftTmp).^2;
angleFft = BangleFft./AangleFft;
[~,indx] = max(angleFft(:)); 
[YIndx, XIndx, zetaIndx, xiIndx] = ind2sub([fftLen1,fftLen2,fftLenZeta,fftLenXi],indx);

%%
fzeta = (zetaIndx-1)/fftLenZeta;  %k/N
if fzeta > 1/2 
    fzeta = fzeta - 1;
end    
zetaFft = fzeta;

fxi = (xiIndx-1)/fftLenXi;  %k/N
if fxi > 1/2
    fxi = fxi - 1;
end
xiFft =  fxi/zcCpLen;

fx = (XIndx-1)/fftLen2;  %k/N
if fx > 1/2
    fx = fx - 1;
end
if fx == 0.5 
    fx = fx - eps;
end
phiFft = acos(2*fx);

fy = (YIndx-1)/fftLen1;  %k/N
if fy > 1/2
    fy = fy - 1;
end    
thetaFft = real(acos(-2*fy/sin(phiFft))); 

%%
WXi = zeros(Mr*numTraining,M);
for nn = 1 : numTraining
    WXi((nn-1)*Mr+(1:Mr),:) = W((nn-1)*Mr+(1:Mr),:)*exp(1j*2*pi*xiFft*(nn-1)*zcCpLen);
end

ARISFft = generate_ar(dx, dz, P, Q, thetaFft, phiFft);
aRISFft = reshape(ARISFft.',[],1);
dZeta = exp(-1j*2*pi*zetaFft*(-zcLen/2:zcLen/2-1).');
betaTilde = (aRISFft'*WXi'*YTilde*dZeta)/(zcLen*norm(W*aRISFft)^2);

%% SAGE 当前径的接收端信号
YTildeSage = betaTilde*WXi*aRISFft*exp(1j*2*pi*zetaFft*(-zcLen/2:zcLen/2-1));
end