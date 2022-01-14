clear;
%%
snrDbSet = 0;  
snrDbSetLen = length(snrDbSet);
numMonte = 100;
IterMax = 20;
maxSageIter = 15; 
maxNewtonIter = 20;
maxFftIter = 20;
minNewtonStepSize = 1e-8;
fftLenZeta = 1024;
fftLenXi = 64;
fftLen1 = 32;
fftLen2 = 32;
fftLenGamma = 32;
numTrainingSet = [2,3,4,5,6,8,10]; 
numTrainingSetLen = length(numTrainingSet);
U = 4;

%% RIS参数设置
P = 32;                %RIS行数
Q = 32;                %RIS列数
M = P*Q;               %RIS端阵列数
Mr = 6;                %基站天线数，相当于原场景的基站端射频数
Mt = 4;                %UE天线数
p0 = Q;                %RIS列数
c = 3e8;               %光速 
fcHz = 28e9;           %载波频率
lambda = c/fcHz;
d0 = 1;                %RIS到基站距离
dx = 1/2;              % normalized dx
dz = 1/2;              % normalized dz
dg = 4; 
mu = 0.5;              %反射系数
GOri = channel_generation_RIS2BS(M,Mr,d0,mu,p0,lambda,dg,dx);
G = GOri.';            %Mr*M 
%% ZC序列参数
zcLen = 600;           % 截取低频部分的信号长度
lTilde = 1024;         % ZC原长度-lTilde/2至lTilde/2-1，偶数
cpLen = 64;            % 大于时延 
zcCpLen = cpLen + lTilde; 
r = 1;
rTilde = r/lTilde;
%% 升余弦函数（pulse shaper）
roCoeff = 0.3;  % 滚降系数
upSampRate = 1; % 上采样，插零
Tp = 15;        % 升余弦函数的半长

%% 理论发射信号
zcSeq = exp(1j*pi*rTilde*(-zcLen/2:zcLen/2-1).^2); 
zcRootOriTmp = exp(1j*pi*rTilde*(-lTilde/2:lTilde/2-1).^2);  
zcRootOri = [zcRootOriTmp(end-cpLen+1:end),zcRootOriTmp];  

%% 提前分配数组
mseTauTmp = zeros(1,numMonte); mseTau = zeros(numTrainingSetLen,snrDbSetLen); mseXiTmp = zeros(1,numMonte); mseXi = zeros(numTrainingSetLen,snrDbSetLen);
mseThetaTmp = zeros(1,numMonte); mseTheta = zeros(numTrainingSetLen,snrDbSetLen); msePhiTmp = zeros(1,numMonte); msePhi = zeros(numTrainingSetLen,snrDbSetLen);
mseBetaTmp = zeros(1,numMonte); mseBeta = zeros(numTrainingSetLen,snrDbSetLen);mseZetaTmp = zeros(1,numMonte); mseZeta = zeros(numTrainingSetLen,snrDbSetLen);

HResTmp = zeros(1,numMonte); HRes = zeros(numTrainingSetLen,snrDbSetLen); rmseHResSum = zeros(numTrainingSetLen,snrDbSetLen);
HPreResTmp = zeros(1,numMonte); HPreRes = zeros(numTrainingSetLen,snrDbSetLen); rmseHPreResSum = zeros(numTrainingSetLen,snrDbSetLen);
HPreResNewTmp = zeros(1,numMonte); HPreResNew = zeros(numTrainingSetLen,snrDbSetLen); rmseHPreResNewSum = zeros(numTrainingSetLen,snrDbSetLen);
HPreResNew2Tmp = zeros(1,numMonte); HPreResNew2 = zeros(numTrainingSetLen,snrDbSetLen); rmseHPreResNew2Sum = zeros(numTrainingSetLen,snrDbSetLen);
HPreResCompareTmp = zeros(1,numMonte); HPreResCompare = zeros(numTrainingSetLen,snrDbSetLen); rmseHPreResCompareSum = zeros(numTrainingSetLen,snrDbSetLen);
HPreResNewCompareTmp = zeros(1,numMonte); HPreResNewCompare = zeros(numTrainingSetLen,snrDbSetLen); rmseHPreResNewCompareSum = zeros(numTrainingSetLen,snrDbSetLen);

%% 产生信道H
[HfSum,HfPreSum,HfPreNewSum,HfPreNew2Sum,zcRoot1,aRIS,tauList,betaList,thetaList,phiList,xiList] = multi_signal_sameBeta(U,upSampRate,roCoeff,Tp,zcCpLen,cpLen,lTilde,zcRootOri,dx,dz,P,Q,M); % rx是未加噪声的接收信号
for tt = 1 : numTrainingSetLen
    tt
    numTraining = numTrainingSet(tt);
%% UE2RIS2BS
%%%%%%%%%% 移相器
W = zeros(Mr*numTraining,M);
for nn = 1 : numTraining 
    D = diag(exp(1j*2*pi/4*randi([1,4],M,1))); % 2bit精度移相器
    W((nn-1)*Mr+(1:Mr),:) = G*D;
end

W1 = W';
W23D = reshape(W1,Q,P,numTraining*Mr);
TFft1 = fft(W23D, fftLen1, 1); 
TFft2 = fft(TFft1, fftLen2, 2);
TFft = reshape(TFft2,fftLen1*fftLen2,numTraining*Mr); % 现fftLen1*fftLen2=M
AangleFft1 = (abs(TFft2)).^2;
AangleFft = sum(AangleFft1,3);

%%%%%%%% BS接收经RIS的信号
YSig = zeros(Mr*numTraining,lTilde);
for nn = 1 : numTraining
    gXi = exp(1j*2*pi*xiList*zcCpLen*(nn-1));
    YSig((nn-1)*Mr+(1:Mr),:) = W((nn-1)*Mr+(1:Mr),:)*aRIS*diag(betaList)*diag(gXi)*(zcRoot1.*exp(1j*2*pi*kron(xiList.',(-lTilde/2:lTilde/2-1))));
end

%%
for mm = 1 : numMonte
    if mod(mm, 10) == 0
        fprintf(['\nThe Monte Carol Run is %d\n'], mm);
    end
    for ss = 1 : snrDbSetLen
        snrDb = snrDbSet(ss);        
        addNoise1 = sqrt(1/2)*(randn(Mr*numTraining,lTilde)+randn(Mr*numTraining,lTilde)*1j)*(db2mag(-snrDb));
        Y1Ori = YSig + addNoise1;
        Y1 = Y1Ori(:,end/2+(-zcLen/2+1:zcLen/2));  % 低频部分
        YTilde = Y1*diag(conj(zcSeq));

        %% SAGE      
         sigPerPath = zeros(Mr*numTraining,zcLen,U);
         HEstSagePerPath = zeros(M,2*Tp+1,U); 
         HPreEstSagePerPath = zeros(M,2*Tp+1,U);
         HPreEstNewSagePerPath = zeros(M,2*Tp+1,U);
         HPreEstNew2SagePerPath = zeros(M,2*Tp+1,U);
         tauNewtonPerPath = zeros(1,U);
         zetaNewtonPerPath = zeros(1,U);
         xiNewtonPerPath = zeros(1,U);
         thetaNewtonPerPath = zeros(1,U);
         phiNewtonPerPath = zeros(1,U);
         aRISNewtonPerPath = zeros(M,1,U);
     
         %% 反射径
         %%% FFT
         for uu = 1 : U
             YTildeSage = YTilde - sum(sigPerPath,3) + sigPerPath(:,:,uu);
             [YTildeSage,zetaNewton,xiNewton,thetaNewton,phiNewton,WXi,aRISNewton] = allD_est_fft_sage(TFft,AangleFft,W,M,Mr,numTraining,zcLen,zcCpLen,YTildeSage,fftLenXi,fftLenZeta,fftLen1,fftLen2,dx,dz,P,Q);
             sigPerPath(:,:,uu) = YTildeSage;
             zetaNewtonPerPath(uu) = zetaNewton;
             xiNewtonPerPath(uu) = xiNewton;
             thetaNewtonPerPath(uu) = thetaNewton;
             phiNewtonPerPath(uu) = phiNewton;
             aRISNewtonPerPath(:,:,uu) = aRISNewton;
         end
         %%% Newton
         for gg = 1 : maxSageIter
             for uu = 1 : U
                 YTildeSage = YTilde - sum(sigPerPath,3) + sigPerPath(:,:,uu);
                 [HEstSage,HPreEstSage,HPreEstNewSage,HPreEstNew2Sage,YTildeSage,zetaNewtonPerPath(uu),xiNewtonPerPath(uu),aRISNewtonPerPath(:,:,uu),thetaNewtonPerPath(uu),phiNewtonPerPath(uu),tauNewtonPerPath(uu),beta] = ...
                     allD_est_newton_sage(zetaNewtonPerPath(uu),xiNewtonPerPath(uu),aRISNewtonPerPath(:,:,uu),thetaNewtonPerPath(uu),phiNewtonPerPath(uu),YTildeSage,...
                     zcRootOri,W,zcLen,zcCpLen,cpLen,lTilde,P,Q,dx,dz,Mr,numTraining,Tp,maxNewtonIter,IterMax,minNewtonStepSize,upSampRate,roCoeff);
                 sigPerPath(:,:,uu) = YTildeSage;
                 HEstSagePerPath(:,:,uu) = HEstSage;
                 HPreEstSagePerPath(:,:,uu) = HPreEstSage;
                 HPreEstNewSagePerPath(:,:,uu) = HPreEstNewSage;
                 HPreEstNew2SagePerPath(:,:,uu) = HPreEstNew2Sage;
             end
             
             YTildeSageTmp(gg) = norm(sum(sigPerPath,3)-YTilde, 'fro');
             %%% SAGE收敛条件
             if YTildeSageTmp(gg)./norm(sum(sigPerPath,3),'fro')<1e-8
                 break;
             end
                     
         end

         HEstfNew = fft(sum(HEstSagePerPath,3),lTilde,2);
         HPreEstfNew = fft(sum(HPreEstSagePerPath,3),lTilde,2);
         HPreEstf3New = fft(sum(HPreEstNewSagePerPath,3),lTilde,2);
         HPreEstf4New = fft(sum(HPreEstNew2SagePerPath,3),lTilde,2);
         
         %% H_UE2RIS
         HResSub = norm(HEstfNew-HfSum, 'fro')^2;
         HResTmp(mm) = HResSub/(M*lTilde);
         HRes(tt,ss) = HRes(tt,ss) + HResTmp(mm);
         
         HPreResSub = norm(HPreEstfNew-HfPreSum, 'fro')^2;
         HPreResTmp(mm) = HPreResSub/(M*lTilde);
         HPreRes(tt,ss) = HPreRes(tt,ss) + HPreResTmp(mm);
         
         HPreResNewSub = norm(HPreEstf3New-HfPreNewSum, 'fro')^2;
         HPreResNewTmp(mm) = HPreResNewSub/(M*lTilde);
         HPreResNew(tt,ss) = HPreResNew(tt,ss) + HPreResNewTmp(mm);
         
         HPreResNew2Sub = norm(HPreEstf4New-HfPreNew2Sum, 'fro')^2;
         HPreResNew2Tmp(mm) = HPreResNew2Sub/(M*lTilde);
         HPreResNew2(tt,ss) = HPreResNew2(tt,ss) + HPreResNew2Tmp(mm);
         
         HPreResCompareSub = norm(HPreEstfNew-HfSum, 'fro')^2;
         HPreResCompareTmp(mm) = HPreResCompareSub/(M*lTilde);
         HPreResCompare(tt,ss) = HPreResCompare(tt,ss) + HPreResCompareTmp(mm);
         
         HPreResNewCompareSub = norm(HPreEstf3New-HfSum, 'fro')^2;
         HPreResNewCompareTmp(mm) = HPreResNewCompareSub/(M*lTilde);
         HPreResNewCompare(tt,ss) = HPreResNewCompare(tt,ss) + HPreResNewCompareTmp(mm);

    end
end
%%% H
rmseHResSum(tt,:) = sqrt(HRes(tt,:)/numMonte);
rmseHPreResSum(tt,:) = sqrt(HPreRes(tt,:)/numMonte);
rmseHPreResNewSum(tt,:) = sqrt(HPreResNew(tt,:)/numMonte);
rmseHPreResNew2Sum(tt,:) = sqrt(HPreResNew2(tt,:)/numMonte);
rmseHPreResCompareSum(tt,:) = sqrt(HPreResCompare(tt,:)/numMonte);
rmseHPreResNewCompareSum(tt,:) = sqrt(HPreResNewCompare(tt,:)/numMonte);

end

%% 观测次数
figure;
semilogy(numTrainingSet, rmseHResSum(:,1), 'r-o','lineWidth',1.5);hold on 
xlabel('The numbers of observation (K)','FontSize',14)
ylabel('{RMSE of channel estimation}','FontSize',14);
title('SNR=0dB,U=4','FontSize',14); 
set(gca,'FontSize',12);
grid on;
