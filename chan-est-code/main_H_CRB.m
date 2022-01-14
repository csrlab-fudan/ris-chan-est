clear;
%%
snrDbSet = -20:6:10;  
snrDbSetLen = length(snrDbSet);
numMonte = 100; 
IterMax = 20;
maxSageIter = 15; 
maxNewtonIter = 20;
maxFftIter = 20;
minNewtonStepSize = 1e-8;1e-10;
fftLenZeta = 1024;
fftLenXi = 64;
fftLen1 = 32;
fftLen2 = 32;
fftLenGamma = 32; % >=Mr
numTrainingSet = 4;
numTrainingSetLen = length(numTrainingSet);
U = 1; 

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
zcLen = 600;   % 截取低频部分的信号长度
lTilde = 1024; % ZC原长度-lTilde/2至lTilde/2-1，偶数
cpLen = 64;    % 大于时延 
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
TFft = reshape(TFft2,fftLen1*fftLen2,numTraining*Mr); 
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
        Y1 = Y1Ori(:,end/2+(-zcLen/2+1:zcLen/2)); 
        YTilde = Y1*diag(conj(zcSeq));

        %% SAGE      
         sigPerPath = zeros(Mr*numTraining,zcLen,U);
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
             end
             
             YTildeSageTmp(gg) = norm(sum(sigPerPath,3)-YTilde, 'fro');
             %%% SAGE收敛条件
             if YTildeSageTmp(gg)./norm(sum(sigPerPath,3),'fro')<1e-8
                 break;
             end
                     
         end

        %% 参数rmse
         mseTauTmp(mm) = (tauNewtonPerPath - tauList)^2 ;
         mseTau(tt,ss) = mseTau(tt,ss) + mseTauTmp(mm);
         mseXiTmp(mm) = (xiNewtonPerPath - xiList)^2 ;
         mseXi(tt,ss) = mseXi(tt,ss) + mseXiTmp(mm);
         mseThetaTmp(mm) = ((thetaNewtonPerPath - thetaList)/ pi * 180)^2 ;
         mseTheta(tt,ss) = mseTheta(tt,ss) + mseThetaTmp(mm);
         msePhiTmp(mm) = ((phiNewtonPerPath - phiList(1))/ pi * 180)^2 ;
         msePhi(tt,ss) = msePhi(tt,ss) + msePhiTmp(mm);
    end
end

end

%% CRB
[crbTau,crbXi,crbTheta,crbPhi,crbBetaR,crbBetaI,crbBeta] = crb_compute_multipath(U,snrDbSet,numTraining,roCoeff,Tp,dx,dz,P,Q,Mr,lTilde,zcLen,zcCpLen,W,zcRoot1,zcRootOri,aRIS,tauList,xiList,thetaList,phiList,betaList);
figure;
subplot(2,2,1)
semilogy(snrDbSet,sqrt(crbTau(1,:)));hold on
semilogy(snrDbSet,sqrt(mseTau(1,:)/numMonte),'-o');
xlabel('SNR(dB)');
ylabel('RMSE(Ts)')
legend('CRB','Our Algorithm')
title(['\tau(Ts) = ',num2str(tauList(1))])
grid on;

subplot(2,2,2)
semilogy(snrDbSet,sqrt(crbXi(1,:)));hold on
semilogy(snrDbSet,sqrt(mseXi(1,:)/numMonte),'-o')
xlabel('SNR(dB)');
ylabel('RMSE(1/Ts)')
legend('CRB','Our Algorithm')
title(['\xi(1/Ts) = ',num2str(xiList(1))])
grid on;

subplot(2,2,3)
semilogy(snrDbSet,sqrt(crbTheta(1,:))/pi*180);hold on
semilogy(snrDbSet,sqrt(mseTheta(1,:)/numMonte),'-o')
xlabel('SNR(dB)');
ylabel('RMSE(\circ)')
legend('CRB','Our Algorithm')
title(['\theta(\circ) = ',num2str(thetaList(1)/pi*180)])
grid on;

subplot(2,2,4)
semilogy(snrDbSet,sqrt(crbPhi(1,:))/pi*180);hold on
semilogy(snrDbSet,sqrt(msePhi(1,:)/numMonte),'-o')
xlabel('SNR(dB)');
ylabel('RMSE(\circ)')
legend('CRB','Our Algorithm')
title(['\phi(\circ) = ',num2str(phiList(1)/pi*180)])
grid on;
