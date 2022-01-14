function [crbTau,crbXi,crbTheta,crbPhi,crbBetaR,crbBetaI,crbBeta] = crb_compute_multipath(U,snrDbSet,numTraining,roCoeff,Tp,dx,dz,P,Q,Mr,lTilde,zcLen,zcCpLen,W,zcRoot1,zcRootOri,aRIS,tauList,xiList,thetaList,phiList,betaList)

%% 预设
crbTau = zeros(U,length(snrDbSet));
crbXi = zeros(U,length(snrDbSet));
crbTheta = zeros(U,length(snrDbSet));
crbPhi = zeros(U,length(snrDbSet));
paramEstNum = 6; % 待估参数个数
y = zeros(numTraining*Mr*zcLen,paramEstNum*U); 

%% 求导
for uu = 1:U
    tTmp = (-Tp-tauList(uu)):(Tp-tauList(uu)); 
    zcLcTmp = zcRoot1(uu,end/2+(-zcLen/2+1:zcLen/2));
    zcCpTp = conv(-rrc_d1(tTmp,roCoeff),zcRootOri);
    zcCp = zcCpTp(end/2+(-zcCpLen/2+1:zcCpLen/2));
    zcTmp = zcCp(zcCpLen-lTilde+1:end);
    zcLcTmpd1 = zcTmp(end/2+(-zcLen/2+1:zcLen/2));
    dXiTmp = exp(1j*2*pi*xiList(uu)*(-zcLen/2:zcLen/2-1));
    dXiTmpd1 = dXiTmp.*(1j*2*pi*(-zcLen/2:zcLen/2-1));
    for tt = 1 : numTraining
        WXiTmp((tt-1)*Mr+(1:Mr),:) = W((tt-1)*Mr+(1:Mr),:)*exp(1j*2*pi*xiList(uu)*(tt-1)*zcCpLen);
        WXiTmpd1((tt-1)*Mr+(1:Mr),:) = WXiTmp((tt-1)*Mr+(1:Mr),:)*(1j*2*pi*(tt-1)*zcCpLen);
    end
    [aThetaTmpd1,aPhiTmpd1] = aAngle_d1( dx, dz, P, Q, thetaList(uu), phiList(uu) );
    
    yTauTmpd1 = betaList(uu)*kron( (zcLcTmpd1.*dXiTmp).',(WXiTmp*aRIS(:,uu)) );
    yXiTmpd1 = betaList(uu)*kron( (zcLcTmp.*dXiTmpd1).',(WXiTmp*aRIS(:,uu)) )+...
        betaList(uu)*kron( (zcLcTmp.*dXiTmp).',(WXiTmpd1*aRIS(:,uu)) );
    yThetaTmpd1 = betaList(uu)*kron( (zcLcTmp.*dXiTmp).',(WXiTmp*aThetaTmpd1) );
    yPhiTmpd1 = betaList(uu)*kron( (zcLcTmp.*dXiTmp).',(WXiTmp*aPhiTmpd1) );
    yBetaRTmpd1 = kron( (zcLcTmp.*dXiTmp).',(WXiTmp*aRIS(:,uu)) );
    yBetaITmpd1 = 1j*yBetaRTmpd1;
    y(:,(uu-1)*paramEstNum+(1:paramEstNum)) = [yTauTmpd1,yXiTmpd1,yThetaTmpd1,yPhiTmpd1,yBetaRTmpd1,yBetaITmpd1];
end
rx = y; % 全频段时

%% fisher 矩阵
Fisher = 2*real(rx'*rx);
FisherInv = inv(Fisher);
for uu = 1:U
    crbTau(uu,:) = FisherInv((uu-1)*paramEstNum+1,(uu-1)*paramEstNum+1).*db2pow(-snrDbSet);
    crbXi(uu,:) = FisherInv((uu-1)*paramEstNum+2,(uu-1)*paramEstNum+2).*db2pow(-snrDbSet);  
    crbTheta(uu,:) = FisherInv((uu-1)*paramEstNum+3,(uu-1)*paramEstNum+3).*db2pow(-snrDbSet);    
    crbPhi(uu,:) = FisherInv((uu-1)*paramEstNum+4,(uu-1)*paramEstNum+4).*db2pow(-snrDbSet);
    crbBetaR(uu,:) = FisherInv((uu-1)*paramEstNum+5,(uu-1)*paramEstNum+5).*db2pow(-snrDbSet);
    crbBetaI(uu,:) = FisherInv((uu-1)*paramEstNum+6,(uu-1)*paramEstNum+6).*db2pow(-snrDbSet);
    crbBeta(uu,:) = crbBetaR(uu,:).^2 + crbBetaI(uu,:).^2;
end

end



