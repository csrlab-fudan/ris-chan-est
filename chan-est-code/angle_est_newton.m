function [thetaNewton,phiNewton,aRISNewton] = angle_est_newton(zetaNewton,WXi,W,zcLen,YTilde,P,Q,dx,dz,thetaFft,phiFft,maxNewtonIter,minNewtonStepSize)

dZeta = exp(-1j*2*pi*zetaNewton*(-zcLen/2:zcLen/2-1).');
RA = W'*W; 
r = WXi'*YTilde*dZeta;
RB = r*r';
for nn = 1:maxNewtonIter
%% URA阵列及其导数列向量化       
    rRange = -1j*2*pi*dx*(0:Q-1);
    cRange = 1j*2*pi*dz*(0:P-1)';
    Ar = generate_ar(dx, dz, P, Q, thetaFft, phiFft);
    a = reshape(Ar.',[],1); 
    
    thetaPartd = -ones(P,1)*rRange*sin(phiFft)*sin(thetaFft);
    phiPartd = ones(P,1)*rRange*cos(phiFft)*cos(thetaFft)-cRange*ones(1,Q)*sin(phiFft);
    ArThetad1 = Ar.*thetaPartd;
    ArPhid1 = Ar.*phiPartd;
    ArThetaThetad2 = ArThetad1.*thetaPartd - Ar.*(ones(P,1)*rRange*sin(phiFft)*cos(thetaFft));
    ArPhiPhid2 = ArPhid1.*phiPartd - Ar.*(ones(P,1)*rRange*sin(phiFft)*cos(thetaFft)+cRange*ones(1,Q)*cos(phiFft));
    ArThetaPhid2 = ArPhid1.*thetaPartd - Ar.*(ones(P,1)*rRange*cos(phiFft)*sin(thetaFft));
    
    aThetad1 = reshape(ArThetad1.',[],1); 
    aPhid1 = reshape(ArPhid1.',[],1); 
    aThetaThetad2 = reshape(ArThetaThetad2.',[],1); 
    aPhiPhid2 = reshape(ArPhiPhid2.',[],1); 
    aThetaPhid2 = reshape(ArThetaPhid2.',[],1); 

%% 目标函数求导
    % 分子及其各阶偏导
    BZC1 = norm(r'*a)^2; 
    BThetaZC1d1 = 2*real(a'*RB*aThetad1);
    BPhiZC1d1 = 2*real(a'*RB*aPhid1);
    BThetaZC1d2 = 2*real(a'*RB*aThetaThetad2) + 2*norm(r'*aThetad1)^2;
    BPhiZC1d2 = 2*real(a'*RB*aPhiPhid2) + 2*norm(r'*aPhid1)^2;
    BThetaPhiZC1d2 = 2*real(a'*RB*aThetaPhid2) + 2*real(aThetad1'*RB*aPhid1);
    % 分母及其各阶偏导
    AZC1 = norm(W*a)^2; 
    AThetaZC1d1 = 2*real(a'*RA*aThetad1);
    APhiZC1d1 = 2*real(a'*RA*aPhid1);
    AThetaZC1d2 = 2*real(a'*RA*aThetaThetad2) + 2*norm(W*aThetad1)^2;
    APhiZC1d2 = 2*real(a'*RA*aPhiPhid2) + 2*norm(W*aPhid1)^2;
    AThetaPhiZC1d2 = 2*real(a'*RA*aThetaPhid2) + 2*real(aThetad1'*RA*aPhid1);
    % 目标函数各阶偏导
    fThetaZC1d1 = (BThetaZC1d1*AZC1-BZC1*AThetaZC1d1)/(AZC1'*AZC1);
    fPhiZC1d1 = (BPhiZC1d1*AZC1-BZC1*APhiZC1d1)/(AZC1'*AZC1);
    fThetaZC1d2 = (BThetaZC1d2*AZC1 - BZC1*AThetaZC1d2 - 2*fThetaZC1d1*AZC1*AThetaZC1d1)/(AZC1'*AZC1);
    fPhiZC1d2 = (BPhiZC1d2*AZC1 - BZC1*APhiZC1d2 - 2*fPhiZC1d1*AZC1*APhiZC1d1)/(AZC1'*AZC1);
    fThetaPhiZC1d2 = (BThetaPhiZC1d2*AZC1 + BThetaZC1d1*APhiZC1d1 - BPhiZC1d1*AThetaZC1d1 - BZC1*AThetaPhiZC1d2 - 2*fThetaZC1d1*AZC1*APhiZC1d1)/(AZC1'*AZC1);

    f = BZC1/AZC1;
    g = [fThetaZC1d1;fPhiZC1d1];
    H = [fThetaZC1d2,fThetaPhiZC1d2;fThetaPhiZC1d2,fPhiZC1d2];

%% 回溯直线优化步长的牛顿迭代
    alpha = 0.1;
    roughDegree = 0.9;
    Loop = 100;
    jacMatrix = g;
    hesMatrix = H;
    [V,E] = eig(hesMatrix); % V特征向量E特征值;特征值分解，若牛顿迭代的hessian矩阵是不定的（特征值有负数），则应将负数特征值置零
    e = diag(E);
    V1 = V(:,e<0); % 将V中小于0的数给V1
    e1 = e(e<0);
    ss = V1*((V1'*jacMatrix)./e1); % 特征值分解求逆ss = hesMatrix\jacMatrix;
    step= 1; 
    fTmp = [];
    for kk = 1:Loop
        thetaFftTmp = thetaFft - step*ss(1);
        phiFftTmp = phiFft - step*ss(2);
        ArTmp = generate_ar(dx, dz, P, Q, thetaFftTmp, phiFftTmp);
        aTmp = reshape(ArTmp.',[],1);
        BTmp = norm(r'*aTmp)^2;
        ATmp = norm(W*aTmp)^2;       
        fTmp(kk) = BTmp/ATmp;
        if fTmp(kk) < f - alpha*step*jacMatrix.'*ss
            step = roughDegree*step;
        else
            break
        end
        if step < 1e-5
            break
        end
    end
    thetaFft = thetaFft - step*ss(1);
    phiFft = phiFft - step*ss(2);  
        
    if norm(step*ss) < minNewtonStepSize
        break
    end
    
end

thetaNewton = thetaFft; 
phiNewton = phiFft;
ARISNewton = generate_ar(dx, dz, P, Q, thetaNewton, phiNewton);
aRISNewton = reshape(ARISNewton.',[],1); 

end

