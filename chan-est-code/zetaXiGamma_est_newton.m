function [zetaNewton,xiNewton,gammaNewton,cNewton,BXi] = zetaXiGamma_est_newton(YTilde,Mr,numTraining,zcLen,zcCpLen,zetaNewton,xiNewton,gammaNewton,maxNewtonIter,minNewtonStepSize,dg)

%% 牛顿迭代优化
zetaRange = -1j*2*pi*(-zcLen/2:zcLen/2-1).';
xiRange = -1j*2*pi*zcCpLen*(0:numTraining-1).';
gammaRange = -1j*2*dg*pi*(0:Mr-1).';
for nn = 1:maxNewtonIter
%% 参数函数求导
    %%% gamma
    cNewton = exp(gammaRange*cos(gammaNewton)); 
    gammaPartd = -gammaRange*sin(gammaNewton);
    gammaPartd2 = -gammaRange*cos(gammaNewton);
    cGammad1 = cNewton.*gammaPartd;
    cGammad2 = cNewton.*(gammaPartd.*gammaPartd+gammaPartd2);% cGammad1.*gammaPartd;    
    %%% zeta      
    dZeta = exp(zetaRange*zetaNewton);
    dZetad1 = dZeta.*zetaRange;
    dZetad2 = dZetad1.*zetaRange;
    %%% xi
    pXi = exp(xiRange*xiNewton);
    pXid1 = pXi.*xiRange;  
    pXid2 = pXid1.*xiRange;
    
%% 目标函数求导(无分母)
    BXi = kron(exp(1j*2*pi*xiNewton*(0:numTraining-1).'*zcCpLen),eye(Mr));
    BXid1 = kron(exp(1j*2*pi*xiNewton*(0:numTraining-1).'*zcCpLen).*(1j*2*pi*(0:numTraining-1).'*zcCpLen),eye(Mr));

    %%% gamma
    rGamma = BXi'*YTilde*dZeta;
    RGamma = rGamma*rGamma';    
    fGammad1 = 2*real(cNewton'*RGamma*cGammad1); 
    fGammad2 = 2*real(cNewton'*RGamma*cGammad2) + 2*norm(rGamma'*cGammad1)^2;

    rrGammaZeta = BXi'*YTilde*dZetad1;
    rrGammaXi = BXid1'*YTilde*dZeta;
    RRGammaZeta = rrGammaZeta*rGamma';
    RRGammaXi = rrGammaXi*rGamma';
    fGammaZetad2 = 2*real(cGammad1'*RRGammaZeta*cNewton) + 2*real(cNewton'*RRGammaZeta*cGammad1);
	fGammaXid2 = 2*real(cGammad1'*RRGammaXi*cNewton) + 2*real(cNewton'*RRGammaXi*cGammad1);
    
    %%% zeta
    rZeta = cNewton'*BXi'*YTilde;
    RZeta = rZeta'*rZeta;
    fZetad1 = 2*real(dZeta'*RZeta*dZetad1);
    fZetad2 = 2*real(dZeta'*RZeta*dZetad2) + 2*abs(rZeta*dZetad1)^2;   

    rrZeta = cNewton'*BXid1'*YTilde;   
    RRZeta = rZeta'*rrZeta;
    fZetaXid2 = 2*real(dZetad1'*RRZeta*dZeta) + 2*real(dZeta'*RRZeta*dZetad1);
    
    %%% xi
    rXi = zeros(1,numTraining);
    for kk = 1 : numTraining
        dTmp = YTilde((kk-1)*Mr+(1:Mr),:)*dZeta;
        rXi(kk) = cNewton'*dTmp;
        ddTmp = YTilde((kk-1)*Mr+(1:Mr),:)*dZetad1;
        rrXi(kk) = cNewton'*ddTmp;
    end
    RXi = rXi'*rXi;
    RRXi = rXi'*rrXi;
    fXid1 = 2*real(pXi'*RXi*pXid1);
    fXid2 = 2*real(pXi'*RXi*pXid2) + 2*abs(rXi*pXid1)^2;
        
    fZetaGammad2 = fGammaZetad2;
    fXiGammad2 = fGammaXid2;
    fXiZetad2 = fZetaXid2;
    
    f = norm(rGamma'*cNewton)^2;
    g = [fGammad1;fZetad1;fXid1];
    H = [fGammad2,fGammaZetad2,fGammaXid2;...
        fZetaGammad2,fZetad2,fZetaXid2;...
        fXiGammad2,fXiZetad2,fXid2;];

%% 回溯直线
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
        gammaNewtonTmp = gammaNewton - step*ss(1);
        zetaNewtonTmp = zetaNewton - step*ss(2);
        xiNewtonTmp = xiNewton - step*ss(3); 
        cTmp = exp(-1j*2*dg*pi*(0:Mr-1).'*cos(gammaNewtonTmp)); 
        BXiTmp = kron(exp(1j*2*pi*xiNewtonTmp*(0:numTraining-1).'*zcCpLen),eye(Mr));
        dZetaTmp = exp(zetaRange*zetaNewtonTmp);
        rTmp = BXiTmp'*YTilde*dZetaTmp;     
        fTmp(kk) = norm(rTmp'*cTmp)^2; 
        if fTmp(kk) < f - alpha*step*jacMatrix.'*ss
            step = roughDegree*step;
        else
            break
        end
        if step < 1e-5
            break
        end
    end
    gammaNewton = gammaNewton - step*ss(1);  
    zetaNewton = zetaNewton - step*ss(2);
    xiNewton = xiNewton - step*ss(3);  
    cNewton = exp(-1j*2*dg*pi*(0:Mr-1).'*cos(gammaNewton)); 
    BXi = kron(exp(1j*2*pi*xiNewton*(0:numTraining-1).'*zcCpLen),eye(Mr));
    
    if norm(step*ss) < minNewtonStepSize
        break
    end
    
end

end

