function [zetaNewton,xiNewton,WXi] = zetaXi_est_newton(W,YTilde,Mr,numTraining,zcLen,zcCpLen,zetaNewton,xiNewton,aRISNewton,maxNewtonIter,minNewtonStepSize)

%%
zetaRange = -1j*2*pi*(-zcLen/2:zcLen/2-1).';
xiRange = -1j*2*pi*zcCpLen*(0:numTraining-1).';
for nn = 1:maxNewtonIter
    for tt = 1 : numTraining
        WXi((tt-1)*Mr+(1:Mr),:) = W((tt-1)*Mr+(1:Mr),:)*exp(1j*2*pi*xiNewton*(tt-1)*zcCpLen);
        WXid1((tt-1)*Mr+(1:Mr),:) = WXi((tt-1)*Mr+(1:Mr),:)*(1j*2*pi*(tt-1)*zcCpLen);
    end
       
    dZeta = exp(zetaRange*zetaNewton);
    dZetad1 = dZeta.*zetaRange;
    dZetad2 = dZetad1.*zetaRange;
    rZeta = aRISNewton'*WXi'*YTilde;
    rrZeta = aRISNewton'*WXid1'*YTilde;
    RZeta = rZeta'*rZeta;
    RRZeta = rZeta'*rrZeta;
    fZetad1 = 2*real(dZeta'*RZeta*dZetad1);
    fZetad2 = 2*real(dZeta'*RZeta*dZetad2) + 2*abs(rZeta*dZetad1)^2;
    
    rXi = zeros(1,numTraining); 
    for kk = 1 : numTraining
        bTmp = W((kk-1)*Mr+(1:Mr),:)*aRISNewton;
        dTmp = YTilde((kk-1)*Mr+(1:Mr),:)*dZeta;
        rXi(kk) = bTmp'*dTmp;
        ddTmp = YTilde((kk-1)*Mr+(1:Mr),:)*dZetad1;
        rrXi(kk) = bTmp'*ddTmp;
    end
    RXi = rXi'*rXi;
    RRXi = rXi'*rrXi;
    pXi = exp(xiRange*xiNewton);
    pXid1 = pXi.*xiRange;
    pXid2 = pXi.*xiRange.*xiRange;
    fXid1 = 2*real(pXi'*RXi*pXid1);
    fXid2 = 2*real(pXi'*RXi*pXid2) + 2*abs(rXi*pXid1)^2;  
%     fZetaXid2 = 2*real(dZeta'*YTilde'*(2*real(WXid1*(aRISNewton*aRISNewton')*WXi'))*YTilde*dZetad1);
    fZetaXid2 = 2*real(dZetad1'*RRZeta*dZeta) + 2*real(dZeta'*RRZeta*dZetad1);
    fXiZetad2 = fZetaXid2; % 2*real(pXid1'*RRXi*pXi) + 2*real(pXi'*RRXi*pXid1);
    
%     f = abs(aRISNewton'*WXi'*YTilde*dZeta)^2;
    f = norm(rZeta*dZeta)^2;
    g = [fZetad1;fXid1];
    H = [fZetad2,fZetaXid2;fXiZetad2,fXid2];

    %% 回溯直线
    alpha = 0.1;
    roughDegree = 0.9;
    Loop = 100;
    jacMatrix = g;
    hesMatrix = H;
    [V,E] = eig(hesMatrix); % V特征向量E特征值;特征值分解，若牛顿迭代的hessian矩阵是不定的（特征值有负数），则应将负数特征值置零
    e = diag(E);
    V1 = V(:,e<0); % 最大化问题；将V中小于0的数给V1
    e1 = e(e<0);
    ss = V1*((V1'*jacMatrix)./e1); % 特征值分解求逆ss = hesMatrix\jacMatrix;
    step= 1; 
    fTmp = [];
    for kk = 1:Loop
        zetaNewtonTmp = zetaNewton - step*ss(1);
        xiNewtonTmp = xiNewton - step*ss(2);
        for tt = 1 : numTraining
            WXiTmp((tt-1)*Mr+(1:Mr),:) = W((tt-1)*Mr+(1:Mr),:)*exp(1j*2*pi*xiNewtonTmp*(tt-1)*zcCpLen);
        end
        dZetaTmp = exp(zetaRange*zetaNewtonTmp);
        fTmp(kk) = abs(aRISNewton'*WXiTmp'*YTilde*dZetaTmp)^2;
        if fTmp(kk) < f - alpha*step*jacMatrix.'*ss
            step = roughDegree*step;
        else
            break
        end
        if step < 1e-5
            break
        end
    end
    zetaNewton = zetaNewton - step*ss(1);
    xiNewton = xiNewton - step*ss(2); 
   
    if norm(step*ss) < minNewtonStepSize
        break
    end
    
end

end

