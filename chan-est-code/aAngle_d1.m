function [aThetad1,aPhid1] = aAngle_d1(dx, dz, P, Q, theta, phi)
% first derivative of {theta,phi}

rRange = -1j*2*pi*dx*(0:Q-1);
cRange = 1j*2*pi*dz*(0:P-1)';
Ar = generate_ar(dx, dz, P, Q, theta, phi);
% a = reshape(Ar.',[],1);
thetaPartd = -ones(P,1)*rRange*sin(phi)*sin(theta);
phiPartd = ones(P,1)*rRange*cos(phi)*cos(theta)-cRange*ones(1,Q)*sin(phi);
ArThetad1 = Ar.*thetaPartd;
ArPhid1 = Ar.*phiPartd;
aThetad1 = reshape(ArThetad1.',[],1);
aPhid1 = reshape(ArPhid1.',[],1);
end

