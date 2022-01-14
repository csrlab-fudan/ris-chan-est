function [ response ] = rrc_d1( tVec ,alpha )
% first derivative of rrc
% a is roll factor beta

b = 1./(1-4*alpha^2*tVec.^2);
c = (cos(pi*tVec)-sin(pi*tVec)./(pi*tVec))./tVec;
d = -pi*alpha*b.*sin(pi*alpha*tVec) + 8*alpha^2*tVec.*b.^2.*cos(pi*alpha*tVec);
diff = c.*cos(pi*alpha*tVec).*b + sin(pi*tVec)./(pi*tVec).*d;

t1 = find(tVec==0);
if isempty(t1)~=1 % 确定数组是否为空，空则返回逻辑值1（true）
    diff(t1) = 0;
end
t2 = find(tVec==1/(2*alpha));
if isempty(t2)~=1
    diff(t2) = pi*alpha*cos(pi/(2*alpha))/2;
end
t3 = find(tVec==-1/(2*alpha));
if isempty(t3)~=1
    diff(t3) = -pi*alpha*cos(pi/(2*alpha))/2;
end
response = diff.';
% hold on;plot(response,'b')
end


