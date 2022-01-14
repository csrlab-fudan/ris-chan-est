function A_RIS = generate_ar(dx, dz, P, Q, theta, phi)
A_RIS = zeros(P, Q);
for p = 1 : P
    for q = 1 : Q
        A_RIS(p,q) = exp(-1j*2*pi*(q-1)*dx*sin(phi)*cos(theta)+1j*2*pi*(p-1)*dz*cos(phi));
    end
end
end
% B = blkdiag(A1,...,AN) ����ͨ���� B �ĶԽ��߶���������� A1,...,AN �����ķֿ�ԽǾ���