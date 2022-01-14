function [H] = channel_generation_RIS2BS(M,Mr,d0,mu,m,lambda,dg,dx)
% ***************************************
%  generate channel from RIS to BS
%  author - Xuemeng Zhou
%  input- M: RIS element number
%         Mr: BS element number
%         d0: distance between RIS and BS
%         mu: reflection efficiency
%         m: the row of RIS
%         lambda£ºwavelength
%         dg£ºinter-element distance for receiver
%         dx£ºinter-element distance for transmitter
%  output-H: channel between BS to RIS
%
%  copyright - CSRL@Fudan,2021
%  ************************************
n = M/m; 
H = zeros(M,Mr);
for ii = 1:Mr
    x1 = -(Mr-1)/2*dg*lambda+(ii-1)*dg*lambda;
    jj = 1;
    for mm=1:m
        for nn=1:n
            x2=-(n-1)/2*dx*lambda+(nn-1)*dx*lambda;
            y2=(m-1)/2*dx*lambda-(mm-1)*dx*lambda;
            d=sqrt((x2-x1)^2+y2^2+d0^2);
            H(jj,ii)=exp(1i*2*pi*d/lambda)*sqrt(mu/M);
            jj = jj+1;
        end
    end
end

end