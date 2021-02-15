% NOTE: these values are based on 'sorted' eigenvalues, eigenvectors which
% are then truncated. To compare against complex_analysis.m or simran.m
% this sorting needs to be the same.

vl = phi_left';
vr = phi_right';
len=size(vl,1);
zl=zeros(len,len);
zr=zeros(len,len);
for i=1:len
  zl(:,i) = vl(:,2*i-1) + j.*vl(:,2*i);
  zr(:,i) = vr(:,2*i-1) + j.*vr(:,2*i);
end
a = Dggev_A';
b = Dggev_B';

Ra = zl'*a*zr;
RRa = abs(Ra);
Rb = zl'*b*zr;
RRb = abs(Rb);
format long g;
eps=1e-5;
spy(abs(RRb)>eps);
RRb
RRa

