function [S]=S_NonNeg_Lasso_11(N,Yh,Ym,A,S_old,D,g,B,iterS,lambda2,eta,mode)
L = size(Ym, 2);
Lh =size(Yh, 2);
r2 = L / Lh;
ym = Ym(:);
yh = Yh(:);
x=reshape(S_old,N*r2,Lh);%%%%%%%%%%%%%%%%%%%%%%改dimension
%%%%%%%%%%%%%%%v = zeros(N * r2, Lh);%%%%%%%%%%%%%%%%%%%%%%%%%%改dimension
u = zeros(N * r2, Lh);
K = kron(g', A);
KER2DAT = (kron(eye(r2), D * A))';
C_bar =( [ K', KER2DAT ] )';
P =inv ((C_bar)' * C_bar + eta * eye(r2*N));%%inv(CT * C + eta * Ir2N)
%%%%%%%%%%%%%%%s = zeros(Lh*r2*N,1);
a=r2*N;
m=P * K';
M=size(m,2);
nu =P * KER2DAT;
Bandr2=size(nu,2);
ONES = ones(N * r2, Lh);%%%%%%%%%%%%%%%%%改dimension
XU = reshape((x-u), N*r2, Lh);%%%%%%%%%%%%%%%%%新增

YH = reshape(yh, M, Lh);%%%%%%%%%%%%%%%%%%%%%新增
YM = reshape(ym, Bandr2, Lh);%%%%%%%%%%%%%%%%%%%%%%新增
for j = 1:iterS
   s_new = P * XU + m * YH + nu * YM;%%490*900
   x = max(s_new + u - (lambda2/eta) * ONES,0);
   u = u + s_new - x;
   XU=x-u;
end
S=reshape(x,N,L);
return;