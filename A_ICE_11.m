function [A]=A_ICE_11(N,Yh,Ym,S,A_old,D,g,B,iterA,lambda1,tildeeta,PtP,mode)
M = size(Yh,1);
z=reshape(A_old,[],1);
vtilde = zeros(M*N,1);
e=speye(M*N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_3d=reshape(S,N,size(g,1),[]);
S_perm=permute(S_3d,[1,3,2]);
S_re=reshape(S_perm,[],size(g,1));
Sgv=reshape(S_re*g,N,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C2ty=reshape(Yh*Sgv',[],1)+reshape(D'*Ym*S',[],1);
for j=0:iterA
    a = inv(lambda1*PtP+tildeeta*e+kron((Sgv*Sgv'),eye(M)) ...
        +kron((S*S'),(D'*D)))*(C2ty+tildeeta*z-vtilde);
    z = max(a+vtilde/tildeeta,0);
    vtilde = vtilde + tildeeta*(a-z);
end
A=reshape(z,M,N);
return;