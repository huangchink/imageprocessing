%======================================================================
%  Input
%  X is M-by-L data matrix, where M is the number of spectral bands and L is the number of pixels.
%  N is the number of endmembers.
%----------------------------------------------------------------------
%  Output
%  A_est is M-by-N mixing matrix whose columns are estimated endmember signatures.
%  S_est is N-by-L source matrix whose rows are estimated abundance maps.
%  time is the computation time (in secs).
%========================================================================

function [A_est, S_est, time] = HyperCSI(X,N)
t0 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_of_columns = sum(X, 2);                  %%% sum_of_columns[224, 1]
average_of_columns = sum_of_columns./10000;
U = zeros(224, 10000);
for i=1:1:10000
    U(:, i) = X(:, i) - average_of_columns;  %%% U = X - d, U[224, 10000]
end
UUT = U * transpose(U);                      %%% UUT[224, 224]
[C, ~] = eigs(UUT, 224);                     %%% C ¬°eigenvector, C[224, 2]
C=[C(:,2),C(:,1)];                           
X_bar = transpose(C) * U;                    %%% X_bar[2*10000]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DR=[X_bar;ones(1,10000)];                    %%% DR  [3*10000]
purest_pixel=zeros(2,3);

for j=1:3
    %%%%%%%%%%% Maximization
    Norm=zeros(1,10000);
    for i=1:10000
        Norm(:,i)=norm(DR(:,i));           %%%[1*10000]
    end

    [Maxvalue,li]=max(Norm);
    purest_pixel(:,j) = X_bar(:,li);
    %%%%%%%%%%% Projection
    p= eye(3)-(1/(Maxvalue)^2)*DR(:,li)*transpose(DR(:,li));
    DR = p*DR;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find b~ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_tilde = zeros(2,3);

i=[2 3 1 2];
for j=1:3
    
    p=purest_pixel(:,i(j))-purest_pixel(:,i(j+1));
    b_tilde(:,j) =( eye(2) - p / (transpose(p)*p) * transpose(p) ) * (purest_pixel(:,i(j+1))-purest_pixel(:,j));
    b_tilde(:,j) =  b_tilde(:,j)/norm(b_tilde(:,j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find active pixel %%%%%%%%%%%%%%%%%%%%%%
active_pixel=[];

                      %%%%%%%%% alpha1 %%%%%%%%%% 
r = norm(purest_pixel(:,2)-purest_pixel(:,3))/2;
R=zeros(2,10000);
distance =zeros(1,10000);

for i=1:10000
   if (norm(X_bar(:,i)-purest_pixel(:,2)) < r )
       R(:,i)=X_bar(:,i);
       distance(:,i)=transpose(b_tilde(:,2))*X_bar(:,i);
   end
end
[~,l]=max(distance);
active_pixel=[active_pixel,R(:,l)];
R=zeros(2,10000);
distance =zeros(1,10000);
for i=1:10000
   if (norm(X_bar(:,i)-purest_pixel(:,3)) < r )
        R(:,i)=X_bar(:,i);
       distance(:,i)=transpose(b_tilde(:,2))*X_bar(:,i);
   end
end
[~,l]=max(distance);
active_pixel=[active_pixel,R(:,l)];

               %%%%%%%%%%%% alpha2 %%%%%%%%%%%%%% 
r = norm(purest_pixel(:,1)-purest_pixel(:,3))/2;
R=zeros(2,10000);
distance =zeros(1,10000);
for i=1:10000
   if (norm(X_bar(:,i)-purest_pixel(:,1)) < r )
       R(:,i)=X_bar(:,i);
       distance(:,i)=transpose(b_tilde(:,2))*X_bar(:,i);
   end
end
[~,l]=max(distance);
active_pixel=[active_pixel,R(:,l)];
R=zeros(2,10000);
distance =zeros(1,10000);
for i=1:10000
   if (norm(X_bar(:,i)-purest_pixel(:,3)) < r )
       R(:,i)=X_bar(:,i);
       distance(:,i)=transpose(b_tilde(:,2))*X_bar(:,i);
   end
end
[~,l]=max(distance);
active_pixel=[active_pixel,R(:,l)];

              %%%%%%%%%%%%%%% alpha3 %%%%%%%%%%%%%%%%%% 
r = norm(purest_pixel(:,1)-purest_pixel(:,2))/2;
R=zeros(2,10000);
distance =zeros(1,10000);

for i=1:10000
   if (norm(X_bar(:,i)-purest_pixel(:,1)) < r )
       R(:,i)=X_bar(:,i);
       distance(:,i)=transpose(b_tilde(:,2))*X_bar(:,i);
   end
end
[~,l]=max(distance);
active_pixel=[active_pixel,R(:,l)];
R=zeros(2,10000);
distance =zeros(1,10000);

for i=1:10000
   if (norm(X_bar(:,i)-purest_pixel(:,2)) < r )
       R(:,i)=X_bar(:,i);
       distance(:,i)=transpose(b_tilde(:,2))*X_bar(:,i);
   end
end
[~,l]=max(distance);
active_pixel=[active_pixel,R(:,l)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find b_hat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b_hat=zeros(2,3);
COS=dot((active_pixel(:,1)-active_pixel(:,2)),b_tilde(:,1))/(norm(active_pixel(:,1)-active_pixel(:,2))*norm(b_tilde(:,1)));
b_hat(:,1)=(sqrt(1-COS^2))*b_tilde(:,1);

COS=dot((active_pixel(:,3)-active_pixel(:,4)),b_tilde(:,2))/(norm(active_pixel(:,3)-active_pixel(:,4))*norm(b_tilde(:,2)));
b_hat(:,2)=(sqrt(1-COS^2))*b_tilde(:,2);

COS=dot((active_pixel(:,5)-active_pixel(:,6)),b_tilde(:,3))/(norm(active_pixel(:,5)-active_pixel(:,6))*norm(b_tilde(:,3)));
b_hat(:,3)=(sqrt(1-COS^2))*b_tilde(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find h_hat   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_hat=zeros(1,3);
for i=1:3
    h = zeros(1, 10000);
    for j=1:10000
        h(:, j) = transpose(b_hat(:,i))*X_bar(:,j);
    end
    h_hat(:,i)=max(h);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find c   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1 = transpose(b_hat);
B1(1,:) =[];
B2 = transpose(b_hat);
B2(2,:) =[];
B3 = transpose(b_hat);
B3(3,:) =[];

h1 = transpose(h_hat);
h1(1,:) =[];
h2 = transpose(h_hat);
h2(2,:) =[];
h3 = transpose(h_hat);
h3(3,:) =[];

c_bar = [C*(transpose(B1)*h1),C*(transpose(B2)*h2),C*(transpose(B3)*h3)];
for i=1:224
    c_bar(i,:) = -1.*c_bar(i,:)./average_of_columns(i);
end
c = max(1,max(max(c_bar)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  alpha_hat A_est  %%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_hat = [B1 \ (h1./c),B2 \ (h2./c),B3 \ (h3./c)];

A_est = C*alpha_hat + [average_of_columns,average_of_columns,average_of_columns];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  S_est  %%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_est=ones(N,10000);
for i=1:3
    for n=1:10000
        S_est(i,n)=max((h_hat(i)-transpose(b_hat(:,i))*X_bar(:,n))/(h_hat(i)-transpose(b_hat(:,i))*alpha_hat(:,i)),0);
    end
end
time = etime(clock,t0);