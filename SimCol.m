function Uc=SimCol(U,c)
% from the code of Anh-Huy Phan (Tensorbox)

R=size(U,2);
K=ones(R,R)*c;
for k=1:R
    K(k,k)=1;
end
C=chol(K);

[Q,~] = qr(U,0);    
Uc=Q*C;