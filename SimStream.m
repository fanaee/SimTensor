function [X0,U1]=SimStream(U0,I,R,v)

% ***********************************
% The idea is borrowed from parafac_sdt
% by Dimitri Nion and Nikos D. Sidiropoulos
%
% U0: factor matrices (no temporal mode)
% I: size of factor
% R: number of components
% v: variation parameter
% ***********************************

for i=1:length(I)-1
    U1{i}=U0{i};
end
X0=[];
for t=1:I(length(I))
    for i=1:length(I)-1
        U1{i} = (1-v) * U1{i} + v * randn(I(i),R(i));
    end
    U1{length(I)} = U0{3}(t,:); 
    XR=full(ktensor(U1));
    sz=I(1:length(I)-1);
    
    idx = repmat({':'}, 1, length(I)-1);
    idx{length(I)}=t;
    X0(idx{:})=reshape(XR,sz);
    
end  


