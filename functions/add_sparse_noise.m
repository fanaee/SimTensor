function [X1]=add_sparse_noise(X0,SparseRatio,typ)
I=size(X0);

SparseOmega = randsample(prod(I), round(prod(I)*SparseRatio));
S = zeros(I);

if typ==1
	S(SparseOmega) = max(X0(:)) * (2*rand(length(SparseOmega),1)-1);
	X1 = X0 + S;
else
	S(SparseOmega) = 10*std(X0(:))*(2*rand(length(SparseOmega),1)-1);
	X1 = X0 + S;
end	