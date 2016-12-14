function [U,Un]=SimFacGen(method,I,R,ns,normz)

% alternatively use 
% U0 = arrayfun(@(x) rand(x,R),I,'uni',0); 

  
if nargin<5
    normz=0;
end    
if nargin<4
    normz=0;
    ns=0;
end    

U = [];
%{
gamma
multi_normal_dist
rand
randn
orthogonal
stochastic
binary
multiple_normal_mu_sigma
%}
switch method
        case 'gamma' % code of Piyush Rai (ECML-PKDD 2015 paper)
            U=gamrnd(1e-1+1e-1*abs(randn(I,R)),1e-2);
        case 'multi_normal_dist' % code of Qibin Zhao (TPAMI 2014 paper)
            U =  gaussSample(zeros(R,1), eye(R), I);
        case 'rand' % uniform on [0,1]
            U = rand(I,R);
        case 'randn' % standard normal distribution
            U = randn(I,R);
        case 'orthogonal' % from create_problem in Tensortoolbox 2.6   
            tmp = matrandorth(max(I,R));
            U = tmp(1:I,1:R);
        case 'stochastic' % uniform on [0,1] with column sums rescaled to 1, from create_problem Tensortoolbox 2.6
            tmp = rand(I, R);
            S = sum(tmp,1);
            U = tmp * diag(1./S);
        case 'binary' 
            % code of Yasuko Matsubara (KDD 2012 paper)    
            for k=1:I
                a=zeros(1,R);
                a(1,floor(rand*R)+1)=1;
                a = a./sum(a);
                U=[U; a];
            end    
        case 'multiple_normal_mu_sigma' % code of Giorgio Tomasi (Chemometrics and Intelligent Laboratory Systems 2005 paper)
            T = randperm(max(I,R));
            for k=1:R
                [~,U(k,:)] = gacu([1 I],I,T(k),I*0.15*(rand(1)+0.1));
            end
            U=U';
        case 'mvnrnd' %under construction
            data = [normrnd(0,1,5000,1),normrnd(0,1,5000,1)]; 
            MU = mean(data,1);
            SIGMA = cov(data);
            r = mvnrnd(MU,SIGMA,5000);

end

if ns>0
    U=U+ns*randn(size(U));
end

   
if normz>0
    switch normz
        case 2
            Un=U./(ones(size(U,1),1)*sqrt(sum(U.^2)));
        case 3
            Un = bsxfun(@rdivide,U,sqrt(sum(U.^2))); % code of Phan Anh-Huy (tensorbox)
        case 4
            Fict = sqrt(diag(sum((U').^2)));Un = U'*Fict^-1; Un=Un'; % code of Giorgio Tomasi
        otherwise
            Un=U./repmat(sum(U),I,1); % code of Yasuko Matsubara & Piyush Rai & Qibin Zhao
    end
end     