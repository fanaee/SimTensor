function lambda=SimLamGen(i_struct,R,method,topn,w)


if strcmp(i_struct,'Tucker') %Tucker
    
    switch method
        case  'rand'
            rndt=rand(R);
        case 'randn'
            rndt=randn(R);
        otherwise
            rndt=ones(R);
    end
    
    for i=1:length(R)
        idx{i}=1:round(topn*R(i)/100);
    end
    rndt(idx{:})=rndt(idx{:})*w;
    lambda= reshape(rndt,prod(R),1);


else % CP
    R=R(1);
    switch method
        case 'gamma'
            % from the code of Piyush Rai (ECML-PKDD 2015 paper)
            xr=10*randperm(500,R);
            a=abs(2+randn(1,R));
            lambda=gamrnd(xr,a)';
        case  'rand'
            lambda= rand(R,1);
        case 'randn'
            lambda= randn(R,1);
        otherwise
            lambda=ones(R,1);
    end   
    pp=round(topn*R/100);
    lambda(1:pp)=lambda(1:pp)*w;

end    

lambda=lambda';