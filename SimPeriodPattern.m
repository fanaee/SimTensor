function [x,xn]=SimPeriodPattern(N,pat,prd,ns,gr)

% N=length of vector
% pat: value for cycle( 0.5 for summer and 0.25 for winter)
% prd: length of each cycle (e.g. for summer 93 days, for )
% ns: noise level
% gr: growth year in the second life cycle

if nargin<6
    cp=[];
end    

if nargin<5
    cp=[];
    gr=0;
end

if nargin<4
    cp=[];
    gr=0;
    ns=0;
end

cycle_total=sum(prd);
x1=[];
for i=1:length(pat)
    ps=repmat ( pat(i),1,prd(i));
    x1=[x1 ps];
end

x=x1;

if N>cycle_total
    for i=1:ceil(N/cycle_total)
         x=[x (1+i*gr)*x1];
    end
end    
x=x(1:N);

x=x+ns*randn(size(x));


xn=zscore(x);
