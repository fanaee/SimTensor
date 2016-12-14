function [x,xn]=SimPeriodWave(typ,N,p,frq,ns)

% typ: cos,sin,sawtooth,square
% N:length of signal
% p: number of periods
% frq: frequency in Hz
% ns: noise level

T = p*(1/frq);
dt = T/N;
t = 0:dt:T-dt;
switch typ
    case 'cos'
        x=cos(2*pi*frq*t)+ns*randn(size(t));
    case 'sin'
        x=sin(2*pi*frq*t)+ns*randn(size(t));
    case 'saw'
        x = sawtooth(2*pi*frq*t)+ns*randn(size(t));
    case 'sq'
        x = square(2*pi*frq*t)+ns*randn(size(t));
end

t= linspace(1,mean(x),N);
yinitial = exp(-t);
x=x+yinitial;


xn=zscore(x);