function [x,y]=gacu(Int,Len,Mu,Sigma)
%Generates a gaussian curve of mean Mu and standard deviation Sigma
%The x values (Len values between Int(1) and Int(2)) are given in x.
%y contains the values of the curve.
x = linspace(Int(1),Int(2),Len);
y = sqrt(2 * pi * Sigma)^-1 * exp(-(x - Mu).^2 / (2 * Sigma^2));
