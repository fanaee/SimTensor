function [X,O,Omega]=remove_random_elements(X0,ObsRatio,typ)
if nargin<3
	typ=0;
end	
if typ==0
	% ObsRatio: Observation Ratio
	I=size(X0);
	Omega = randperm(prod(I)); 
	Omega = Omega(1:round(ObsRatio*prod(I)));
	O = zeros(I); 
	O(Omega) = 1;
	X= O.*X0;
else
	I=size(X0);
    X0=double(X0);
	Omega = randperm(prod(I)); 
	Omega = Omega(1:round((1-ObsRatio)*prod(I)));
	X0(Omega)=NaN;
	O = ones(I); 
	O(Omega) = 0;
	Omega=find(O==1);
	X=X0;
end	