function i = zscore(i,DIM)
% ZSCORE removes the mean and normalizes the data 
% to a variance of 1. 
%
% z = zscore(x,DIM)
%   calculates the z-score of x along dimension DIM
%   it removes the 
%
% DIM	dimension
%	1: STATS of columns
%	2: STATS of rows
%	default or []: first DIMENSION, with more than 1 element
%
% features:
% - can deal with NaN's (missing values)
% - dimension argument 
% - compatible to Matlab and Octave
%
% see also: SUMSKIPNAN, MEAN, STD, DETREND
%
% REFERENCE(S):
% [1] http://mathworld.wolfram.com/z-Score.html

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


%	Copyright (C) 2000-2003 by  Alois Schloegl	<a.schloegl@ieee.org>	
%	$Revision: 1.14 $
%	$Id: zscore.m,v 1.14 2003/10/09 09:00:50 schloegl Exp $


if any(size(i)==0); return; end;

if nargin < 2,
	[S,N,SSQ] = sumskipnan(i);		% sum
else
	[S,N,SSQ] = sumskipnan(i,DIM);		% sum
end;

M = S./N;
i = i - repmat(M,size(i)./size(S));		% remove mean
i = i.*repmat(sqrt((SSQ-real(S).*real(M)-imag(S).*imag(M)).\max(N-1,0)),size(i)./size(S));	 % normalize by STD


