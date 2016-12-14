function chkinputdatatype(varargin)
%CHKINPUTDATATYPE Check for unsupported single precision inputs

%   Copyright 2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/08/11 15:48:34 $


for n = 1:nargin
    if isa(varargin{n},'single')
        error(generatemsgid('NotSupported'),'Input arguments cannot be ''single''.');
    end
end



% [EOF]
