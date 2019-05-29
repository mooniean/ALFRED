function requireR2015b

% Copyright 2015 The MathWorks, Inc.
% Steve Eddins
% Edited by Nuno Nobre (02/01/2019)

if verLessThan('matlab','8.6')
   throwAsCaller(MException('se:ig:RequiresR2015b', ...
      'MATLAB R2015b is required.'));
end
