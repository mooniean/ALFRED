function [z,index] = findValue(structurePassed, valueToSearch)
% Finds a certain value in a structure.


% Squeeze the struct into a cell
C = squeeze(struct2cell(structurePassed))';

% Find the value in the cell cell. @(m) means an anonymous function, which
% will be whatever we want.
% found = cellfun(@(m) any(strcmp(m, valueToSearch)), C);
[index,~] = find( cellfun(@(x)isequal(x,valueToSearch),C) );

% Count how many, return that value
z = length(index);
end