function image = concatCells(cellarray)
% Concatenates cells from a cell array.
L=size(cellarray,1);
image = [];
for i = 1:L
    image = cat(3,image,cellarray{i,1});
end

end