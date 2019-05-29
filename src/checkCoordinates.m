function [x,y,nodenumber] = checkCoordinates(oldX,oldY,setofcoordinates,radius)
% Input: Position X and Y of the selected user point, the set of coordinates that compose the graph, radius of search.
% Output: The closest coordinates of the graph to selected point, the node number to which those coordinates belong to.
% The function will try to find, within a reasonable radius, the closest (Euclidean distance) coordinates from the set of coordinates to the X and Y of the graphical selected user point on the UI display.


oldX = floor(oldX);
oldY = floor(oldY);
startX = max(oldX - radius,0);
endX   = oldX + radius;
startY = max(oldY - radius,0);
endY   = oldY + radius;
minDist = Inf;
x=0;y=0;
for i = startX:endX
    % Checking whether the X exists, cuts the time of a Y cycle
    if size(find(table2array(setofcoordinates(:, 1)) == i),1) == 0
        continue
    end
    for j = startY:endY
        if size(find(table2array(setofcoordinates(:, 1)) == i & table2array(setofcoordinates(:, 2)) == j ),1) == 0
            continue
        end
        tempDist = sqrt((i-oldX)^2+(j-oldY)^2);
        if tempDist < minDist
            minDist = tempDist;
            x=i;
            y=j;
        end
    end
end
[nodenumber,~] = find(table2array(setofcoordinates(:, 1)) == x & table2array(setofcoordinates(:, 2)) == y);
end