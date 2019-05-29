function distance = checkDistance(x1,x2,y1,y2)
% Input: Coordinates x1,x2, y1 and y2.
% Output: Distance.
% Measures the Euclidean distance between point (x1,y1) and point (x2,y2).


distance = sqrt((x1-x2)^2+(y1-y2)^2);
end