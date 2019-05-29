function partialLength = lengthCalculation(X,Y)
% Input: Points X and Y.
% Output: Partial length
% This function calculates the Euclidean distances between two points X and Y using their matrix coordinates (i,j), i.e., checks whether the next neighbour is on the same
% row or column, or if they are diagonally connected. 
% The values of distance are as follows:
% - 1 - same row or same column;
% - 2 - diagonally connected, therefore the distance is sqrt(2)

% Basically, if the totalDiff is 1, it means that the points are on the
% same row or the same column (distance of 1) but if totalDiff is 2, it
% means they are diagonally connected, which has the distance of sqrt(2)
% (approx 1.41)

    partialLength = 0;
    xDiff=diff(X);
    yDiff=diff(Y);
    totalDiff = abs(xDiff) + abs(yDiff);
    for i = 1:length(xDiff)
        if totalDiff(i) == 1
            partialLength = partialLength + 1;
        else
            partialLength = partialLength + sqrt(2);
        end
    end
end