k=1;
distances=[];
while k<=sizeStraightLines
    try
        if (straightLines(k).rho == straightLines(k+1).rho) && (straightLines(k).theta == straightLines(k+1).theta)
            temp = tempROIs(9).straightlineDistances(k);
            values = [k];
            while (k+1<=sizeStraightLines && straightLines(k).rho == straightLines(k+1).rho) && (straightLines(k).theta == straightLines(k+1).theta)
                %                         distances = [distances finalROIs(i).straightlineDistances(k)+finalROIs(i).straightlineDistances(k+1)];
                %                             temp = temp + handles.finalROIs(i).straightlineDistances(k+1);
                
                k=k+1;
                values = [values k];
            end
            coordinates = [];
            for tempindex = values
                coordinates = [coordinates; straightLines(tempindex).point1;straightLines(tempindex).point2];
            end
            temp = checkDistance(coordinates(1,1),coordinates(end,1),coordinates(1,2),coordinates(end,2));
            distances = [distances temp];
        else
            distances = [distances tempROIs(9).straightlineDistances(k)];
        end
    catch
        distances = [distances tempROIs(9).straightlineDistances(k)];
    end
    k=k+1;
end