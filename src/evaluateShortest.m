function [shortest,endNumber,minSize] = evaluateShortest(g,nodenumber,beginningPoint,donePoints)
% Input: Image graph, list of all node numbers and the beginning node, as well as a list of the points that were already done.
% Output: Shortest path, number of the ending node, path length.
% 
% This function will calculate the shortest path between beginning node and all theavailable nodes in the node numbers list but not on the done points. 
% It also tries to circumvent points that are not necessarily connected, but rather need to have a straight line to the closest forward node in the graph, and then the distance between that new
% node and all the other node numbers. This will be calculated using an Euclidean distance, and we often call it "bridge gap".



minSize = Inf;
endNumber = 1;
shortest=[];
for i = 1:size(nodenumber,2)
    if ismember(i,donePoints)
        continue
    end
    lengthPathTemp = 0;
    path = shortestpath(g,beginningPoint,nodenumber(i));
    if isempty(path)
        % Code to find the closest point in the connected component of
        % the next point. BASICALLY, if there can't be a path, try to
        % find the closest possible point to have a path
        connectednodes = conncomp(g);
        bin = connectednodes(nodenumber(i));
        newNodes = find(connectednodes==bin);
        setofcoordinates=g.Nodes(newNodes,:); %#ok<*FNDSB>
        beginningPointX = g.Nodes(beginningPoint,:).x;
        beginningPointY = g.Nodes(beginningPoint,:).y;
        %radius = round(sqrt((beginningPointX-g.Nodes(nodenumber(i),:).x)^2+(beginningPointY-g.Nodes(nodenumber(i),:).y)^2));
        %[x,y,~] = checkCoordinates(beginningPointX,beginningPointY,setofcoordinates,radius) ;
        [~,index]=min(sum((table2array(setofcoordinates(:,1:2))-[beginningPointX beginningPointY]).^2,2));
        x=table2array(setofcoordinates(index, 1)); y=table2array(setofcoordinates(index, 2));
        [nodenumbertemp,~] = find(table2array(g.Nodes(:, 1)) == x & table2array(g.Nodes(:, 2)) == y);
        path = shortestpath(g,nodenumbertemp,nodenumber(i));
        lengthPathTemp = sqrt((beginningPointX-x)^2+(beginningPointY-y)^2);
    end
    bigx=g.Nodes(path(1,:),:).x;
    bigy=g.Nodes(path(1,:),:).y;
    lengthPathTemp = lengthPathTemp + lengthCalculation(bigx,bigy);
    if lengthPathTemp<minSize
        shortest = path;
        minSize = lengthPathTemp;
        endNumber = i;
    end
end
end