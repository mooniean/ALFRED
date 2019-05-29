function [totalSize,path,bigx,bigy] = lengthFinding(hObject,handles)
% Function to calculate the length of the image object. The user chooses
% a certain number of points (numPoints) and it tries to find the shortest
% path between the points (even if they weren't clicked consecutively).
% Furthermore, if there is a gap, it has a solution for such as it tries to
% find the shortest path to the closest connected component of the next
% point.

% Input: UI variables hObject and handles.
% Output: Total length, list of numbers constituting the new path, X coordinates of new path, Y coordinates of new path.
% This function finds and calculates the path between the two starting points (beginning and ending of the path) that the user has chosen.
% It heavily uses the function evaluateShortest, and returns the updated and concluded values. Most of the needed input arrives from the handles object.

totalSize=0;
if isempty(handles.nodenumber)
    bigx=handles.graph.Nodes(:,:).x;
    bigy=handles.graph.Nodes(:,:).y;
    path = [];
    totalSize = 0;
else
    nStops = length(handles.nodenumber);
    
    
    if nStops > 2
        beginningPoint=handles.nodenumber(1);
        path = [];
        donePoints = [1 2];
        while 1
            if ismember(beginningPoint,handles.linematrix)
                % Estou a fazer isto mal... pensa bem o que é o beginningPoint
                % e o que é o endNumber que é para andares com índices
                % correctos! o que tu queres é ver se o nodenumber
                % (beginningpoint) existe na matriz da linha. se existir, então
                % esse ponto está feito (adicionar o numero do ponto aos done
                % points) e saltas para o próximo nó na matriz, que vai ser
                % handles.linematrix(find(handles.linematrix==beginningPoint)+1)
                
                % endNumber = index of the point in the nodenumber matrix where
                % this path ends.
                endNumber = find(handles.nodenumber==handles.linematrix(find(handles.linematrix==beginningPoint)+1));
                % Add pathsize = distance between the two nodes
                pathsize = sqrt((handles.graph.Nodes(beginningPoint,:).x-handles.graph.Nodes(handles.nodenumber(endNumber),:).x)^2+(handles.graph.Nodes(beginningPoint,:).y-handles.graph.Nodes(handles.nodenumber(endNumber),:).y)^2);
                %             pathsize = sqrt((handles.g.Nodes(beginningPoint,:).x-handles.g.Nodes(handles.nodenumber(endNumber),:).x)^2+(handles.g.Nodes(beginningPoint,:).y-handles.g.Nodes(handles.nodenumber(endNumber),:).y)^2);
                % Add path = straight line between them, does this mean it's
                path = [path beginningPoint handles.nodenumber(endNumber)];
                
                % just the two points that I need to add, and when printing, we
                % print the line?
            else
                [shortest,endNumber,pathsize] = evaluateShortest(handles.graph,handles.nodenumber,beginningPoint,donePoints);
                path = [path shortest]; %#ok<*AGROW>
            end
            
            donePoints = [donePoints endNumber];
            beginningPoint = handles.nodenumber(endNumber);
            totalSize = totalSize + pathsize;
            if size(donePoints,2)>=size(handles.nodenumber,2)
                donePoints(2)=[];
                [shortest,~,pathsize] = evaluateShortest(handles.graph,handles.nodenumber,beginningPoint,donePoints);
                path = [path shortest];
                totalSize = totalSize + pathsize;
                break
            end
            
        end
    elseif nStops>0
        donePoints = 1;
        [path,~,totalSize] = evaluateShortest(handles.graph,handles.nodenumber,handles.nodenumber(1),donePoints);
    end
    bigx=handles.graph.Nodes(path(1,:),:).x;
    bigy=handles.graph.Nodes(path(1,:),:).y;
end
    
end