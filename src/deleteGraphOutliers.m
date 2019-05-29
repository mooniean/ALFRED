function [newNodes,setofcoordinates,temp] = deleteGraphOutliers(graph,node1,node2)
% Input: Image graph, and node number for node1 and node2.
% Output: Updated list of node numbers, set of coordinates and temporary graph.
% 
% This is a function used whenever there is a need to select just the area between two user input points. 
% It works as follows:
% 1. See the node neighbours of both nodes in the graph.
% 2. For each neighbour of node1, see if there is a connected path between them and node2. 
%    If that path includes node1, we remove add the node to the removal list.
% 3. Repeat for the neighbours of node2 with regards to node1.
% 4. Remove nodes in the removal list from the graph. 
%    Note: this will update the node numbers. As such, there is a need to return from this function the updated values of node numbers, coordinates and the new graph.



%     temp = graph;
    nodeOneNeighbors = neighbors(graph,node1);
    nodeTwoNeighbors = neighbors(graph,node2);
    removeNodes = [];
    for i = 1:size(nodeOneNeighbors,1)
        if ~isempty(find(shortestpath(graph,nodeOneNeighbors(i),node2)==node1)) % não deveria estar a checkar o temp, em vez do graph?
            removeNodes = [removeNodes nodeOneNeighbors(i)]; %#ok<*AGROW>
        end
    end

    for i = 1:size(nodeTwoNeighbors,1)
        if ~isempty(find(shortestpath(graph,nodeTwoNeighbors(i),node1)==node2)) %#ok<*EFIND>
            removeNodes = [removeNodes nodeTwoNeighbors(i)];
        end
    end
    temp = rmnode(graph,removeNodes);
    newNode1x = graph.Nodes(node1,1).x;
    newNode1y = graph.Nodes(node1,2).y;
    newNode1 = find(table2array(temp.Nodes(:, 1)) == newNode1x & table2array(temp.Nodes(:, 2)) == newNode1y);
    connectednodes   = conncomp(temp); %aposto que podes determinar só a connected component do newnode1 directamente
    bin              = connectednodes(newNode1); %#ok<FNDSB>
    newNodes         = find(connectednodes==bin);
    temp = subgraph(temp, newNodes);
    setofcoordinates = temp.Nodes(:,:);
end