function [donelines,xy_long, donepath,donepoints,distances] = houghTest(skel,numpeaks,recurse)
% Input: Image skeleton and number of peaks.
% Output: The straightest lines in the skeleton and coordinates of the longest line.
% Using a Hough Transform, there is a first parametrisation of the skeleton.
% Using that, we select the Hough Peaks using number of peaks.
% We then parametrise back in order to obtain the straightest lines corresponding to the brightest peaks.
% All of this is plotted on the display presenting the skeleton.


[H,T,R] = hough(skel,'RhoResolution',3); % default RhoResolution is 1



% old = H;
% figure
% % H =  (old-min(old(:))) ./ (max(old(:)-min(old(:))));
% % imshow(H/255,'XData',T,'YData',R,'InitialMagnification','fit');
% imagesc(H,'XData',T,'YData',R)
% title('Hough transform');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% colormap(gca,hot);
% size(H)


P = houghpeaks(H,numpeaks,'Threshold',ceil(0.3*max(H(:))));% default Threshold is half numpeaks, 15,'NHoodSize'

% x = T(P(:,2));
% y = R(P(:,1));
% plot(x,y,'s','color','white');


% Length of lines gives us the number of merged line segments found
lines = houghlines(skel,T,R,P,'FillGap',5,'MinLength',50);


%NOT THE EUCLIDEAN BETWEEN POINTS BUT THE ACTUAL PATH FROM SKEL - PERHAPS
%CHECK THE GRAPHS AND USE THEM?
graph = binaryImageGraph(skel);
% HAVE A TEMPORARY MATRIX THAT IS THE SAME AS THE PATH
% WHENEVER SOME LINES ARE DONE, BECOME ZERO
% iF POINTS NOT IN TEMPORARY MATRIX, SKIP

% if ~recurse
% figure, imshow(skel)
% end
hold on
max_len = 0;
donepath = [];
donepoints =[];
donelines = [];
xy_long = [];
distances = [];
for k = 1:length(lines)
    
    oldxy = [lines(k).point1; lines(k).point2];
%     tempLine = line(xy(:,1),xy(:,2));
    
    [x1,y1,nodenumber1] = checkCoordinates(oldxy(1,1),oldxy(1,2),graph.Nodes,20);
    [x2,y2,nodenumber2] = checkCoordinates(oldxy(2,1),oldxy(2,2),graph.Nodes,20);
    xy = [x1 y1; x2 y2];
    nodenumber = [nodenumber1 nodenumber2];
   
    x = xy(:,1);
    y = xy(:,2);

    p = [lines(k).rho lines(k).theta];
    flag = true;
%     donepoints
    [shortest,endNumber,sizePath] = evaluateShortest(graph,[nodenumber1 nodenumber2],nodenumber1,[1]);
    if sum([donepath==nodenumber2 donepath==nodenumber1])>0
        if (length(intersect(shortest,donepath))/length(shortest))>0.1
            continue
        end
    end
    distances = [distances; sizePath];
    donepoints = [donepoints endNumber];
    donepath = unique([donepath shortest]);
    plotx=graph.Nodes(shortest(1,:),:).x;
    ploty=graph.Nodes(shortest(1,:),:).y;
    if flag
%         doneLines = [doneLines; p k];
        plot(plotx,ploty,'LineWidth',2,'Color','green');
%         p
        donelines = [donelines lines(k)];
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red'); 
        
        % Determine the endpoints of the longest line segment
        len = norm(lines(k).point1 - lines(k).point2);
        
        if ( len > max_len)
            max_len = len;
            xy_long = xy;
        end
        
    end
end
% 
if length(donepath)/size(graph.Nodes,1)<1 && ~isempty(donepath)
    if length(donepoints)<numpeaks
        tempGraph = rmnode(graph,donepath);
        tempPath=zeros(size(skel));
        for i = 1:length(tempGraph.Nodes.x)
            tempPath(tempGraph.Nodes.y(i),tempGraph.Nodes.x(i))=1;
        end
        [templines,~,tempdonepath,tempdonepoints,tempdistances] = houghTest(tempPath,numpeaks,true);
        donelines = [donelines templines];
        donepath = unique([donepath tempdonepath]);
        donepoints = [donepoints tempdonepoints];
        distances = [distances; tempdistances];
    end
end

% highlight the longest line segment
% plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
% hold off
end