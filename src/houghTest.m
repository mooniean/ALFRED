function [lines,xy_long] = houghTest(skel,numpeaks)
% Input: Image skeleton and number of peaks.
% Output: The straightest lines in the skeleton and coordinates of the longest line.
% Using a Hough Transform, there is a first parametrisation of the skeleton.
% Using that, we select the Hough Peaks using number of peaks.
% We then parametrise back in order to obtain the straightest lines corresponding to the brightest peaks.
% All of this is plotted on the display presenting the skeleton.


[H,T,R] = hough(skel,'RhoResolution',8); % default RhoResolution is 1
old = H;

figure
% H =  (old-min(old(:))) ./ (max(old(:)-min(old(:))));
% imshow(H/255,'XData',T,'YData',R,'InitialMagnification','fit');
imagesc(H,'XData',T,'YData',R)
title('Hough transform');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(gca,hot);
size(H)
P = houghpeaks(H,numpeaks,'Threshold',ceil(0.3*max(H(:))));% default Threshold is half numpeaks, 15,'NHoodSize'

x = T(P(:,2));
y = R(P(:,1));
plot(x,y,'s','color','white');pause


% Length of lines gives us the number of merged line segments found
lines = houghlines(skel,T,R,P,'FillGap',5,'MinLength',50);


%NOT THE EUCLIDEAN BETWEEN POINTS BUT THE ACTUAL PATH FROM SKEL - PERHAPS
%CHECK THE GRAPHS AND USE THEM?
% HAVE A TEMPORARY MATRIX THAT IS THE SAME AS THE PATH
% WHENEVER SOME LINES ARE DONE, BECOME ZERO
% iF POINTS NOT IN TEMPORARY MATRIX, SKIP

figure, imshow(skel), hold on
hold on
max_len = 0;
doneLines = [];
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
%     tempLine = line(xy(:,1),xy(:,2));
    flag = true;
    
    
    
    x = [lines(k).point1(1) lines(k).point2(1)];
    y = [lines(k).point1(2) lines(k).point2(2)];
%     p=polyfit(x,y,1
    p = [lines(k).rho lines(k).theta];
    flag = true;
    for index = 1:size(doneLines,1)
%         if (p(1)==doneLines(index,1) && p(2)==doneLines(index,2))
            lineIndex = doneLines(index,3);
            %Checking if X is in the list of X between the points and if the X
            %calculated with the y=mx+b from up there is actually well
            %calculated
            x_int = linspace(lines(lineIndex).point1(1),lines(lineIndex).point2(1),abs(lines(lineIndex).point1(1)-lines(lineIndex).point2(1))+1);
            y_int = linspace(lines(lineIndex).point1(2),lines(lineIndex).point2(2),abs(lines(lineIndex).point1(2)-lines(lineIndex).point2(2))+1);
            if (ismember(x(1),x_int) || ismember(x(2),x_int))
                if (ismember(y(1),y_int) || ismember(y(2),y_int))
                    %             if sum(ismember(round(polyval([doneLines(index,1) doneLines(index,2)],x)),linspace(lines(lineIndex).point1(2),lines(lineIndex).point2(2),abs(lines(lineIndex).point1(2)-lines(lineIndex).point2(2))+1)))>0
%                     flag = false;
                    break
                end
            end
%         end
    end
    if flag
        doneLines = [doneLines; p k];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        p
        
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');pause
        
        % Determine the endpoints of the longest line segment
        len = norm(lines(k).point1 - lines(k).point2);
        
        if ( len > max_len)
            max_len = len;
            xy_long = xy;
        end
        
    end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
hold off
end