function [lines,xy_long] = houghTest(skel,numpeaks)
% Input: Image skeleton and number of peaks.
% Output: The straightest lines in the skeleton and coordinates of the longest line.
% Using a Hough Transform, there is a first parametrisation of the skeleton. 
% Using that, we select the Hough Peaks using number of peaks. 
% We then parametrise back in order to obtain the straightest lines corresponding to the brightest peaks. 
% All of this is plotted on the display presenting the skeleton.


[H,T,R] = hough(skel,'RhoResolution',1); % default RhoResolution is 1
% old = H;

% figure
% % H =  (old-min(old(:))) ./ (max(old(:)-min(old(:))));
% % imshow(H/255,'XData',T,'YData',R,'InitialMagnification','fit');
% imagesc(H,'XData',T,'YData',R)
% title('Hough transform');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% colormap(gca,hot);
P = houghpeaks(H,9,'Threshold',0);% default Threshold is half numpeaks

% x = T(P(:,2));
% y = R(P(:,1));
% plot(x,y,'s','color','black');


% Length of lines gives us the number of merged line segments found
lines = houghlines(skel,T,R,P,'FillGap',5,'MinLength',7);

% figure, imshow(skel), hold on
hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

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
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
hold off
end