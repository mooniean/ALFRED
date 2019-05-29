% clear; close all; clc;
% INTERACTIVE=true;

function [dTdt, dsdt, xt, yt, dsdT] = curvatureGauss(coordX,coordY,sigma)


%For the test images, uncomment the intended from the following lines:
% %obtain image skeleton
% image=sum(imread('Slide1.tiff'),3);
% skeleton=bwmorph(imrotate(image>250,0),'skel',Inf);
% TestImages=imread('TestImagesWhiteOnBlack.tiff');
% skeleton=bwmorph(imrotate(TestImages(50:1210,1250:2440)>250,0),'skel',Inf);

%Otherwise, the skeleton is saved in the handles object of the GUI

%build a binary image out of the coordinates
offset=20; %zero-pad image so I have like a frame around it
skeleton = zeros(max(coordX)+offset,max(coordY)+offset);
for i = 1:length(coordX)
    skeleton(coordX(i)+offset/2,coordY(i)+offset/2) = 1;
end

%find branchpoints
brPoints=bwmorph(skeleton,'branchpoints');

%eliminate branchpoints and their 24 nearest neighbours from the skeleton
[rowBr,colBr] = find(brPoints);
for k=1:length(rowBr)
    for i=-2:2
        for j=-2:2
            skeleton(rowBr(k)+i,colBr(k)+j)=0;
        end
    end
end

%clean groups of up to 10 pixels that are isolated in the skeleton
skeleton=bwareaopen(skeleton,10);

%find boundaries for objects not holes (better performance)
% boundaries=bwboundaries(skeleton,'noholes');
boundaries=bwboundaries(skeleton);

%storage for x(t) and y(t), where the curvature computation is valid
xt=cell(length(boundaries),1);
yt=cell(length(boundaries),1);

%storage for x'(t) and y'(t), where the curvature computation is valid
dxdt=cell(length(boundaries),1);
dydt=cell(length(boundaries),1);

%storage for derivatives of T (the curve's tangent vector) and s (the curve
%'s arc length) for each boundary
dTdt=cell(length(boundaries),1);
dsdt=cell(length(boundaries),1);

for b=1:length(boundaries)
    
    x=boundaries{b}(:,1).';y=boundaries{b}(:,2).';
    
    %for line-style objects, the boundaries will include a back and forth
    %trace, which produces a mirrored path: here we eliminate one of them
    start=1;
    for i=3:length(x)
        if (y(i)==y(i-2) && x(i)==x(i-2))
            start=i;break;
        end
    end
    if (start~=1)
        x=x(start-1:start-1+(length(x)+1)/2-1);
        y=y(start-1:start-1+(length(y)+1)/2-1);
    end
    
    %t is the parameter in relation to which the curve is parametrized xD
    %     t=1:length(x);
    
    %here we use a gaussian filter to compute the first and second
    %derivatives, x'(t) x''(t) and y'(t) y''(t), of the parametrizations, x(t) and y(t)
    
    %the std dev of the gaussian
    if nargin<3; sigma=25; end
    
    %the length of the filter: radius=2*sigma
    fLen=4*sigma+1;
    
    %the gaussian filter
    %Gf=fspecial('gauss',[1 fLen], sigma);
    Gf = @(x) 1/sigma/sqrt(2*pi)*exp(-1/2*(x/sigma).^2);
    G = Gf(-(fLen-1)/2:(fLen-1)/2);
    factor=sum(G);
    G=G/factor;
    
    %smooth x(t) and y(t)
    xt{b}=conv(x,G,'valid');
    yt{b}=conv(y,G,'valid');
    
    %compute the first derivative of the gaussian filter and scale it so
    %that the result of its application is rightly scaled
    %G1=gradient(Gf);
    Gf1 = @(x) -x/sigma^2.*Gf(x);
    G1 = Gf1(-(fLen-1)/2:(fLen-1)/2);
    factor=sum(-G1.*((1:fLen)-(fLen+1)/2));
    G1=G1/factor;
    
    %use the above filter to estimate x'(t) and y'(t)
    dxdt{b}=conv(x,G1,'valid');
    dydt{b}=conv(y,G1,'valid');
    
    %compute the second derivative of the gaussian filter, force the 0th
    %term to annihilate the others' contributions to f(x) in the Taylor
    %Series and scale it so that the result of its application is rightly
    %scaled
    %G2 = gradient(gradient(Gf));
    Gf2 = @(x) (x.^2-sigma^2)/sigma^4.*Gf(x);
    G2 = Gf2(-(fLen-1)/2:(fLen-1)/2);
    G2((fLen+1)/2)=-sum(G2)+G2((fLen+1)/2);
    factor=sum(G2.*(((1:fLen)-(fLen+1)/2).^2))/2;
    G2=G2/factor;
    
    %use the above filter to estimate x''(t) and y''(t)
    d2xdt2=conv(x,G2,'valid');
    d2ydt2=conv(y,G2,'valid');
    
    %computing dT/dt=dT/ds*ds/dt=k*ds/dt and ds/dt
    dTdt{b}=(dydt{b}.*d2xdt2-d2ydt2.*dxdt{b})./(dxdt{b}.^2+dydt{b}.^2);
    dsdt{b}=sqrt(dxdt{b}.^2+dydt{b}.^2);    
end

%from here on, we have all the data, this is just a fancy way of showing
%the data we have... it's better if we make this interactive
%and to make it interactive you HAVE to output the positions as well
dsdT = [];
for b=1:length(boundaries)
    %compute all the available radii for the current boundary
    %Matthias hates this, so we shouldn't use this anymore! (:
    radii=dsdt{b}./abs(dTdt{b});
    dsdT = [dsdT radii]; %#ok<*AGROW>
end

% % if (~INTERACTIVE)
%     for b=1:length(boundaries)
%
%         %compute all the available radii for the current boundary
%         radii=dsdt{b}./dTdt{b};
%         results = [results radii]; %#ok<*AGROW>
%
%         %draw the current boundary region under inspection
%         region=zeros(size(skeleton));
%         for j=1:length(radii); region(xt{b}(j),yt{b}(j))=1; end
%
%         %compute pointwise curvature and osculating circle
%         for i=1:3:length(radii)
%
%             %centre of the osculating circle
%             centre=round([xt{b}(i) yt{b}(i)]...
%                 +radii(i)*[dydt{b}(i) -dxdt{b}(i)]/dsdt{b}(i));
%
%             %let's draw the osculating circle: creating the canvas
%             oscCircle=zeros(size(skeleton));
%
%             %for the following points, draw on (2l+1)^2 pixels
%             l=2;
%
%             %draw the current point under inspection
%             oscCircle(xt{b}(i)-l:xt{b}(i)+l,yt{b}(i)-l:yt{b}(i)+l)=1;
%
%             %draw the the circle's center
%             if (centre(1)-l>=1 && centre(1)+l<=size(skeleton,1) &&...
%                     centre(2)-l>=1 && centre(2)+l<=size(skeleton,2))
%                 oscCircle(centre(1)-l:centre(1)+l,centre(2)-l:centre(2)+l)=1;
%             end
%
%             %draw the circle itself
%             for k=0:.01:2*pi
%                 if (centre(1)+round(radii(i)*cos(k))>=1 &&...
%                         centre(1)+round(radii(i)*cos(k))<=size(skeleton,1)&&...
%                         centre(2)+round(radii(i)*sin(k))>=1 &&...
%                         centre(2)+round(radii(i)*sin(k))<=size(skeleton,2))
%                     oscCircle(centre(1)+round(radii(i)*cos(k)),...
%                         centre(2)+round(radii(i)*sin(k)))=1;
%                 end
%             end
%
%             %join and show the current region, the skeleton and the osculating
%             %circle
%             imshow(cat(3,region,0.5*skeleton,oscCircle));
%
%             %write the image sequence to a GIF
%             %capture the plot as an image
%             im = frame2im(getframe(gcf));
%             [imind,cm] = rgb2ind(im,256);
%             %write to the GIF file
%             if i == 1
%                 imwrite(imind,cm,'curvature.gif','DelayTime',0,...
%                     'WriteMode','overwrite');
%             else
%                 imwrite(imind,cm,'curvature.gif','DelayTime',0,...
%                     'WriteMode','append');
%             end
%         end
%     end
% end


% %example of an interactive implementation
% if (INTERACTIVE)
%     %show image and make it mouse interactable
%     imshow(skeleton); hold on;
%     dcm_obj = datacursormode(gcf);
%     set(dcm_obj,'Enable','on','UpdateFcn',@displayRadius)
% end


% guidata(hObject,handles)
end


% function txt=displayRadius(~,event_obj,handles)
%
% %get cursor information
% curPos = get(event_obj,'Position');
%
% %get the coordinates of the point under inspection
% y=curPos(1);
% x=curPos(2);
%
% %detect the boundary that contains said point
% %range of search: (2l+1)^2 pixels
% l=2;
% foundPoint=false;
% for b=1:length(boundaries)
%     for j=1:length(xt{b})
%         if (xt{b}(j)>=x-l && xt{b}(j)<=x+l &&...
%                 yt{b}(j)>=y-l && yt{b}(j)<=y+l)
%             foundPoint=true;break;
%         end
%     end
%     if foundPoint; break; end
% end
%
% %if point falls outside a boundary or lacks a curvature computation
% %the latter might be caused by a large sigma
% if ~foundPoint
%     txt={'Point not available','Try reducing sigma'};
%     imshow(skeleton); return;
% end
%
% %compute the radius
% radius=handles.dsdt{b}(j)/handles.dTdt{b}(j);
% txt=num2str(abs(radius));
%
% %compute pointwise curvature and osculating circle
% %centre of the osculating circle
% centre=round([xt{b}(j) yt{b}(j)]...
%     +radius*[dydt{b}(j) -dxdt{b}(j)]/handles.dsdt{b}(j));
%
% %let's draw the osculating circle: creating the canvas
% oscCircle=zeros(size(skeleton));
%
% %for the following point, draw on (2l+1)^2 pixels
% l=2;
%
% %draw the the circle's center
% if (centre(1)-l>=1 && centre(1)+l<=size(skeleton,1) &&...
%         centre(2)-l>=1 && centre(2)+l<=size(skeleton,2))
%     oscCircle(centre(1)-l:centre(1)+l,centre(2)-l:centre(2)+l)=1;
% end
%
% %draw the circle itself
% for k=0:.01:2*pi
%     if (centre(1)+round(radius*cos(k))>=1 &&...
%             centre(1)+round(radius*cos(k))<=size(skeleton,1)&&...
%             centre(2)+round(radius*sin(k))>=1 &&...
%             centre(2)+round(radius*sin(k))<=size(skeleton,2))
%         oscCircle(centre(1)+round(radius*cos(k)),...
%             centre(2)+round(radius*sin(k)))=1;
%     end
% end
%
% %draw the current boundary region under inspection
% region=zeros(size(skeleton));
% for i=1:length(xt{b}); region(xt{b}(i),yt{b}(i))=1; end
%
% %join and show the current region, the skeleton and the osculating
% %circle
% imshow(cat(3,region,0.5*skeleton,oscCircle));
%
% end