function [curvatures,media,numPoints,straightness] = compCurv(binaryImage,objIndex,windowSize,junctions,~)

warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

polyDegree=2; %está a fittar a um polinómio de grau 2 (parábola)
boundaries = bwboundaries(binaryImage);

x = boundaries{objIndex}(:,2);
y = boundaries{objIndex}(:,1);

assignin('base','x',x)
assignin('base','y',y)
halfWidth = floor(windowSize/2);   %half window size
curvatures = zeros(size(x));

for k = halfWidth+1 : length(x) - halfWidth
    theseX = x(k-halfWidth:k+halfWidth);
    theseY = y(k-halfWidth:k+halfWidth);
    
    try
        if sum(sum([junctions(theseY,theseX),junctions(theseY+1,theseX),junctions(theseY-1,theseX),junctions(theseY,theseX+1),...
                junctions(theseY,theseX-1),junctions(theseY+1,theseX+1),junctions(theseY-1,theseX-1),junctions(theseY+1,theseX-1),junctions(theseY-1,theseX+1)]))>0
            continue
        end
    catch
        continue
    end
    
    if ~issorted(theseX) && ~issorted(flip(theseX)); curvatures(k)=1000; continue; end; % to be filtered ahead. reason: x not monotonic
    [uniqueX,i]=unique(theseX);uniqueY=theseY(i);
    
    if length(uniqueY)<0.5*length(theseY);
        [uniqueX,i]=unique(theseY);uniqueY=theseX(i);
    end
    
    if length(uniqueX)<=polyDegree; curvatures(k)=1000; continue; end % to be filtered ahead. reason: not enough points
    
    if polyDegree==2 %this is just scales, does the fit, and returns to non-scaled form... It avoids the warning (I'm now supressing it anyway), but the results seem to be the same
        [coefficientsScaled,~,mu]=polyfit(uniqueX, uniqueY, polyDegree);
        coefficientsScaled=flip(coefficientsScaled);
        coefficients=[coefficientsScaled(3)*(mu(1)/mu(2))^2-coefficientsScaled(2)*mu(1)/mu(2)+coefficientsScaled(1),...
            -2*mu(1)*coefficientsScaled(3)/mu(2)^2+coefficientsScaled(2)/mu(2),...
            coefficientsScaled(3)/mu(2)^2];
    else
        coefficients= flip(polyfit(uniqueX, uniqueY, polyDegree));
    end
    allTerms=0:polyDegree;
    firstDer=sum(coefficients.*allTerms.*(x(k).^(allTerms-1)));
    secondDer=sum(coefficients.*allTerms.*(allTerms-1).*(x(k).^(allTerms-2)));
    curvatures(k) = secondDer/(1+firstDer^2)^(3/2);

    
    dreal=0;
    for i = 2:length(uniqueX)
        dreal = dreal + sqrt((uniqueX(i)-uniqueX(i-1))^2+(uniqueY(i)-uniqueY(i-1))^2);
    end
    
    dmin = sqrt((uniqueX(end)-uniqueX(1))^2+(uniqueY(end)-uniqueY(1))^2);
    
    %     
%     if curvatures(k)<1e-10
%         straight=straight+1;
%     end
%     tempImage=double(binaryImage);
%     for i = 1:length(uniqueX)
%         tempImage(uniqueY(i),uniqueX(i)) = 0.5;
%     end
    
%     if nargin==4
%         subplot(121),imshow(tempImage);set(gcf,'units','normalized','outerposition',[0 0 1 1])
%         title(['Curvature: ',num2str(curvatures(k))])
%         pause(0.2)
%     end
end

curvatures=curvatures(halfWidth+1 : length(x) - halfWidth);
curvatures(abs(curvatures) > 1) = [];
media=mean(abs(curvatures));
numPoints = length(x);
% straightness = (dmin/dreal)^2;
straightness=1;
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');

end