function allK = curvature(origImage)
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

% imshow(origImage)

allBoundaries=bwboundaries(origImage);
avgK=[];avgR=[];
allK=[];allR=[];

for d=3:21

    ks=[];s=[];
    
    for b=1:length(allBoundaries)
        
        if size(allBoundaries{b},1)<d+1; continue; end  
        x=allBoundaries{b}(:,1).';y=allBoundaries{b}(:,2).';t=1:length(x);
        
        [~,~,~,weightedK,ds]=polyXY(t,x,y,d);
%         diameter = max(max(resX))-min(min(resX))
        ks(end+1)=sum(weightedK)-1/2*(weightedK(1)+weightedK(end));
        s(end+1)=sum(ds);
        
%         %%%% This can be deleted... useful for debugging purposes only
%         bImage=zeros(size(origImage)); rImage=bImage;
%         for i=t
%             bImage(x(i),y(i))=i/t(end);
%             j=find(u==i,1);
%             if j; rImage(round(resX(i)),round(resY(i)))=weightedK(j)+0.5; end
%         end
%         figure(1);imshow(origImage);
%         figure(2);imshow(bImage);
%         figure(3);imshow(rImage);
%         figure(4),subplot(121),plot(weightedK,'o-'),subplot(122),plot(ds,'-o')
        %%%
        
    end
    
    %mean(ks./s)
    %mean(s./ks)
    allK=[allK ks./s];
    allR=[allR s./ks];
    avgK(end+1)=sum(ks)/sum(s);
    avgR(end+1)=sum(s)/sum(ks);
    
end
% % 
% mean(avgK)  
% mean(avgR)
% 
end