coloredSkeletonStack=zeros(0,0,0);
tolVec=10.^(-3:.1:4);

for tol=tolVec
    curvatureMatthiasLS;
    coloredSkeletonStack(:,:,end+1)=[coloredSkeleton;zeros(500-size(coloredSkeleton,1),size(coloredSkeleton,2))];
end
analyticSkel=[analyticSkel;zeros(500-size(analyticSkel,1),size(analyticSkel,2))];

figure()

for i=1:length(tolVec)
    imagesc(coloredSkeletonStack(:,:,i)+(coloredSkeletonStack(:,:,i)==0).*analyticSkel,[cMin cMax]);
    daspect([1 1 1]); colormap(gca,'jet'); %cheating: fixing the colormap to that of the analytic curve
    h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
    
    title({['a = ',num2str(a),', b = ',num2str(b),...
    ', ',num2str(tMin/pi),'\pi \leq \theta \leq ',num2str(tMax/pi),'\pi, tol = ',...
    num2str(tolVec(i)),', ',num2str(sections),' sections']})
    
    if i==1; gif(['lSpiral_Spline_TolAnim.gif'],'frame',gcf,'DelayTime',.5); else; gif; end
    pause(eps)
end
