function histogramprinting(tempROIs, roinumber, expectedCurvatureOne, expectedCurvatureTwo, firstname)

a = abs(tempROIs(roinumber(1)).individualCurvatures);
b = abs(tempROIs(roinumber(2)).individualCurvatures);
c = abs(tempROIs(roinumber(3)).individualCurvatures);
d = abs(tempROIs(roinumber(4)).individualCurvatures);
% e = abs(tempROIs(6).individualCurvatures);
% f = abs(tempROIs(7).individualCurvatures);
% g = abs(tempROIs(8).individualCurvatures);
% h = abs(tempROIs(9).individualCurvatures);
% 
% expectedCurvatureOne = 0.0307;
% expectedCurvatureTwo = 0.0029;
% nameA=tempROIs(2).imageName;
% nameB=tempROIs(3).imageName;
% nameC=tempROIs(4).imageName;
% nameD=tempROIs(5).imageName;
% nameE=tempROIs(6).imageName;
% nameF=tempROIs(7).imageName;
% nameG=tempROIs(8).imageName;
% nameH=tempROIs(9).imageName;
% a = a(a>0.01);
% b = b(b>0.01);
% c = c(c>0.01);
% d = d(d>0.01);

% Histogram printing
% colorA = [230 97 1]./255;
% colorB = [253 184 99]./255;
% colorC = [178 171 210]./255;

colorA = [215 48 39]./255;
colorB = [244 109 67]./255;
colorC = [253 174 97]./255;
colorD = [254 224 144]./255;
% colorE = [224 243 248]./255;
% colorF = [171 217 233]./255;
% colorG = [116 173 209]./255;
% colorH = [69 117 180]./255;
figure()
% set(gca,'color','black')
% set(gca,'color',[220 224 212]./255) % poster colour beige
% set(gcf,'color',[220 224 212]./255) % poster colour beige
% set(gcf,'color','none')
hold on
histogram(a,200,'facecolor',colorA,'facealpha',.7,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex','FontName','Calibri Light');ylabel('$Frequency$','Interpreter','latex');title('Individual Curvature Histogram','FontName','Calibri Light');
histogram(b,200,'facecolor',colorB,'facealpha',.7,'edgecolor','none','normalization','probability')
histogram(c,200,'facecolor',colorC,'facealpha',.7,'edgecolor','none','normalization','probability')
histogram(d,200,'facecolor',colorD,'facealpha',.7,'edgecolor','none','normalization','probability')
ylim([0 0.09]);
% xlim([0 0.1])
yL = get(gca,'YLim');
line([expectedCurvatureOne expectedCurvatureOne],yL,'Color','r','LineStyle','--');
if expectedCurvatureTwo > 0.0
    line([expectedCurvatureTwo expectedCurvatureTwo],yL,'Color','r','LineStyle','--');
end
% 
% histogram(e,200,'facecolor',colorE,'facealpha',.7,'edgecolor','none','normalization','probability')
% histogram(f,200,'facecolor',colorF,'facealpha',.7,'edgecolor','none','normalization','probability')
% histogram(g,200,'facecolor',colorG,'facealpha',.7,'edgecolor','none','normalization','probability')
% histogram(h,200,'facecolor',colorH,'facealpha',.7,'edgecolor','none','normalization','probability')
hold off
lgd = legend({firstname,'Waves','Non-rotated Crossed Waves','Rotated Crossed Waves'});
lgd.FontSize = 14;
lgd.FontName = 'Calibri Light';
% legend({'Ellipses','Waves','Non-rotated Crossed Waves','Rotated Crossed Waves'},'TextColor','w')
% legend({nameA,nameB,nameC,nameD,nameE,nameF,nameG,nameH})
% figure()
% subplot(221);set(gca,'color','black');histogram(a,200,'facecolor',colorA,'facealpha',.5,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex');ylabel('$Frequency$','Interpreter','latex');title('$Gaussian  Curvature  Squares$','Interpreter','latex');
% subplot(222);set(gca,'color','black');histogram(b,200,'facecolor',colorB,'facealpha',.5,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex');ylabel('$Frequency$','Interpreter','latex');title('$Gaussian  Curvature  Lines$','Interpreter','latex');
% subplot(223);set(gca,'color','black');histogram(c,200,'facecolor',colorC,'facealpha',.5,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex');ylabel('$Frequency$','Interpreter','latex');title('$Gaussian  Curvature  Non Rotated Crossed Lines$','Interpreter','latex');
% subplot(224);set(gca,'color','black');histogram(d,200,'facecolor',colorD,'facealpha',.5,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex');ylabel('$Frequency$','Interpreter','latex');title('$Gaussian  Curvature  Rotated Crossed Lines$','Interpreter','latex');
end