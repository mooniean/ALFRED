% Ellipse curvature

% clc; clear; close all;




colorA = [215 48 39]./255;
colorB = [191 191 191]./255;

% a = abs(tempROIs(3).individualCurvatures);
figure()
hold on
% histogram(a,200,'facecolor',colorA,'facealpha',.7,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex','FontName','Calibri Light');ylabel('$Frequency$','Interpreter','latex');title('Individual Curvature Histogram','FontName','Calibri Light');
% xlim([0 0.04])
% ylim([0 0.04])
% yL = get(gca,'YLim');


a=5;
b=2;

x=0:1e-5:a;

dydx=-b*x/a^2./sqrt(1-(x/a).^2);
d2ydx2=-b/a^2./sqrt(1-(x/a).^2).^3;

k=abs(d2ydx2)./(1+dydx.^2).^(3/2);


histogram(k,200,'facecolor',colorB,'facealpha',.8,'edgecolor','none','normalization','probability');xlabel('$Curvature (px^{-1})$','Interpreter','latex','FontName','Calibri Light');ylabel('$Frequency$','Interpreter','latex');title('Individual Curvature Histogram','FontName','Calibri Light');

% curveOne = b/a^2;
% curveTwo = a/b^2;
% line([curveOne curveOne],yL,'Color','r','LineStyle','--');
% line([curveTwo curveTwo],yL,'Color','r','LineStyle','--');
