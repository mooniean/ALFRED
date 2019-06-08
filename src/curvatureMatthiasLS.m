clc; clear; close all;
figure('Units','normalized','OuterPosition',[0 0 1 1]);
offset=20;

% Logarithmic Spiral parameters

a=3;
b=.1;

% analytic expressions for the spiral parameterization

rt = @(t) a*exp(b*t);
xt = @(t) rt(t).*cos(t);
yt = @(t) rt(t).*sin(t);

% "analytic" mean radius
tMax=14*pi;

dxdt=matlabFunction(diff(sym(xt),1));
d2xdt2=matlabFunction(diff(sym(xt),2));
dydt=matlabFunction(diff(sym(yt),1));
d2ydt2=matlabFunction(diff(sym(yt),2));

kdldt= @(t) abs(dxdt(t).*d2ydt2(t)-dydt(t).*d2xdt2(t))./(dxdt(t).^2+dydt(t).^2);
dldt = @(t) sqrt(dxdt(t).^2+dydt(t).^2);
k = @(t) kdldt(t)./dldt(t);
dldtk = @(t) dldt(t)./k(t);

analyticalRadius=integral(dldtk,0,tMax)/integral(dldt,0,tMax)
radiusDensity=matlabFunction(sym(dldt)/diff(sym(k),1)*...
    sym(k)^2/integral(dldt,0,tMax));

%%%%%----------%%%%%

% fp-precision digitization

dtP=1e-3; tP=0:dtP:tMax; % the spacing affects what's the right p... with 1e-5 you should use p=0.(9) - some nines...
xP=xt(tP);
yP=yt(tP);
rP=1./k(tP);

%%%%%----------%%%%%

% fp-precision skeleton

minxP=min(xP); maxxP=max(xP);
minyP=min(yP); maxyP=max(yP);
analyticSkel=zeros(round(maxxP-minxP)+offset,round(maxyP-minyP)+offset);

colorX=round(xP-minxP)+1+offset/2;
colorY=round(yP-minyP)+1+offset/2;
for i=1:length(xP)
    analyticSkel(colorX(i)-1:colorX(i)+1,colorY(i)-1:colorY(i)+1)=...
        log10(rP(i));
end

subplot(141)
imagesc(analyticSkel); daspect([1 1 1]); colormap(gca,'jet')
h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
title(['[1] Logarithmic spiral: a = ',int2str(a),', b = ',num2str(b)]);
cMin=min(min(analyticSkel));cMax=max(max(analyticSkel));

%%%%%----------%%%%%

% spline the fp-precision curve

for p=1%10^-3.456%10.^[-0.001:-.2:-6]
    pause(eps)
%p=1e-3;%.9999999999;
ppX = csaps(tP, xP, p);
ppY = csaps(tP, yP, p);

xtP=@(t) fnval(ppX,t);
ytP=@(t) fnval(ppY,t);

dxdtP=@(t) fnval(fnder(ppX),t);
dydtP=@(t) fnval(fnder(ppY),t);

d2xdt2P=@(t) fnval(fnder(fnder(ppX)),t);
d2ydt2P=@(t) fnval(fnder(fnder(ppY)),t);

kdldtP=@(t) abs(dydtP(t).*d2xdt2P(t)-d2ydt2P(t).*dxdtP(t))./(dxdtP(t).^2+dydtP(t).^2);
dldtP=@(t) sqrt(dxdtP(t).^2+dydtP(t).^2);
kP = @(t) kdldtP(t)./dldtP(t);
dldtkP = @(t) dldtP(t)./kP(t);

radiusDensitySpline=@(t) (dldtP(t)./gradient(kP(t)).*kP(t).^2)...
    ./trapz(dldtP(t));

%%%%%----------%%%%%

% spline of the fp-precison curve's skeleton

tS=tP;
xS=xtP(tS);
yS=ytP(tS);
rS=1./kP(tS);

minxS=min(xS);
minyS=min(yS);
splineAnalyticSkel=zeros(size(analyticSkel));

colorX=round(xS-minxS)+1+offset/2;
colorY=round(yS-minyS)+1+offset/2;
for i=1:length(xS)  
    splineAnalyticSkel(colorX(i)-1:colorX(i)+1,colorY(i)-1:colorY(i)+1)=...
        log10(rS(i));
end
%assert(isequal(size(analyticSkel),size(splineAnalyticSkel)));

subplot(142),cla
imagesc(splineAnalyticSkel,[cMin cMax]); daspect([1 1 1]); colormap(gca,'jet'); %cheating: fixing the colormap to that of the analytic curve
h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
title(['[2] Smoothing spline of the analytic curve: p = ',num2str(p)]);

%%%%%----------%%%%%

% pixel-precision digitization

xD=round(xP-min(xP))+1;
yD=round(yP-min(yP))+1;

skeleton = zeros(max(xD),max(yD));
for i = 1:length(xD)
    skeleton(xD(i),yD(i)) = 1;
end

skeleton=bwmorph(skeleton,'skel','Inf');
[x,y]=find(skeleton);

% % give this to Fourier fit function
% 
% [~, ~, fourierRadius] = curvatureFourier(x,y,0);
% 
% % give this to Gaussian convolution function
% 
% sigma=18;
% [dTdt, dsdt, xt, yt, ~] = curvatureGauss(x,y,sigma);
% 
% dt=1;
% xt=round(xt{1});
% yt=round(yt{1});

% give this to Spline fit function

p=2e-5;
[dTdt, dsdt, xt, yt, dsdT] = curvatureSpline(x,y,p);

dt=dtP;t=20:dt:length(x)-20; % cheating: completely empirical
dTdt=dTdt{1}(t.');
dsdt=dsdt{1}(t.');
xt=round(xt{1}(t.'));
yt=round(yt{1}(t.'));

%%%%%----------%%%%%

% colorize curvature on the logarithmic spiral per constant-length-sector
sections=350;

circumference=dt*trapz(dsdt);
secLen=circumference/sections;

deltaS=dt*cumtrapz(dsdt);
breakpoints=ones(1,sections+1);
for i=1:sections-1; [~,breakpoints(i+1)]=min(abs(deltaS-secLen*i)); end
breakpoints(end)=length(dsdt);

meanGaussRadius=trapz(abs(dsdt./dTdt.*dsdt))/trapz(dsdt)
meanGaussRadii=zeros(1,sections);
actualSectionLens=zeros(1,sections);
coloredSkeleton=zeros(size(analyticSkel));
for i=1:sections
    meanGaussRadii(i)=trapz(abs(dsdt(breakpoints(i):breakpoints(i+1))...
        ./dTdt(breakpoints(i):breakpoints(i+1))...
        .*dsdt(breakpoints(i):breakpoints(i+1))))...
        /trapz(dsdt(breakpoints(i):breakpoints(i+1)));
    actualSectionLens(i)=dt*trapz(dsdt(breakpoints(i):breakpoints(i+1)));
    for j=breakpoints(i):breakpoints(i+1)-1
        colorX=xt(j);colorY=yt(j);
        coloredSkeleton(colorX-1:colorX+1,colorY-1:colorY+1)=log10(abs(meanGaussRadii(i)));
    end
end

%%%%%----------%%%%%

% plotting

subplot(143),cla
imagesc(coloredSkeleton,[cMin cMax]); daspect([1 1 1]); colormap(gca,'jet'); %cheating: fixing the colormap to that of the analytic curve
h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
% title(['Gaussian Filter: \sigma = ',num2str(sigma),...
%     '. Segments: ',int2str(sections)])
title(['[3] Smoothing spline of the pixel-digitized curve: p = ',num2str(p),...
    '. Segments: ',int2str(sections)])

subplot(144),cla
yLim=0.35;xlabel('Radius (px)'),ylabel('Relative frequency');

% setup histogram
binWidth=15;nBins=round((a*exp(b*tMax)*sqrt(1+b^2)-a*sqrt(1+b^2))/binWidth);
binWidth=(a*exp(b*tMax)*sqrt(1+b^2)-a*sqrt(1+b^2))/nBins;
meanGaussRadii(meanGaussRadii>a*exp(b*tMax)*sqrt(1+b^2))=a*exp(b*tMax)*sqrt(1+b^2); % cheating: everything to the right of last bin, goes to the last bin
meanGaussRadii(meanGaussRadii<a*sqrt(1+b^2))=a*sqrt(1+b^2); % cheating: everything to the left of first bin, goes to the first bin
H=histogram(meanGaussRadii,a*sqrt(1+b^2):binWidth:a*exp(b*tMax)*sqrt(1+b^2),'Normalization','probability');

% start plotting auxiliary data
hold on; line([a*sqrt(1+b^2) a*sqrt(1+b^2)],[0 yLim],'Color','r');
line([a*exp(b*tMax)*sqrt(1+b^2) a*exp(b*tMax)*sqrt(1+b^2)],[0 yLim],'Color','r');
plot(dldt(tP(2:end-1))./kdldt(tP(2:end-1)),abs(radiusDensity(tP(2:end-1)))*binWidth,'Linewidth',11);
line([analyticalRadius analyticalRadius],[0 yLim],'Color','k');

plot(dldtP(tS)./kdldtP(tS),abs(radiusDensitySpline(tS))*binWidth,'.');

line([meanGaussRadius meanGaussRadius],[0 yLim],'Color','g');
%line([fourierRadius fourierRadius],[0 yLim],'Color','b');

%%%%%----------%%%%%

%manually made and potentially more precise histogram

barPlot=zeros(1,nBins);
for i=1:length(meanGaussRadii)
    if (ceil((meanGaussRadii(i)-a*sqrt(1+b^2))/binWidth)>nBins); barPlot(nBins)=barPlot(nBins)+actualSectionLens(i)/circumference; continue; end
    if (ceil((meanGaussRadii(i)-a*sqrt(1+b^2))/binWidth)<1); barPlot(1)=barPlot(1)+actualSectionLens(i)/circumference; continue; end
    barPlot(ceil((meanGaussRadii(i)-a*sqrt(1+b^2))/binWidth))=...
        barPlot(ceil((meanGaussRadii(i)-a*sqrt(1+b^2))/binWidth))+...
        actualSectionLens(i)/circumference;
end

%%%%%----------%%%%%


% legend({'Gaussian section radius','Minimum radius','Maximum radius',...
%     'Gaussian radius (average)',...%'Fourier radius',...
%     '"Analytic" (average)','Probability density (\times binWidth)'})
plots=get(gca, 'Children');
legend(plots([6:-1:1,7]),...
    {'Minimum radius [1]','Maximum radius[1]',...
    'Analytic probability density (\times binWidth) [1]','Analytic average [1]',...
    'Spine probability density (\times binWidth) [2]',...
    'Spline radius (average) [3]','Spline section radius [3]',...%'Fourier radius'
    })
xlim([0,250]);ylim([0 yLim]);
%if p>.9; gif(['lSpiral_Spline_',int2str(sections),'.gif'],'frame',gcf,'DelayTime',2); else; gif; end
end