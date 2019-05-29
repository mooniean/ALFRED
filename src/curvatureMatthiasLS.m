clc; clear; close all;

% Logarithmic Spiral parameters

a=3;
b=.1;

% analytic expressions for the spiral parameterization

rt = @(t) a*exp(b*t);
xt = @(t) rt(t).*cos(t);
yt = @(t) rt(t).*sin(t);

% "analytic" mean radius
tMax=14*pi;

x1=matlabFunction(diff(sym(xt),1));
x2=matlabFunction(diff(sym(xt),2));
y1=matlabFunction(diff(sym(yt),1));
y2=matlabFunction(diff(sym(yt),2));

kdldt= @(t) abs(x1(t).*y2(t)-y1(t).*x2(t))./(x1(t).^2+y1(t).^2);
dldt = @(t) (x1(t).^2+y1(t).^2).^(1/2);
k = @(t) kdldt(t)./dldt(t);
dldtk = @(t) dldt(t)./k(t);

analyticalRadius=integral(dldtk,0,tMax)/integral(dldt,0,tMax)
radiusDensity=matlabFunction(sym(dldt)/diff(sym(k),1)*...
    sym(k)^2/integral(dldt,0,tMax));

%%%%%----------%%%%%

% fp-precision digitization

dtP=1e-3; tP=0:dtP:tMax;
xP=xt(tP);
yP=yt(tP);

p=1;
ppX = csaps( tP, xP.', p);
ppY = csaps( tP, yP.', p);

coeffvaluesX=ppX.coefs;
coeffvaluesY=ppY.coefs;

xtP=@(t) coeffvaluesX(floor(t/dtP)+1,1).*(t-dtP*floor(t/dtP)).^3+coeffvaluesX(floor(t/dtP)+1,2).*(t-dtP*floor(t/dtP)).^2+...
    coeffvaluesX(floor(t/dtP)+1,3).*(t-dtP*floor(t/dtP))+coeffvaluesX(floor(t/dtP)+1,4);
ytP=@(t) coeffvaluesY(floor(t/dtP)+1,1).*(t-dtP*floor(t/dtP)).^3+coeffvaluesY(floor(t/dtP)+1,2).*(t-dtP*floor(t/dtP)).^2+...
    coeffvaluesY(floor(t/dtP)+1,3).*(t-dtP*floor(t/dtP))+coeffvaluesY(floor(t/dtP)+1,4);

dxdt=@(t) 3*coeffvaluesX(floor(t/dtP)+1,1).*(t-dtP*floor(t/dtP)).^2+2*coeffvaluesX(floor(t/dtP)+1,2).*(t-dtP*floor(t/dtP))+...
    coeffvaluesX(floor(t/dtP)+1,3);
dydt=@(t) 3*coeffvaluesY(floor(t/dtP)+1,1).*(t-dtP*floor(t/dtP)).^2+2*coeffvaluesY(floor(t/dtP)+1,2).*(t-dtP*floor(t/dtP))+...
    coeffvaluesY(floor(t/dtP)+1,3);

d2xdt2=@(t) 6*coeffvaluesX(floor(t/dtP)+1,1).*(t-dtP*floor(t/dtP))+2*coeffvaluesX(floor(t/dtP)+1,2);
d2ydt2=@(t) 6*coeffvaluesY(floor(t/dtP)+1,1).*(t-dtP*floor(t/dtP))+2*coeffvaluesY(floor(t/dtP)+1,2);

kdldtP=@(t) abs(dydt(t).*d2xdt2(t)-d2ydt2(t).*dxdt(t))./(dxdt(t).^2+dydt(t).^2);
dldtP=@(t) sqrt(dxdt(t).^2+dydt(t).^2);
kP = @(t) kdldtP(t)./dldtP(t);
dldtkP = @(t) dldtP(t)./kP(t);

radiusDensitySpline=@(t) (dldt(t)./gradient(k(t)).*k(t).^2)...
    ./trapz(dldt(t));

%%%%%----------%%%%%

% pixel-precision digitization

tD=0:1e-5:tMax;
xD=round(xt(tD));
yD=round(yt(tD));

xD=xD-min(xD)+1;
yD=yD-min(yD)+1;

skeleton = zeros(max(xD),max(yD));
for i = 1:length(xD)
    skeleton(xD(i),yD(i)) = 1;
end

skeleton=bwmorph(skeleton,'skel','Inf');
[x,y]=find(skeleton);

% give this to Fourier fit function

[~, ~, fourierRadius] = curvatureFourier(x,y,0);

% give this to Gaussian convolution function

sigma=18;
[dTdt, dsdt, xt, yt, ~] = curvatureGauss(x,y,sigma);

dt=1;
xt{1}=round(xt{1});
yt{1}=round(yt{1});

% give this to Spline fit function

p=2e-5;
[dTdt, dsdt, xt, yt, dsdT] = curvatureSpline(x,y,p);

dt=1e-2;t=2*sigma+1:dt:length(x)-2*sigma;
dTdt{1}=dTdt{1}(t.');
dsdt{1}=dsdt{1}(t.');
xt{1}=round(xt{1}(t.'));
yt{1}=round(yt{1}(t.'));

%%%%%----------%%%%%

% colorize curvature on the ellipse per constant-length-sector
sections=350;

circumference=dt*trapz(dsdt{1});
secLen=circumference/sections;

deltaS=dt*cumtrapz(dsdt{1});
breakpoints=ones(1,sections+1);
for i=1:sections-1; [~,breakpoints(i+1)]=min(abs(deltaS-secLen*i)); end
breakpoints(end)=length(dsdt{1});

meanGaussRadius=trapz(abs(dsdt{1}./dTdt{1}.*dsdt{1}))/trapz(dsdt{1})
meanGaussRadii=zeros(1,sections);
actualSectionLens=zeros(1,sections);
offset=20; coloredSkeleton=zeros(max(xt{1})+offset,max(yt{1})+offset);
for i=1:sections
    meanGaussRadii(i)=trapz(abs(dsdt{1}(breakpoints(i):breakpoints(i+1))...
        ./dTdt{1}(breakpoints(i):breakpoints(i+1))...
        .*dsdt{1}(breakpoints(i):breakpoints(i+1))))...
        /trapz(dsdt{1}(breakpoints(i):breakpoints(i+1)));
    actualSectionLens(i)=dt*trapz(dsdt{1}(breakpoints(i):breakpoints(i+1)));
    for j=breakpoints(i):breakpoints(i+1)-1
        colorX=xt{1}(j)+offset/2;colorY=yt{1}(j)+offset/2;
        coloredSkeleton(colorX-1:colorX+1,colorY-1:colorY+1)=log10(abs(meanGaussRadii(i)));
    end
end
origSkeleton=[zeros(offset,3/2*offset+size(skeleton,2));...
    zeros(size(skeleton,1),offset),skeleton,zeros(size(skeleton,1),offset/2);...
    zeros(offset/2,3/2*offset+size(skeleton,2))];
coloredSkeleton=[coloredSkeleton;zeros(size(origSkeleton,1)-size(coloredSkeleton,1),size(coloredSkeleton,2))];
coloredSkeleton=[coloredSkeleton,zeros(size(origSkeleton,1),size(origSkeleton,2)-size(coloredSkeleton,2))];

%%%%%----------%%%%%

% plotting
figure('Units','normalized','OuterPosition',[0 0 1 1]);

subplot(131)
imagesc(origSkeleton); daspect([1 1 1]); colormap(gca,'gray');
h=colorbar('southoutside');
title(['Logarithmic spiral: a = ',int2str(a),', b = ',num2str(b)]);

subplot(132)
imagesc(coloredSkeleton); daspect([1 1 1]); colormap(gca,'jet')
h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
% title(['Gaussian Filter: \sigma = ',num2str(sigma),...
%     '. Segments: ',int2str(sections)])
title(['Smoothing spline: p = ',num2str(p),...
    '. Segments: ',int2str(sections)])

subplot(133),yLim=0.35;ylim([0 yLim]);xlabel('Radius (px)'),ylabel('Relative frequency');
binWidth=15;nBins=round((a*exp(b*tMax)*sqrt(1+b^2)-a*sqrt(1+b^2))/binWidth);
binWidth=(a*exp(b*tMax)*sqrt(1+b^2)-a*sqrt(1+b^2))/nBins;
meanGaussRadii(meanGaussRadii>a*exp(b*tMax)*sqrt(1+b^2))=a*exp(b*tMax)*sqrt(1+b^2);
meanGaussRadii(meanGaussRadii<a*sqrt(1+b^2))=a*sqrt(1+b^2);
H=histogram(meanGaussRadii,a*sqrt(1+b^2):binWidth:a*exp(b*tMax)*sqrt(1+b^2),'Normalization','probability');

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

hold on; line([a*sqrt(1+b^2) a*sqrt(1+b^2)],[0 yLim],'Color','r'); line([a*exp(b*tMax)*sqrt(1+b^2) a*exp(b*tMax)*sqrt(1+b^2)],[0 yLim],'Color','r');
line([meanGaussRadius meanGaussRadius],[0 yLim],'Color','g');
%line([fourierRadius fourierRadius],[0 yLim],'Color','b');
line([analyticalRadius analyticalRadius],[0 yLim],'Color','k');
plot(dldt(tD(2:end-1))./kdldt(tD(2:end-1)),abs(radiusDensity(tD(2:end-1)))*binWidth,'Linewidth',5);
plot(dldtP((dtP:dtP:tMax-5*dtP).')./kdldtP((dtP:dtP:tMax-5*dtP).'),abs(radiusDensitySpline((dtP:dtP:tMax-5*dtP).'))*binWidth,'Linewidth',2); hold off;
% legend({'Gaussian section radius','Minimum radius','Maximum radius',...
%     'Gaussian radius (average)',...%'Fourier radius',...
%     '"Analytic" (average)','Probability density (\times binWidth)'})
legend({'Spline section radius','Minimum radius','Maximum radius',...
    'Spline radius (average)',...%'Fourier radius',...
    '"Analytic" (average)','Probability density (\times binWidth)','Spine probability density (\times binWidth)'})

if ~0; gif(['lSpiral_Spline_',int2str(sections),'.gif'],'frame',gcf,'DelayTime',2); else; gif; end

%%%%%----------%%%%%

