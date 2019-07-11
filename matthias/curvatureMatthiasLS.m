clc; clearvars -except f radiusRMS; close all;
figure('Units','normalized','OuterPosition',[0 0 1 1]);
offset=20;

% Logarithmic Spiral parameters

a=3;
b=.1;
alpha=0*pi/180;

% analytic expressions for the spiral parameterization

rt =  @(t) a*exp(b*t);
xtU = @(t) rt(t).*cos(t);
ytU = @(t) rt(t).*sin(t);

xt = @(t) xtU(t)*cos(alpha)-ytU(t)*sin(alpha);
yt = @(t) xtU(t)*sin(alpha)+ytU(t)*cos(alpha);

% "analytic" mean radius
tMin=8*pi; tMax=14*pi;

dxdt=matlabFunction(diff(sym(xt),1));
d2xdt2=matlabFunction(diff(sym(xt),2));
dydt=matlabFunction(diff(sym(yt),1));
d2ydt2=matlabFunction(diff(sym(yt),2));

kdldt= @(t) abs(dxdt(t).*d2ydt2(t)-dydt(t).*d2xdt2(t))./(dxdt(t).^2+dydt(t).^2);
dldt = @(t) sqrt(dxdt(t).^2+dydt(t).^2);
l = @(t) rt(t)*sqrt(1+b^2)/b;
k = @(t) kdldt(t)./dldt(t);
dldtk = @(t) dldt(t)./k(t);

analyticalRadius=integral(dldtk,tMin,tMax)/integral(dldt,tMin,tMax)
radiusDensity=matlabFunction(sym(dldt)/diff(sym(k),1)*...
    sym(k)^2/integral(dldt,tMin,tMax));

%%%%%----------%%%%%

% fp-precision digitization

dtP=1e-3; tP=tMin:dtP:tMax;
tP1=linspace(tMin,tMax,75);
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
title({'[1] Logarithmic spiral:',['a = ',num2str(a),', b = ',num2str(b),...
    ', ',num2str(tMin/pi),'\pi \leq \theta \leq ',num2str(tMax/pi),'\pi']});
cMin=min(min(analyticSkel));cMax=max(max(analyticSkel));

%%%%%----------%%%%%

% spline the fp-precision curve

p=1;
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

%errX=sqrt(mean((xP-xS).^2))
%errY=sqrt(mean((yP-yS).^2))

splineAnalyticSkel=zeros(size(analyticSkel));

colorX=round(xS-minxP)+1+offset/2;
colorY=round(yS-minyP)+1+offset/2;
for i=1:length(xS)
    splineAnalyticSkel(colorX(i)-1:colorX(i)+1,colorY(i)-1:colorY(i)+1)=...
        log10(rS(i));
end
%assert(isequal(size(analyticSkel),size(splineAnalyticSkel)));

subplot(142),cla
imagesc(splineAnalyticSkel,[cMin cMax]); daspect([1 1 1]); colormap(gca,'jet'); %cheating: fixing the colormap to that of the analytic curve
h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
title({'[2] Smoothing spline of the analytic curve:',['p = ',num2str(p)]});

%%%%%----------%%%%%

% pixel-precision digitization 
nKnots=150; knots=zeros(1,nKnots);
knots(1)=tMin; knots(end)=tMax;
deltaL=(l(tMax)-l(tMin))/(nKnots-1);

for i=1:nKnots-2; knots(i+1)=fzero(@(t) l(t)-l(tMin)-deltaL*i,tMin); end

%f=1e4; %resolution increase factor
xD=round(f*(xt(knots)-minxP+1));  rms(xD-(xt(knots)-minxP)-1)
yD=round(f*(yt(knots)-minyP+1));  rms(yD-(yt(knots)-minyP)-1)

x=xD.';
y=yD.';
tD=knots;

% give this to Spline fit function (imported here, because I would have to mess up the function way too much)

%I'm using a predefined tolerance for now, I can perhaps allow the
%caller to set it... we'll see
tol=0;
L=nKnots; N=L; t=1:N;
ppX = spaps( t, x/f, tol*L);
ppY = spaps( t, y/f, tol*L);

xtR=@(t) fnval(ppX,t);
ytR=@(t) fnval(ppY,t);

dxdtR=@(t) fnval(fnder(ppX),t);
dydtR=@(t) fnval(fnder(ppY),t);

d2xdt2R=@(t) fnval(fnder(fnder(ppX)),t);
d2ydt2R=@(t) fnval(fnder(fnder(ppY)),t);

kdldtR=@(t) abs(dydtR(t).*d2xdt2R(t)-d2ydt2R(t).*dxdtR(t))./(dxdtR(t).^2+dydtR(t).^2);
dldtR=@(t) sqrt(dxdtR(t).^2+dydtR(t).^2);
kR = @(t) kdldtR(t)./dldtR(t);
dldtkR = @(t) dldtR(t)./kR(t);

dt=dtP;t=2:dt:N-1; % cheating: completely empirical; not so much now, just taking out the first and last splines
dTdt=kdldtR(t);
dsdt=dldtR(t);%dsdt=dldt(tP);
xt=round(xtR(t));%xt=round(xP-minxP)+10;
yt=round(ytR(t));%yt=round(yP-minyP)+10;
dsdT=1./kR(t);%dsdT=rP;

radiusRMS(end+1)=rms(1./kR(2:N-1)-1./k(knots(2:N-1)))

%%%%%----------%%%%%

% colorize curvature on the logarithmic spiral per constant-length-sector
sections=450; % try 35, 75, 350, 450, 750 with step=80

circumference=dt*trapz(dsdt);
secLen=circumference/sections;

deltaS=dt*cumtrapz(dsdt);
breakpoints=ones(1,sections+1);
for i=1:sections-1; [~,breakpoints(i+1)]=min(abs(deltaS-secLen*i)); end
breakpoints(end)=length(dsdt);

meanEstimatedRadius=trapz(abs(dsdT.*dsdt))/trapz(dsdt)
err=100*abs(meanEstimatedRadius-analyticalRadius)/analyticalRadius

meanEstimatedRadii=zeros(1,sections);
actualSectionLens=zeros(1,sections);
coloredSkeleton=zeros(size(analyticSkel));
for i=1:sections
    meanEstimatedRadii(i)=trapz(abs(dsdT(breakpoints(i):breakpoints(i+1))...
        .*dsdt(breakpoints(i):breakpoints(i+1))))...
        /trapz(dsdt(breakpoints(i):breakpoints(i+1)));
    actualSectionLens(i)=dt*trapz(dsdt(breakpoints(i):breakpoints(i+1)));
    for j=breakpoints(i):breakpoints(i+1)-1
        colorX=xt(j);colorY=yt(j);
        coloredSkeleton(offset/2+(colorX-1:colorX+1),offset/2+(colorY-1:colorY+1))=log10(abs(meanEstimatedRadii(i)));
        if mod(j-1+19/dt,25/dt)==0
            %coloredSkeleton(colorX-3:colorX+3,colorY-3:colorY+3)=1;
        end
    end
end

%%%%%----------%%%%%

% plotting

subplot(143),cla
imagesc(coloredSkeleton,[cMin cMax]); daspect([1 1 1]); colormap(gca,'jet'); %cheating: fixing the colormap to that of the analytic curve
h=colorbar('southoutside'); h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
% title(['Gaussian Filter: \sigma = ',num2str(sigma),...
%     '. Segments: ',int2str(sections)])
title({'[3] Smoothing spline of the pixel-digitized curve:',['Error tol. = ',num2str(tol),...
    ', segments: ',int2str(sections)]})

subplot(144),cla
xlabel('Radius (px)'),ylabel('Relative frequency');

% setup histogram
rMinT=a*exp(b*tMin)*sqrt(1+b^2);
rMaxT=a*exp(b*tMax)*sqrt(1+b^2);
rMinE=a*exp(b*min(tD))*sqrt(1+b^2);
rMaxE=a*exp(b*max(tD))*sqrt(1+b^2);

nBins=10;
binWidth=(rMaxE-rMinE)/nBins;
%meanEstimatedRadii(meanEstimatedRadii>rMaxE)=rMaxE; % cheating: everything to the right of last bin, goes to the last bin
%meanEstimatedRadii(meanEstimatedRadii<rMinE)=rMinE; % cheating: everything to the left of first bin, goes to the first bin
H=histogram(meanEstimatedRadii,rMinE:binWidth:rMaxE,'Normalization','probability');

% start plotting auxiliary data
yLim=2*max(H.Values);
hold on; line([rMinT rMinT],[0 yLim],'Color','r');
line([rMaxT rMaxT],[0 yLim],'Color','r');
plot(dldt(tP(2:end-1))./kdldt(tP(2:end-1)),abs(radiusDensity(tP(2:end-1)))*binWidth,'Linewidth',13);
line([analyticalRadius analyticalRadius],[0 yLim],'Color','k');

plot(dldtP(tS)./kdldtP(tS),abs(radiusDensitySpline(tS))*binWidth,'.');

line([meanEstimatedRadius meanEstimatedRadius],[0 yLim],'Color','g');
%line([fourierRadius fourierRadius],[0 yLim],'Color','b');

%%%%%----------%%%%%

%manually made and potentially more precise histogram

barPlot=zeros(1,nBins);
for i=1:length(meanEstimatedRadii)
    if ((meanEstimatedRadii(i)-rMinE)/binWidth>nBins); barPlot(nBins)=barPlot(nBins)+actualSectionLens(i)/circumference; continue; end
    if ((meanEstimatedRadii(i)-rMinE)/binWidth<1); barPlot(1)=barPlot(1)+actualSectionLens(i)/circumference; continue; end
    barPlot(ceil((meanEstimatedRadii(i)-rMinE)/binWidth))=...
        barPlot(ceil((meanEstimatedRadii(i)-rMinE)/binWidth))+...
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
    ['Spline radius (average) [3] (',num2str(err,'%.2f'),'%)'],'Spline section radius [3]',...%'Fourier radius'
    })
xlim([rMinT-0.1*(rMaxT-rMinT) rMaxT+0.1*(rMaxT-rMinT)]);ylim([0 yLim]);
%if f>0; gif(['lSpiral_Spline_',int2str(f),'_',num2str(nKnots),'.gif'],'frame',gcf,'DelayTime',2); else; gif; end