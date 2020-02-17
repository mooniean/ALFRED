close all force
clc 
clear

load('results/20191101-134043_TestImages_Squares_Original.tif.mat')

gray = [191 191 191]./255;
orange = [247 150 70]./255;
blue = [89 146 217]./255;

%% LogSpiral


a=3;b=.1;theta = 8*pi:1e-3:14*pi;s = a.*sqrt(1+b^2)*exp(b.*theta)/b;k = exp(-b.*theta)/a.*sqrt(1+b^2);
arclength = 1e-3*cumtrapz(tempROIs(2).arcLength);
figure(); 
hold on
plot(abs(s-s(end)),k,'.','Color',gray)
plot(arclength,tempROIs(2).individualCurvatures,'.','Color',blue)
set(gca,'color','none')
hold off

% export_fig -transparent logspiralResult.png

%% Ellipse
offset=20;

% ellipse's axes parameters

a=400;
b=200;
alpha=45*pi/180;

% analytic expressions for the ellipse parameterization

xt = @(t) a.*cos(t)*cos(alpha)-b.*sin(t)*sin(alpha);
yt = @(t) a.*cos(t)*sin(alpha)+b.*sin(t)*cos(alpha);
EI2= @(phi,k) integral(@(t) sqrt(1-k^2*sin(t).^2),0,phi);

% "analytic" mean radius
tMin=0; tMax=2*pi;

dxdt=matlabFunction(diff(sym(xt),1));
d2xdt2=matlabFunction(diff(sym(xt),2));
dydt=matlabFunction(diff(sym(yt),1));
d2ydt2=matlabFunction(diff(sym(yt),2));

kdldt= @(t) abs(dxdt(t).*d2ydt2(t)-dydt(t).*d2xdt2(t))./(dxdt(t).^2+dydt(t).^2);
dldt = @(t) sqrt(dxdt(t).^2+dydt(t).^2);
l = @(t) b*EI2(t,sqrt(1-(a/b)^2));
k = @(t) kdldt(t)./dldt(t);
dldtk = @(t) dldt(t)./k(t);

analyticalRadius=integral(dldtk,tMin,tMax)/integral(dldt,tMin,tMax);
radiusDensity=matlabFunction(sym(dldt)/diff(sym(k),1)*...
    sym(k)^2/integral(dldt,tMin,tMax));

figure(); hold on
arclength = 1e-3*cumtrapz(tempROIs(3).arcLenangth);
tempValue = 0.407407346410207;

plot(abs(arrayfun(l,(tempValue:1e-3:2*pi+tempValue))-l(2*pi+tempValue)),k(tempValue:1e-3:2*pi+tempValue),'.','Color',gray);
plot(arclength,tempROIs(3).individualCurvatures,'.','Color',orange)
set(gca,'color','none')

hold off

% export_fig -transparent ellipseResult.png
