close all; clear; clc

alphaVec=(0:5:90)*pi/180;
nKnotsVec=[1/3,1:8,251+1/3]*75;
nKnotsEffVec=zeros(size(nKnotsVec));

figure(1)

for nKnots=nKnotsVec
    for alpha=0
        curvatureMatthiasLS;
        figure(1)%,cla;
        plot(cumtrapz(dsdt)*dt+l(tD(2))-l(tP(1)),dsdT);hold on;
        %plot(l(tP)-l(tP(1)),1./k(tP),'LineWidth',4,'Color','k')
        
        xlabel('s [{\Delta}p_0]')
        xlim([0,2.1e3])
        ylabel('R(s) [{\Delta}p_0]')
        ylim([0,300])
        grid on
        
        rMinT=a*exp(b*tMin)*sqrt(1+b^2);
        rMaxT=a*exp(b*tMax)*sqrt(1+b^2);
        
        title({['a = ',num2str(a),', b = ',num2str(b),...
            ', ',num2str(tMin/pi),'\pi \leq \theta \leq ',num2str(tMax/pi),'\pi'],...
            ['R_{min} \approx ',num2str(rMinT),'{\Delta}p_0, R_{max} \approx ',num2str(rMaxT),...
            '{\Delta}p_0, R_{avg} \approx ',num2str(analyticalRadius),...
            '{\Delta}p_0, L \approx ',num2str(l(tMax)-l(tMin)),...
            '{\Delta}p_0'],['\alpha = ',num2str(alpha*180/pi),char(176)]})
        
        %if alpha==0; gif('lSpiral_Spline_RadiusAnim.gif','frame',gcf,'DelayTime',.5); else; gif; end
    end
    nKnotsEffVec(nKnots==nKnotsVec)=nKnotsEff;
end

% legend([string([num2str(alphaVec.'*180/pi),repmat(char(176),length(alphaVec),1)]);'Analytical radius'],...
%      'Location','NorthWest')

plot(l(tP)-l(tP(1)),1./k(tP),'LineWidth',4,'Color','k')

legend([string([num2str(nKnotsEffVec.'),repmat(' knots',length(nKnotsEffVec),1)]);'Analytical radius'],...
    'Location','NorthWest')