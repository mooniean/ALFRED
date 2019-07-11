radiusRMS=[];
fVec=10.^(0:.05:6);

for f=fVec
    curvatureMatthiasLS;
end

figure()
loglog(fVec,radiusRMS,'o')

xlabel('Resolution artificial increase factor')
ylabel('RMSE on the knots'' radii')
grid on

title({['a = ',num2str(a),', b = ',num2str(b),...
    ', ',num2str(tMin/pi),'\pi \leq \theta \leq ',num2str(tMax/pi),'\pi, ',...
    int2str(nKnots),' knots, tol = ',num2str(tol),', ',num2str(sections),' sections'],...
    ['R_{min} \approx ',num2str(rMinT),'px, R_{max} \approx ',num2str(rMaxT),...
    'px, R_{avg} \approx ',num2str(analyticalRadius),...
    'px, L \approx ',num2str(l(tMax)-l(tMin)),...
    'px, {\Delta}L \approx ',num2str(deltaL),'px']})