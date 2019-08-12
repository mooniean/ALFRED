radiusRMS=[];
fVec=10.^(0:.1:6);
nKnotsVec=(1:8)*75;

for nKnots=nKnotsVec
    for f=fVec
        curvatureMatthiasLS;
    end
end

radiusRMS=reshape(radiusRMS,length(fVec),length(nKnotsVec));
dpdx=1./(fVec.'*(l(tMax)-l(tMin))./(nKnotsVec-1));

figure()

loglog(dpdx,radiusRMS,'o-')

xlabel('{\Delta}p/{\Delta}x')
ylabel('RMSE on the knots'' radii')
grid on

title({['a = ',num2str(a),', b = ',num2str(b),...
    ', ',num2str(tMin/pi),'\pi \leq \theta \leq ',num2str(tMax/pi),'\pi, tol = ',...
    num2str(tol),', ',num2str(sections),' sections'],...
    ['R_{min} \approx ',num2str(rMinT),'{\Delta}p_0, R_{max} \approx ',num2str(rMaxT),...
    '{\Delta}p_0, R_{avg} \approx ',num2str(analyticalRadius),...
    '{\Delta}p_0, L \approx ',num2str(l(tMax)-l(tMin)),...
    '{\Delta}p_0']})

legend([num2str(nKnotsVec.'),repmat(' knots',length(nKnotsVec),1)])