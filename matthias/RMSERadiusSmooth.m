radiusRMS=[];
tolVec=10.^(-2:.1:3);
nKnotsVec=[1/3,1:8,251+1/3]*75;
nKnotsEffVec=zeros(size(nKnotsVec));

for nKnots=nKnotsVec
    for tol=tolVec
        curvatureMatthiasLS;
    end
    nKnotsEffVec(nKnots==nKnotsVec)=nKnotsEff;
end

radiusRMS=reshape(radiusRMS,length(tolVec),length(nKnotsVec));
tolX=sqrt(tolVec.'*ones(size(nKnotsVec)));

figure()

loglog(tolX,radiusRMS,'o-')

xlabel('Tolerance: estimate of \sigma at each knot')
ylabel('RMSE on the knots'' radii')
grid on

title({['a = ',num2str(a),', b = ',num2str(b),...
    ', ',num2str(tMin/pi),'\pi \leq \theta \leq ',num2str(tMax/pi),'\pi'],...
    ['R_{min} \approx ',num2str(rMinT),'{\Delta}p_0, R_{max} \approx ',num2str(rMaxT),...
    '{\Delta}p_0, R_{avg} \approx ',num2str(analyticalRadius),...
    '{\Delta}p_0, L \approx ',num2str(l(tMax)-l(tMin)),...
    '{\Delta}p_0'],...
    ['Gap_L \approx ',num2str(l(tD(selector(1)))-l(tMin)),'{\Delta}p_0, Gap_R \approx ',...
    num2str(l(tD(selector(end)))-l(tMin)),'{\Delta}p_0']})

legend([num2str(nKnotsEffVec.'),repmat(' knots',length(nKnotsEffVec),1)])