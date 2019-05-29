function [resX,resY,u,weightedK,ds]=polyXY(t,x,y,d)

if nargin<4; d = 9; end
allTerms = d:-1:0;

[Px,Sx] = polyfit(t,x,d);
resX=Px*(repmat(t,length(Px),1).^repmat(allTerms.',1,length(t)));
rmseX=norm(x-resX)/sqrt(Sx.df);

[Py,Sy] = polyfit(t,y,d);
resY=Py*(repmat(t,length(Py),1).^repmat(allTerms.',1,length(t)));
rmseY=norm(y-resY)/sqrt(Sy.df);

if (rmseX+rmseY)>3 && length(t)>=2*(d+1) %test for bad fit in the previous steps and guarantees the fits are possible with polynomials of degree d
    br=floor(length(t)/2);
    [resX1,resY1,u1,weightedK1,ds1]=polyXY(t(1:br),x(1:br),y(1:br),d);
    [resX2,resY2,u2,weightedK2,ds2]=polyXY(t(br+1:end),x(br+1:end),y(br+1:end),d);
    resX=[resX1,resX2];
    resY=[resY1,resY2];
    u=[u1,u2];
    weightedK=[weightedK1,weightedK2];
    ds=[ds1,ds2];
    return;
end

firstDerX=(Px.*allTerms)*(repmat(t,length(Px),1).^repmat(allTerms.'-1,1,length(t)));
secondDerX=(Px.*allTerms.*(allTerms-1))*(repmat(t,length(Px),1).^repmat(allTerms.'-2,1,length(t)));

firstDerY=(Py.*allTerms)*(repmat(t,length(Py),1).^repmat(allTerms.'-1,1,length(t)));
secondDerY=(Py.*allTerms.*(allTerms-1))*(repmat(t,length(Py),1).^repmat(allTerms.'-2,1,length(t)));

ds=sqrt(firstDerX.^2+firstDerY.^2);
weightedK=abs(firstDerX.*secondDerY-firstDerY.*secondDerX)./(firstDerX.^2+firstDerY.^2).^(3/2).*ds;

[maxWK,br]=max(weightedK);
if maxWK>3*mean(weightedK) && br>=d+1 && (length(t)-br)>=d+1 %tries to predict change in MT w/ large curvature and guarantees the fits are possible with polynomials of degree d
    [resX1,resY1,u1,weightedK1,ds1]=polyXY(t(1:br),x(1:br),y(1:br),d);
    [resX2,resY2,u2,weightedK2,ds2]=polyXY(t(br+1:end),x(br+1:end),y(br+1:end),d);
    resX=[resX1,resX2];
    resY=[resY1,resY2];
    u=[u1,u2];
    weightedK=[weightedK1,weightedK2];
    ds=[ds1,ds2];
    return;
end

%%%% clean extremes
n1=round(length(t)*0.1);
n2=round(length(t)*0.2); %keep the middle 70% of the weightedK values
[~,I]=sort(weightedK);
newI=sort(I(n1+1:end-n2));

u=t(newI);
ds=ds(newI);
weightedK=weightedK(newI);
%%%%

end