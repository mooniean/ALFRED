function [N, xt, yt, dsdT, dsdt, dTdt] = curvatureGauss(skeleton,sigma)

%find boundaries for objects not holes (better performance)
% boundaries=bwboundaries(skeleton,'noholes');
boundaries=bwboundaries(skeleton);
lenBoundaries=length(boundaries);

%storage for x(t) and y(t), where the curvature computation is valid
xt=cell(lenBoundaries,1);
yt=cell(lenBoundaries,1);

%storage for derivatives of T (the curve's tangent vector) and s (the curve
%'s arc length) for each boundary
dTdt=cell(lenBoundaries,1);
dsdt=cell(lenBoundaries,1);

%storage for radius at t
dsdT = cell(lenBoundaries,1);

%storage for parametrization length
N = cell(lenBoundaries,1);

for b=1:lenBoundaries
    
    x=boundaries{b}(:,1).';y=boundaries{b}(:,2).';
    
    %for line-style objects, the boundaries will include a back and forth
    %trace, which produces a mirrored path: here we eliminate one of them
    start=1;
    for i=3:length(x)
        if (y(i)==y(i-2) && x(i)==x(i-2))
            start=i;break;
        end
    end
    if (start~=1)
        x=x(start-1:start-1+(length(x)+1)/2-1);
        y=y(start-1:start-1+(length(y)+1)/2-1);
    end
    
    %t is the parameter in relation to which the curve is parametrized xD
    N{b}=length(x);
    
    %here we use a gaussian filter to compute the first and second
    %derivatives, x'(t) x''(t) and y'(t) y''(t), of the parametrizations, x(t) and y(t)
    
    %the std dev of the gaussian
    if nargin<3; sigma=25; end
    
    %the length of the filter: radius=2*sigma
    fLen=4*sigma+1;
    
    %the gaussian filter
    Gf = @(x) 1/sigma/sqrt(2*pi)*exp(-1/2*(x/sigma).^2);
    G = Gf(-(fLen-1)/2:(fLen-1)/2);
    factor=sum(G);
    G=G/factor;
    
    %smooth x(t) and y(t)
    xt{b}=conv(x,G,'valid');
    yt{b}=conv(y,G,'valid');
    
    %compute the first derivative of the gaussian filter and scale it so
    %that the result of its application is rightly scaled
    Gf1 = @(x) -x/sigma^2.*Gf(x);
    G1 = Gf1(-(fLen-1)/2:(fLen-1)/2);
    factor=sum(-G1.*((1:fLen)-(fLen+1)/2));
    G1=G1/factor;
    
    %use the above filter to estimate x'(t) and y'(t)
    dxdt=conv(x,G1,'valid');
    dydt=conv(y,G1,'valid');
    
    %compute the second derivative of the gaussian filter, force the 0th
    %term to annihilate the others' contributions to f(x) in the Taylor
    %Series and scale it so that the result of its application is rightly
    %scaled
    Gf2 = @(x) (x.^2-sigma^2)/sigma^4.*Gf(x);
    G2 = Gf2(-(fLen-1)/2:(fLen-1)/2);
    G2((fLen+1)/2)=-sum(G2)+G2((fLen+1)/2);
    factor=sum(G2.*(((1:fLen)-(fLen+1)/2).^2))/2;
    G2=G2/factor;
    
    %use the above filter to estimate x''(t) and y''(t)
    d2xdt2=conv(x,G2,'valid');
    d2ydt2=conv(y,G2,'valid');
    
    %computing dT/dt=dT/ds*ds/dt=k*ds/dt and ds/dt
    dTdt{b}=(dydt.*d2xdt2-d2ydt2.*dxdt)./(dxdt.^2+dydt.^2);
    dsdt{b}=sqrt(dxdt.^2+dydt.^2);
    
    dsdT{b}=dsdt{b}./abs(dTdt{b});
end