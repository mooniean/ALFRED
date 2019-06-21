function [dTdt, dsdt, xt, yt, N, dsdT] = curvatureSpline(coordX,coordY,tol)

%build a binary image out of the coordinates
offset=20; %zero-pad image so I have like a frame around it
skeleton = zeros(max(coordX)+offset,max(coordY)+offset);
for i = 1:length(coordX)
    skeleton(coordX(i)+offset/2,coordY(i)+offset/2) = 1;
end

%find branchpoints
brPoints=bwmorph(skeleton,'branchpoints');

%eliminate branchpoints and their 24 nearest neighbours from the skeleton
[rowBr,colBr] = find(brPoints);
for k=1:length(rowBr)
    for i=-2:2
        for j=-2:2
            skeleton(rowBr(k)+i,colBr(k)+j)=0;
        end
    end
end

%clean groups of up to 10 pixels that are isolated in the skeleton
skeleton=bwareaopen(skeleton,10);

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
    L=length(x); N{b}=L; t=1:L;
    
    %here we use a smoothing spline to compute the first and second
    %derivatives, x'(t) x''(t) and y'(t) y''(t), of the parametrizations, x(t) and y(t)
    
    ppX = spaps( t, x, tol*L);
    ppY = spaps( t, y, tol*L);
    
    xt{b}=@(t) fnval(ppX,t);
    yt{b}=@(t) fnval(ppY,t);
    
    dxdt=@(t) fnval(fnder(ppX),t);
    dydt=@(t) fnval(fnder(ppY),t);
    
    d2xdt2=@(t) fnval(fnder(fnder(ppX)),t);
    d2ydt2=@(t) fnval(fnder(fnder(ppY)),t);
    
    dTdt{b}=@(t) (dydt(t).*d2xdt2(t)-d2ydt2(t).*dxdt(t))./(dxdt(t).^2+dydt(t).^2);
    dsdt{b}=@(t) sqrt(dxdt(t).^2+dydt(t).^2);
    
    dsdT{b}=@(t) dsdt{b}(t)./abs(dTdt{b}(t));
end