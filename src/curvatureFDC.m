function [dTdt, dsdt, results] = curvatureFDC(coordX,coordY,quarterWinSize)

skeleton = zeros(max(coordX)+10,max(coordY)+10);
for i = 1:length(coordX)
    skeleton(coordX(i),coordY(i)) = 1;
end

skeleton=bwmorph(skeleton,'skel');
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

%storage for x(t) and y(t), where the curvature computation is valid
xt=cell(length(boundaries),1);
yt=cell(length(boundaries),1);

%storage for x'(t) and y'(t), where the curvature computation is valid
dxdt=cell(length(boundaries),1);
dydt=cell(length(boundaries),1);

%storage for derivatives of T (the curve's tangent vector) and s (the curve
%'s arc length) for each boundary
dTdt=cell(length(boundaries),1);
dsdt=cell(length(boundaries),1);

for b=1:length(boundaries)
    
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
    
    %here we use an FDC filter to compute the first and second
    %derivatives, x'(t), y'(t) and x''(t), y''(t), of the parametrizations,
    %x(t) and y(t)
    
    %given we cleaned groups of up to 10 isolated pixels, we'll use a
    %9-point filter to be safe
    
    if nargin<3; quarterWinSize=2; end
    n=4*quarterWinSize+1;
    p=-(n-1)/2:(n-1)/2;
    
    G1=fliplr(finDiffCoeffs(1,p));
    G2=fliplr(finDiffCoeffs(2,p));
    
    %estimate x'(t) and y'(t) where we can apply a central finite
    %difference filter
    dxdt{b}=conv(x,G1,'valid');
    dydt{b}=conv(y,G1,'valid');
    
    %estimate x''(t) and y''(t) where we can apply a central finite
    %difference filter
    d2xdt2=conv(x,G2,'valid');
    d2ydt2=conv(y,G2,'valid');
    
    %estimate x'(t), y'(t), x''(t) and y''(t) where we can't apply a
    %central finite difference filter and have to use sided-approximations
    dxdt{b}=[zeros(1,(n-1)/2),dxdt{b},zeros(1,(n-1)/2)];
    dydt{b}=[zeros(1,(n-1)/2),dydt{b},zeros(1,(n-1)/2)];
    
    d2xdt2=[zeros(1,(n-1)/2),d2xdt2,zeros(1,(n-1)/2)];
    d2ydt2=[zeros(1,(n-1)/2),d2ydt2,zeros(1,(n-1)/2)];
    
    for i=1:(n-1)/2
        dxdt{b}(i)=sum(finDiffCoeffs(1,1-i:n-i).*x(1:n));
        
        dydt{b}(i)=sum(finDiffCoeffs(1,1-i:n-i).*y(1:n));
        
        dxdt{b}(end-(n-1)/2+i)=...
            sum(finDiffCoeffs(1,(1-i:n-i)-(n+1)/2).*x(end-n+1:end));
        
        dydt{b}(end-(n-1)/2+i)=...
            sum(finDiffCoeffs(1,(1-i:n-i)-(n+1)/2).*y(end-n+1:end));
        
        d2xdt2(i)=sum(finDiffCoeffs(2,1-i:n-i).*x(1:n));
        
        d2ydt2(i)=sum(finDiffCoeffs(2,1-i:n-i).*y(1:n));
        
        d2xdt2(end-(n-1)/2+i)=...
            sum(finDiffCoeffs(2,(1-i:n-i)-(n+1)/2).*x(end-n+1:end));
        
        d2ydt2(end-(n-1)/2+i)=...
            sum(finDiffCoeffs(2,(1-i:n-i)-(n+1)/2).*y(end-n+1:end));
    end
    
    %computing dT/dt=dT/ds*ds/dt=k*ds/dt and ds/dt
    dTdt{b}=(dydt{b}.*d2xdt2-d2ydt2.*dxdt{b})./(dxdt{b}.^2+dydt{b}.^2);
    dsdt{b}=sqrt(dxdt{b}.^2+dydt{b}.^2);
    
    %store x(t) and y(t)
    xt{b}=x;
    yt{b}=y;
end

results = [];
for b=1:length(boundaries)
    %compute all the available radii for the current boundary
    radii=dsdt{b}./abs(dTdt{b});
    results = [results radii]; %#ok<*AGROW>
end

function [filter,accuracy]=finDiffCoeffs(m,p)

% desired derivative order
if nargin<1; m=1; end

% points
if nargin<2; p=-4:4; end

% no. of equations
nEq=length(p);

% compute Taylor series for neighbours up to the nEq-th term
taylorCoeffs = @(neighbour,nEq) neighbour.^(0:nEq)./factorial(0:nEq);
allTaylorCoeffs=taylorCoeffs(p.',nEq);

% build matrix A of Ax=b, constrain this
% to a number of equations/rows equal to the no. of vars/columns
Ac=allTaylorCoeffs.';
A=Ac(1:nEq,:);

% build vector b of Ax=b
b=zeros(nEq,1);
b(m+1)=1;

% solve the system and, thus, find the filter
filter=(A\b).';

% compute order of accuracy h^a - I just compute a, h is the grid spacing
if Ac(end,:)*filter.'<1e-12; accuracy=nEq-m+1;
else; accuracy=nEq-m; end