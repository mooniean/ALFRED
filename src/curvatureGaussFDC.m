function [dTdt, dsdt, xt, yt, results] = curvatureGaussFDC(coordX,coordY,sigma)

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
    
    %t is the parameter in relation to which the curve is parametrized xD
    %     t=1:length(x);
    
    %here we use (1) a gaussian filter to smooth the curve and (2) an 
    %FDC filter to compute the first and second
    %derivatives, x'(t) and x''(t), of the parametrizations, x(t) and y(t)
    
    %the std dev of the gaussian
    if nargin<3; sigma=25; end
    
    %the length of the filter: radius=2*sigma
    fGLen=4*sigma+1;
    fDLen=3;
    fLen=fGLen+fDLen-1;
    
    %estimate x'(t) and y'(t)
    dxdt{b}=conv(conv(x,gaussFilter(0,(fGLen-1)/4),'valid'),...
        -finDiffCoeffs(1,-(fDLen-1)/2:(fDLen-1)/2),'valid');
    dydt{b}=conv(conv(y,gaussFilter(0,(fGLen-1)/4),'valid'),...
        -finDiffCoeffs(1,-(fDLen-1)/2:(fDLen-1)/2),'valid');
    
    %estimate x''(t) and y''(t)
    d2xdt2=conv(conv(x,gaussFilter(0,(fGLen-1)/4),'valid'),...
        finDiffCoeffs(2,-(fDLen-1)/2:(fDLen-1)/2),'valid');
    d2ydt2=conv(conv(y,gaussFilter(0,(fGLen-1)/4),'valid'),...
        finDiffCoeffs(2,-(fDLen-1)/2:(fDLen-1)/2),'valid');
    
    %computing dT/dt=dT/ds*ds/dt=k*ds/dt and ds/dt
    dTdt{b}=(dydt{b}.*d2xdt2-d2ydt2.*dxdt{b})./(dxdt{b}.^2+dydt{b}.^2);
    dsdt{b}=sqrt(dxdt{b}.^2+dydt{b}.^2);
    
    %store x(t) and y(t)
    xt{b}=x((fLen+1)/2:end-(fLen-1)/2);
    yt{b}=y((fLen+1)/2:end-(fLen-1)/2);
end

%from here on, we have all the data, this is just a fancy way of showing
%the data we have... it's better if we make this interactive
%and to make it interactive you HAVE to output the positions as well
results = []; 
for b=1:length(boundaries)
    %compute all the available radii for the current boundary
    radii=dsdt{b}./abs(dTdt{b});
    results = [results radii]; %#ok<*AGROW>
end

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

end

function [filter,accuracy]=gaussFilter(m,s)

% desired derivative order
if nargin<1; m=1; end

% sigma
if nargin<2; s=2; end

%the length of the filter: radius=2*sigma
fLen=4*s+1;

%the gaussian filter
Gf=fspecial('gauss',[1 fLen], s);

%compute the first derivative of the gaussian filter and scale it so
%that the result of its application is rightly scaled
G1=-gradient(Gf);
factor=sum(G1.*((1:fLen)-(2*s+1)));
G1=G1/factor;

%compute the second derivative of the gaussian filter, force the 0th
%term to annihilate the others' contributions to f(x) in the Taylor
%Series and scale it so that the result of its application is rightly
%scaled
G2 = gradient(gradient(Gf));
G2(2*s+1)=-sum(G2)+G2(2*s+1);
factor=sum(G2.*(((1:fLen)-(2*s+1)).^2))/2;
G2=G2/factor;

if     m==0; filter=Gf;
elseif m==1; filter=G1;
elseif m==2; filter=G2; end

% compute order of accuracy h^a - I just compute a, h is the grid spacing
if     m==0; accuracy=Inf;
elseif m==1; accuracy=2;
elseif m==2; accuracy=2; end

end