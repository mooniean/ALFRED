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