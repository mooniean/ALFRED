function maxVrendV = vesselness2D(V, sigmas, spacing, tau, whiteondark)
% calculates the vesselness probability map (local tubularity) of a 3d input 
% image
% 
% maxVrendV = vesselness2D(V, sigmas, spacing, tau, whiteondark)
% 
% inputs,
%   V : 2D image
%   sigmas : vector of scales on which the vesselness is computed
%   spacing : input image spacing
%   tau : (between 0.5 and 1) : parameter that controls sensitivity
%   whiteondark: (true/false) : are vessels bright on dark background or dark
%       on dark
%
%
% outputs,
%   maxVrendV: maximum vesselness response over scales sigmas
%
% example:
%   V = vesselness2D(I, 1:5, [1;1;1], 1, false);
% Function is written by T. Jerman, University of Ljubljana (October 2014)

global DEBUG;

if nargin<5
    whiteondark = false; % default mode for 2D is dark vessels compared to the background
end

for j = 1:length(sigmas)
    
    if DEBUG; disp(['Current Filter Sigma: ' num2str(sigmas(j)./spacing') ]); end
    
    [Lambda1, Lambda2] = volumeEigenvalues(V,sigmas(j),spacing,whiteondark); 
    
    % proposed filter
    Lambda3 = Lambda2;
    Lambda3(Lambda3<0 & Lambda3 >= tau .* min(Lambda3(:)))=tau .* min(Lambda3(:));
    VrendV = (abs(Lambda2).*abs(Lambda2).*abs(Lambda3-Lambda2)).* 27 ./ (2.*abs(Lambda2)+abs(Lambda3-Lambda2)).^3;    
    VrendV(Lambda2<Lambda3./2)=1;    
    VrendV(Lambda2>=0) = 0;     
    VrendV(~isfinite(VrendV)) = 0;

    %keep max response
    if(j==1)
        maxVrendV=VrendV;
    else        
        maxVrendV = max(maxVrendV,VrendV);
    end
        
    clear VrendV    
     
end

maxVrendV = maxVrendV ./ max(maxVrendV(:));
    


function [Lambda1, Lambda2] = volumeEigenvalues(V,sigma,spacing,whiteondark)
% calculates the three eigenvalues for each voxel in a volume

% Calculate 2D hessian
[Hxx, Hyy, Hxy] = Hessian2D(V,spacing,sigma);

% Correct for scaling
c=sigma.^2;
Hxx = c*Hxx; 
Hxy = c*Hxy;
Hyy = c*Hyy;

if whiteondark == false
    c=-1;
    Hxx = c*Hxx; 
    Hxy = c*Hxy;
    Hyy = c*Hyy;   
end

% compute vesselness only where needed
B1 = - (Hxx+Hyy);
B2 = Hxx .* Hyy - Hxy.^2;

T = ones(size(B1));
T(B1<0) = 0;
T(B2==0 & B1 == 0) = 0;

clear B1 B2;

indeces = find(T==1);

Hxx = Hxx(indeces);
Hyy = Hyy(indeces);
Hxy = Hxy(indeces);

% Calculate eigen values
[Lambda1i,Lambda2i]=eigvalOfHessian2D(Hxx,Hxy,Hyy);

% Free memory
clear Hxx Hyy Hxy;

Lambda1 = zeros(size(T));
Lambda2 = zeros(size(T));

Lambda1(indeces) = Lambda1i;
Lambda2(indeces) = Lambda2i;

% some noise removal
Lambda1(~isfinite(Lambda1)) = 0;
Lambda2(~isfinite(Lambda2)) = 0;

Lambda1(abs(Lambda1) < 1e-4) = 0;
Lambda2(abs(Lambda2) < 1e-4) = 0;


function [Dxx, Dyy, Dxy] = Hessian2D(Volume,spacing,Sigma)
%  This function Hessian3D filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
% 
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,Sigma)
% 
% inputs,
%   I : The image volume, class preferable double or single
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.
%
% outputs,
%   Dxx, Dyy, Dxy: The 2nd derivatives

if nargin < 3, Sigma = 1; end

if(Sigma>0)
    %F=imbigaussian(Volume,Sigma,0.5);
    F=imgaussian(Volume,Sigma,spacing);
else
    F=Volume;
end

% Create first and second order diferentiations
Dy=gradient2(F,'y');
Dyy=(gradient2(Dy,'y'));
clear Dy;

Dx=gradient2(F,'x');
Dxx=(gradient2(Dx,'x'));
Dxy=(gradient2(Dx,'y'));
clear Dx;

function D = gradient2(F,option)
% Example:
%
% Fx = gradient2(F,'x');

[k,l] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:) = (F(2,:) - F(1,:));
    D(k,:) = (F(k,:) - F(k-1,:));
    % Take centered differences on interior points
    D(2:k-1,:) = (F(3:k,:)-F(1:k-2,:))/2;
case 'y'
    D(:,1) = (F(:,2) - F(:,1));
    D(:,l) = (F(:,l) - F(:,l-1));
    D(:,2:l-1) = (F(:,3:l)-F(:,1:l-2))/2;
otherwise
    disp('Unknown option')
end
        
function I=imgaussian(I,sigma,spacing,siz)
% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an 
% Gaussian filter. This function uses for filtering IMFILTER or if 
% compiled the fast  mex code imgaussian.c . Instead of using a 
% multidimensional gaussian kernel, it uses the fact that a Gaussian 
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: 2D input image
%   SIGMA: The sigma used for the Gaussian kernel
%   SPACING: input image spacing
%   SIZ: Kernel size (single value) (default: sigma*6)
% 
% outputs,
%   I: The gaussian filtered image
%

if(~exist('siz','var')), siz=sigma*6; end

if(sigma>0)

    % Filter each dimension with the 1D Gaussian kernels\
    x=-ceil(siz/spacing(1)/2):ceil(siz/spacing(1)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(1))^2)));
    H = H/sum(H(:));    
    Hx=reshape(H,[length(H) 1]);
    
    x=-ceil(siz/spacing(2)/2):ceil(siz/spacing(2)/2);
    H = exp(-(x.^2/(2*(sigma/spacing(2))^2)));
    H = H/sum(H(:));    
    Hy=reshape(H,[1 length(H)]);
    
    I=imfilter(imfilter(I,Hx, 'same' ,'replicate'),Hy, 'same' ,'replicate');
end

function [Lambda1,Lambda2]=eigvalOfHessian2D(Dxx,Dxy,Dyy)
% This function calculates the eigen values from the
% hessian matrix, sorted by abs value

% Compute the eigenvectors of J, v1 and v2
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);

% Compute the eigenvalues
mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);

Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);

