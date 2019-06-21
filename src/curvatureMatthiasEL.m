clc; clear; close all;
figure('Units','normalized','OuterPosition',[0 0 1 1]);
colormap jet

% ellipse's axes parameters

ax=400;
by=200;

for alphad=0%:10:90
    
    alpha=alphad *pi/180;
    
    % analytic expressions for the tilted ellipse
    
    a = @(x) (sin(alpha)/ax)^2+(cos(alpha)/by)^2;
    b = @(x) 2*x*cos(alpha)*sin(alpha)*(1/ax^2-1/by^2);
    c = @(x) ((cos(alpha)/ax)^2+(sin(alpha)/by)^2)*x.^2-1;
    
    f = @(x) (-b(x)+sqrt(b(x).^2-4*a(x)*c(x)))/2/a(x);
    g = @(x) (-b(x)-sqrt(b(x).^2-4*a(x)*c(x)))/2/a(x);
    
    % "analytic" mean radius
    limX=sqrt((ax*cos(alpha))^2+(by*sin(alpha))^2);
    
    f1=matlabFunction(diff(sym(f),1));
    g1=matlabFunction(diff(sym(g),1));
    f2=matlabFunction(diff(sym(f),2));
    g2=matlabFunction(diff(sym(g),2));
    
    kdldx = @(x) abs(f2(x)./(1+f1(x).^2));
    dldx = @(x) (1+f1(x).^2).^(1/2);
    k = @(x) kdldx(x)./dldx(x);
    dldxk = @(x) dldx(x)./k(x);
    
    analyticalRadius=integral(dldxk,-limX,limX)/integral(dldx,-limX,limX);
    radiusDensity=matlabFunction(sym(dldx)/diff(sym(kdldx)/sym(dldx),1)*...
        (sym(kdldx)/sym(dldx))^2/integral(dldx,-limX,limX));
    
    % digitalization
    
    xD=-limX:1e-4:limX; % D = Digitized, DC = Digitized Complete range
    yDC=[round(f(xD)),round(g(xD))];
    
    xDC=round([xD,xD]);
    xDC(imag(yDC)~=0)=[];
    yDC(imag(yDC)~=0)=[];
    
    xDC=xDC-min(xDC)+1;
    yDC=yDC-min(yDC)+1;
    
    skeleton = zeros(max(xDC),max(yDC));
    for i = 1:length(xDC)
        skeleton(xDC(i),yDC(i)) = 1;
    end
    
    skeleton=bwmorph(skeleton,'skel','Inf');
    [x,y]=find(skeleton);
    
    % give this to Fourier fit function
    
    [~, ~, fourierRadius] = curvatureFourier(x,y,0);
    
    % give this to Gaussian convolution function
    sigma=25;
    [dTdt, dsdt, xt, yt, ~] = curvatureGauss(x,y,sigma);
    
    %small hack just for this case (periodic conditions)
    copyPoint=floor((length(xt{1})-4*sigma)/2)+1;
    copyRange=copyPoint:copyPoint+2*sigma-1;
    boundaries=bwboundaries(skeleton);
    
    xt{1}=boundaries{1}(:,1).';yt{1}=boundaries{1}(:,2).';
    dTdt{1}=[dTdt{1}(copyRange+2*sigma),dTdt{1},dTdt{1}(copyRange)];
    dsdt{1}=[dsdt{1}(copyRange+2*sigma),dsdt{1},dsdt{1}(copyRange)];
    
    % colorize curvature on the ellipse per constant-length-sector
    sections=45;
    
    ellipseCircumference=trapz(dsdt{1});
    secLen=ellipseCircumference/sections;
    
    deltaS=cumtrapz(dsdt{1});
    breakpoints=ones(1,sections+1);
    for i=1:sections-1; [~,breakpoints(i+1)]=min(abs(deltaS-secLen*i)); end
    breakpoints(end)=length(dsdt{1});
    
    meanGaussRadius=trapz(abs(dsdt{1}./dTdt{1}.*dsdt{1}))/trapz(dsdt{1});
    meanGaussRadii=zeros(1,sections);
    offset=20; coloredSkeleton=zeros(max(xt{1})+offset,max(yt{1})+offset);
    for i=1:sections
        meanGaussRadii(i)=trapz(abs(dsdt{1}(breakpoints(i):breakpoints(i+1))...
            ./dTdt{1}(breakpoints(i):breakpoints(i+1))...
            .*dsdt{1}(breakpoints(i):breakpoints(i+1))))...
            /trapz(dsdt{1}(breakpoints(i):breakpoints(i+1)));
        for j=breakpoints(i):breakpoints(i+1)-1
            colorX=xt{1}(j)+offset/2;colorY=yt{1}(j)+offset/2;
            coloredSkeleton(colorX-3:colorX+3,colorY-3:colorY+3)=log10(abs(meanGaussRadii(i)));
        end
    end
    
    lim=3.5;
    coloredSkeleton(1,1)=lim; assert(max(max(coloredSkeleton))==lim);
    
    % plotting
    
    subplot(121)
    imagesc(coloredSkeleton); daspect([1 1 1]); h=colorbar;
    h.TickLabels=num2str(10.^(h.Ticks.'),'%.2f');
    title(['a = ',int2str(ax),', b = ',int2str(by),...
        ', \alpha = ',int2str(alphad),...
        ' \circ. Gaussian Filter: \sigma = ',num2str(sigma),...
        '. Segments: ',int2str(sections)])
    
    subplot(122)
    binWidth=50;
    histogram(meanGaussRadii,0:binWidth:1000,'Normalization','probability');
    yLim=0.35;ylim([0 yLim]);xlabel('Radius (px)'),ylabel('Relative frequency');
    hold on; line([ax^2/by ax^2/by],[0 yLim],'Color','r'); line([by^2/ax by^2/ax],[0 yLim],'Color','r');
    line([meanGaussRadius meanGaussRadius],[0 yLim],'Color','g');
    line([fourierRadius fourierRadius],[0 yLim],'Color','b');
    line([analyticalRadius analyticalRadius],[0 yLim],'Color','k');
    plot(dldx(xD(2:end-1))./kdldx(xD(2:end-1)),2*abs(radiusDensity(xD(2:end-1)))*binWidth); hold off;
    legend({'Gaussian section radius','Minimum radius','Maximum radius',...
        'Gaussian radius (average)','Fourier radius',...
        '"Analytic" (average)','Probability density (\times binWidth)'})
    
    %if ~alphad; gif(['ellipse_',int2str(sections),'.gif'],'frame',gcf,'DelayTime',2); else; gif; end
end