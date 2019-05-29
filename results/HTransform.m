% Hough Transform

clc; clear; close all

n=1e2;
canvas=zeros(n);

m1=1;
b1=0;

m2=.3;
b2=65;

for x=1:n
    canvas(round(m1*x+b1),x)=1;
    canvas(round(m2*x+b2),x)=1;
end

imshow(canvas)
axis square

[I,J]=find(canvas);
lines=[];
rho=[];
theta=-pi/2:1e-2:pi/2-1e-2;

for i=1:length(I)
    rho(end+1,:)=round(J(i)*cos(theta)+I(i)*sin(theta));
    
    for j=i+1:length(I)
        thetaS=atan((J(j)-J(i))/(I(i)-I(j)));
        rhoS=I(j)*sin(thetaS)+J(j)*cos(thetaS);
        lines(end+1,:)=[thetaS,rhoS];
    end
end

d=round(sqrt(2*n^2));

matrix=zeros(2*d+1,round(pi/0.01)+1);

for i=1:size(lines,1)
    matrix(round(lines(i,2)+d)+1,round(round(lines(i,1)+pi/2,2)/0.01)+1)=matrix(round(lines(i,2)+d)+1,round(round(lines(i,1)+pi/2,2)/0.01)+1)+1;
end

figure
imshow(imadjust(rescale(matrix)))

%%%%%%%%

d=round(sqrt(2*n^2));
maxRho=max(max(rho));

matrix=zeros(2*d+1,length(theta));

for i=1:size(rho,1)
    for j=1:length(theta)
        matrix(round(rho(i,j)+d)+1,j)=matrix(round(rho(i,j)+d)+1,j)+1;
    end
end

figure
imshow(imadjust(rescale(matrix)))

figure
imshow(imadjust(rescale(hough(canvas,'RhoResolution',1,'Theta',-90:0.573:90-0.573))))

