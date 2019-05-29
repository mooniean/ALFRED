% explicar matematicamente como calculas a curvatura
% explicar quais os métodos do MATLAB que usas, que devíamos fazer o fit entre branch points e não tudo à grande (mesmo assim, a curvatura dá bastante bem...)
% explicar que ainda tens problemas a contar distâncias repetidas... isto é boundaries que na realidade repetem outras...
% para a curvatura isso não deve ser problema, mas para a distância é...
% explicar como vais pesar a curvatura em cada boundary/branch:
% passas weightedK=sum(K em cada ponto * deltaS do ponto) e a S=sum(delta S do ponto)
% e depois, no fim, fazes a média para todas as boundaries/branches: sum(wightedK)/sum(S)
% podes dar exemplos com figuras (não é preciso esconder o que nós sabemos que tem que ser melhorado)
% podes dizer que esta nova versão é muito mais rápida e evita parametrizações complicadas tipo o tamanho da janela para o ajuste
% basicamente explicar muito bem o que há para fazer em vez de focar em fazer

function [dTdt, dsdt, radiiRot] = curvatureFourier(coordX,coordY,rotations)

%build a binary image out of the coordinates
offset=20; %zero-pad image so I have like a frame around it
skeleton = zeros(max(coordX)+offset,max(coordY)+offset);
for i = 1:length(coordX)
    skeleton(coordX(i)+offset/2,coordY(i)+offset/2) = 1;
end

if nargin<3; rotations=0:5:85; end

noTerms=1;
radiiRot=zeros([noTerms,length(rotations)]);

for ind=1:length(rotations)
    
    rot=rotations(ind);
    skeleton = imrotate(skeleton,rot);
    
    %skeleton = bwmorph(skeleton,'skel');
    
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
    % b=bwboundaries(skeleton,'noholes');
    b=bwboundaries(skeleton);
    
    for terms=1:noTerms
        % prepare Fourier fit
        ft = fittype( ['fourier',int2str(terms)] );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares');
        
        dTdt=[];dsdt=[];
        
        gofT=zeros(1,length(b));
        
        for h=1:length(b)
            x=b{h}(:,1).';y=b{h}(:,2).';
            
            start=1;
            for i=3:length(x) % check periodicity
                if (y(i)==y(i-2) && x(i)==x(i-2))
                    start=i;break;
                end
            end
            if (start~=1)
                x=x(start-1:start-1+(length(x)+1)/2-1);
                y=y(start-1:start-1+(length(y)+1)/2-1);
            end
            t=1:length(x);
            
            [fitresultX, gofX] = fit( t.', x.', ft, opts );
            [fitresultY, gofY] = fit( t.', y.', ft, opts );
            gofT(h)=gofX.rsquare;gofT(h)=gofY.rsquare;
            
            coeffvaluesX=coeffvalues(fitresultX);
            coeffvaluesY=coeffvalues(fitresultY);
            
            resX=coeffvaluesX(1)*ones(size(t));
            resY=coeffvaluesY(1)*ones(size(t));
            
            for i=2:2:length(coeffvaluesX)-1
                resX=resX+coeffvaluesX(i)*cos(i/2*coeffvaluesX(end)*t);
                resY=resY+coeffvaluesY(i)*cos(i/2*coeffvaluesY(end)*t);
                resX=resX+coeffvaluesX(i+1)*sin(i/2*coeffvaluesX(end)*t);
                resY=resY+coeffvaluesY(i+1)*sin(i/2*coeffvaluesY(end)*t);
            end
            
            firstDerX=zeros(size(t));
            firstDerY=zeros(size(t));
            
            for i=2:2:length(coeffvaluesX)-1
                firstDerX=firstDerX-i/2*coeffvaluesX(end)*coeffvaluesX(i)*sin(i/2*coeffvaluesX(end)*t);
                firstDerY=firstDerY-i/2*coeffvaluesY(end)*coeffvaluesY(i)*sin(i/2*coeffvaluesY(end)*t);
                firstDerX=firstDerX+i/2*coeffvaluesX(end)*coeffvaluesX(i+1)*cos(i/2*coeffvaluesX(end)*t);
                firstDerY=firstDerY+i/2*coeffvaluesY(end)*coeffvaluesY(i+1)*cos(i/2*coeffvaluesY(end)*t);
            end
            
            secondDerX=zeros(size(t));
            secondDerY=zeros(size(t));
            
            for i=2:2:length(coeffvaluesX)-1
                secondDerX=secondDerX-(i/2*coeffvaluesX(end))^2*coeffvaluesX(i)*cos(i/2*coeffvaluesX(end)*t);
                secondDerY=secondDerY-(i/2*coeffvaluesY(end))^2*coeffvaluesY(i)*cos(i/2*coeffvaluesY(end)*t);
                secondDerX=secondDerX-(i/2*coeffvaluesX(end))^2*coeffvaluesX(i+1)*sin(i/2*coeffvaluesX(end)*t);
                secondDerY=secondDerY-(i/2*coeffvaluesY(end))^2*coeffvaluesY(i+1)*sin(i/2*coeffvaluesY(end)*t);
            end
            
            dTdtCurBound=abs(firstDerX.*secondDerY-firstDerY.*secondDerX)./(firstDerX.^2+firstDerY.^2).^(2/2);
            dsdtCurBound=sqrt(firstDerX.^2+firstDerY.^2);
            
            dTdt=[dTdt,dTdtCurBound]; %#ok<*AGROW>
            dsdt=[dsdt,dsdtCurBound];
            
        end
        
        % changed by Nuno (28/04/19)
        % compute mean radius weighted by arc length
        
        radiiRot(terms,ind) = trapz(dsdt./dTdt.*dsdt)./trapz(dsdt);
        
    end
end

radiiRot=radiiRot(~isnan(radiiRot));

end