close all; clear; clc;

r=50:25:250;
theta=44.427;

alpha=0:3:90;

for aa=alpha
    M=[cosd(aa) -sind(aa); sind(aa) cosd(aa)];
    
    for rr=r
        xP=-rr*sind(theta):rr*sind(theta);
        xR=round(xP);
        
        yP=-sqrt(rr^2-xP.^2);
        yR=round(yP);
        
        n=length(xP);
        pP=M*[xP;yP];
        
        for i=1:2
            
            if i==1; pR=round(pP);
            elseif i==2; pR=round(M*[xR;yR]);
            end
            
            subplot(144),hold on
            if i==1; plot(aa,rms(pP(1,:)-pR(1,:)),'x','Color',[(rr-min(r))/(max(r)-min(r)),.5,.5])
            elseif i==2; plot(aa,rms(pP(1,:)-pR(1,:)),'o','Color',[(rr-min(r))/(max(r)-min(r)),.5,.5])
            end
            if i==1; plot(aa,rms(pP(2,:)-pR(2,:)),'x','Color',[(rr-min(r))/(max(r)-min(r)),.5,.5])
            elseif i==2; plot(aa,rms(pP(2,:)-pR(2,:)),'o','Color',[(rr-min(r))/(max(r)-min(r)),.5,.5])
            end
            line([0, 90],sqrt(1/12*ones(1,2)),'LineStyle','--','Color','k')
            xlim([0,90]),xlabel('Rotation angle'),ylabel('RMSE')
            hold off
            
            t=1:n;tol=10^-.7;
            ppX = spaps( t, pR(1,:), tol*n);
            ppY = spaps( t, pR(2,:), tol*n);
            
            xtR=@(t) fnval(ppX,t);
            ytR=@(t) fnval(ppY,t);
            
            dxdtR=@(t) fnval(fnder(ppX),t);
            dydtR=@(t) fnval(fnder(ppY),t);
            
            d2xdt2R=@(t) fnval(fnder(ppX,2),t);
            d2ydt2R=@(t) fnval(fnder(ppY,2),t);
            
            kdldtR=@(t) abs(dydtR(t).*d2xdt2R(t)-d2ydt2R(t).*dxdtR(t))./(dxdtR(t).^2+dydtR(t).^2);
            dldtR=@(t) sqrt(dxdtR(t).^2+dydtR(t).^2);
            kR = @(t) kdldtR(t)./dldtR(t);
            dldtkR = @(t) dldtR(t)./kR(t);
            
            subplot(1,4,i+1),hold on
            plot(asind(xP/rr),1./kR(1:n),'Color',[(rr-min(r))/(max(r)-min(r)),.5,.5])
            ylim([0,300]),xlabel('\theta'),ylabel('Radius of curvature')
            line([-theta,theta],[rr,rr],'Color','k'), xlim([-theta,theta])
            line([0, 0],[0, 300],'LineStyle','--','Color','k')
            hold off
            
        end
        
        subplot(141),hold on;
        fill(pR(1,:)+[-.5;.5;.5;-.5],pR(2,:)+[.5;.5;-.5;-.5],.5*ones(1,3))
        plot(pP(1,:),pP(2,:),'.','Color',[(rr-min(r))/(max(r)-min(r)),.5,.5])
        xlim([-(n+1)/2,max(r)]+.5),ylim([-max(r)-.5,(n+1)/2]),axis square
        hold off
    end
    
    subplot(141)
    line([0, pP(1,1)],[0, pP(2,1)],'Color','k')
    line([0, rr*sind(aa)],[0, -rr*cosd(aa)],'LineStyle','--','Color','k')
    line([0, pP(1,end)],[0, pP(2,end)],'Color','k')
    
    if (aa==alpha(end)); break; end
    %pause()
    subplot(141),cla
end

subplot(142),title('Rounded solely after rotating')
%subplot(143),title('Rounded solely before rotating')
subplot(143),title('Rounded before and after rotating')