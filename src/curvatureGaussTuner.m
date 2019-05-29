function [dTdt, dsdt, results] = curvatureGaussTuner(coordX,coordY)

sigma=.5;
results=Inf;
resultsAll={};

while(~isempty(results))
    resultsAll{end+1}=results;
    [~,~,~,~,results]=curvatureGauss(coordX,coordY,sigma);
    sigma=sigma+.5;
end

resultsAll=resultsAll(2:end);

f = figure('IntegerHandle','off',...
    'Name','Curvature - Gaussian filter tuning: tick to finish');
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);

b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
    'value',0.5, 'min',0.5, 'max',length(resultsAll)*.5);
bgcolor = f.Color;

bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
    'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
    'String',num2str(length(resultsAll)*.5),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
    'String','sigma','BackgroundColor',bgcolor);

b.Callback = {@sliderCallback,ax,resultsAll};

c=uicontrol('Parent',f,'Style','checkbox','Position',[520,54,23,23],...
    'String','Confirm');

waitfor(c,'Value');

[dTdt, dsdt,~,~,results]=curvatureGauss(coordX,coordY,round(b.Value/.5)*.5);

end

function sliderCallback(es,~,ax,resultsAll)
histogram(ax,1./resultsAll{round(es.Value/.5)},200,'facecolor',[254 224 144]./255,...
    'facealpha',.7,'edgecolor','none','normalization','probability');
xlim([0,0.1]);
end