function varargout = ALFRED_mainpanel(varargin)
% ALFRED_MAINPANEL MATLAB code for ALFRED_mainpanel.fig
%      ALFRED_MAINPANEL, by itself, creates a new ALFRED_MAINPANEL or raises the existing
%      singleton*.
%      Main panel for individual image analysis.
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ALFRED_mainpanel

% Last Modified by GUIDE v2.5 11-Mar-2019 08:32:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ALFRED_mainpanel_OpeningFcn, ...
    'gui_OutputFcn',  @ALFRED_mainpanel_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ALFRED_mainpanel is made visible.
function ALFRED_mainpanel_OpeningFcn(hObject, ~, handles, varargin)
set ( 0, 'DefaultFigureColor', [0.94 0.94 0.94] )
handles.output = hObject;

handles.imageGroup = varargin{1};           % Group of images
handles.totalNumberImages = varargin{2};    % Total image number in the group
handles.currentImage = varargin{3};         % Current image number
handles.fileName = varargin{4};             % Name of the file(s)
handles.metadata = varargin{5};             % Metadata of the image(s)
% numvalue = metadata{1}.getPixelsPhysicalSizeX(0).value();
%     unit = metadata{1}.getPixelsPhysicalSizeX(0).unit().getSymbol();
handles.imagePanel = false;
if ~iscell(handles.fileName)
    temp = handles.fileName;
    handles.fileName = {temp};
    handles.imagePanel = true;
end

set(handles.figure1, 'Name', ['ALFRED Image Number ' num2str(handles.currentImage) ' from ' char(handles.fileName{handles.currentImage})]);
handles.display=axes('Parent',handles.displaypanel);

% Setting icon of the program
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1,'javaframe');
jIcon=javax.swing.ImageIcon('icons\alfred_50x50.png');
jframe.setFigureIcon(jIcon);

% Checkpoint - only one image was passed from the first GUI
if ~((handles.totalNumberImages ~= 1) || (handles.totalNumberImages ~= length(handles.imageGroup)))
    handles.image = handles.imageGroup;
    set(handles.continuebutton,'Enable','off');
    set(handles.previousButton,'Enable','off');
end
handles.image = handles.imageGroup{handles.currentImage};
if handles.currentImage == 1
    set(handles.previousButton,'Enable','off');
elseif handles.currentImage == handles.totalNumberImages || handles.currentImage == length(handles.imageGroup)
    set(handles.continuebutton,'Enable','off');
    set(handles.previousButton,'Enable','on');
end

handles.channelNumber = size(handles.image,3);
channelOpt = ['All' newline num2str(linspace(1,handles.channelNumber, handles.channelNumber),'%d\n')];
set(handles.channelselectmenu,'String',channelOpt);

handles.subplotdisplay = {};
handles.channelChoice = 0;
handles.regionChoice = 'Rectangle';
handles.virtualScaleValue = 3;

handles.finalROIs = struct([]);
handles.imageIndex = 1;
handles.ROIindex = 1;
handles.baseImageIndex = '0';
handles.conversionFactor = 1;

handles.finalROIs(handles.imageIndex).imageName = handles.fileName{handles.currentImage};
handles.finalROIs(handles.imageIndex).ROInumber =  handles.baseImageIndex;
handles.finalROIs(handles.imageIndex).type = 'Full';
handles.finalROIs(handles.imageIndex).conversionFactor = handles.conversionFactor;
% handles.finalROIs(handles.imageIndex).metadata = handles.metadata{handles.currentImage};


% Variable Initialisation - Possible change into a different file.

set(handles.imagenormalisecheck,'Value',1);

set(handles.textdisplay,'String',['Image Number ', num2str(handles.currentImage) ,' ', num2str(size(handles.image,1)), 'x', num2str(size(handles.image,2)) ]);
set(handles.imageSpecDisplay,'String','1 um');
% set(handles.savebutton,'Enable','off');
set(handles.calculateallbutton,'Enable','off');
% set(handles.virtualisationbutton,'Enable','off');
set(handles.regionbutton,'Enable','off');
set(handles.freehand,'Enable','off');
set(handles.rectangle,'Enable','off');
set(handles.virtualisationbutton,'Enable','off');

% set(handles.restorebutton,'Enable','off');
% set(handles.applygrayscalecontrastchange,'Enable','off');
set(handles.applygrayscalethreshold,'Enable','off');

checkuint(hObject, handles)
handles.selectedRegion = handles.image;
handles.region = handles.image;
cla
if size(handles.image,3)>3
    imshow(handles.image(:,:,1));
    text(0,100,'Channel 1','Color','white');
    handles.oldImage = handles.image(:,:,1);
else
    imshow(handles.image)
    handles.oldImage = handles.image;
end
% Choose default command line output for ALFRED_mainpanel
handles.output = hObject;
set(handles.undobutton,'Enable','off');
set(handles.regionbutton,'Enable','off');


% Update handles structure
guidata(hObject, handles);
% UIWAIT makes ALFRED_mainpanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ALFRED_mainpanel_OutputFcn(~, ~, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ------------------------------------------------------------------------------------------------------------------------------------------------
% Button Press Functions

% --- Executes on button press in virtualisationbutton.
function virtualisationbutton_Callback(hObject, ~, handles)
handles = virtualising(handles);

guidata(hObject, handles);


function applygrayscalecontrastchange_Callback(hObject, ~, handles)
TOL = [handles.GSContrastMinValue handles.GSContrastMaxValue];
handles.adjusted = imadjust(handles.region,stretchlim(handles.region,TOL),[]);
imshow(handles.adjusted)
checkPrint(handles.figure1,handles)
handles.oldImage = handles.adjusted;
guidata(hObject, handles);
function applygrayscalethreshold_Callback(hObject, ~, handles)
%Verifies that the activatemicrotubule button was pressed.
global DEBUG;
if DEBUG; disp('Applying threshold mask'); end
cla(handles.display,'reset');
if isnan(handles.virtualisedImage)
    handles.maskedImage = handles.region>get(handles.grayscalethresholdslider,'Value');
else
    handles.maskedImage = handles.virtualisedImage>(get(handles.grayscalethresholdslider,'Value')/255);
end
handles.grayscalevalue = get(handles.grayscalethresholdslider,'Value');
imshow(handles.maskedImage);
name = handles.fileName{handles.currentImage};
[existInStruct,index] = findValue(handles.finalROIs,name);
if existInStruct > 1
    for i = 1:numel(index)
        roi = str2double(handles.finalROIs(index(i)).ROInumber);
        if roi == 0
            continue
        end
        x1 = handles.finalROIs(index(i)).box(1);
        x2 = handles.finalROIs(index(i)).box(1)+handles.finalROIs(index(i)).box(3);
        y1 = handles.finalROIs(index(i)).box(2);
        y2 = handles.finalROIs(index(i)).box(2)+handles.finalROIs(index(i)).box(4);
        hold on;
        x = [x1, x2, x2, x1, x1];
        y = [y1, y1, y2, y2, y1];
        switch handles.finalROIs(index(i)).type
            case 'Axon'
                handles.plotroi=plot(x, y, 'm--','markers',12);
            case 'MT disorganisation'
                handles.plotroi=plot(x, y, 'g--','markers',12);
        end
        hold off;
    end
    handles.ROIindex = max(index)+1;
end
checkPrint(handles.figure1,handles)
handles.oldImage = handles.maskedImage;
set(handles.freehand,'Enable','on');
set(handles.rectangle,'Enable','on');
set(handles.rectangle,'Value',0);
set(handles.freehand,'Value',0);
guidata(hObject, handles);

function regionbutton_Callback(hObject, ~, handles)
set(hObject, 'Enable', 'off');
handles.imageIndex = handles.imageIndex + 1;
if strcmp(handles.ROICutChoice,'Rectangle')
%     handles.regiontemp = imrect(handles.display);
    handles.regiontemp = drawrectangle(handles.display);
else
%     handles.regiontemp = imfreehand(handles.display);
    handles.regiontemp = drawfreehand(handles.display);
end

handles.roiFlag = 0;
handles.maskFlag = 1;
handles.mask = createMask(handles.regiontemp);


S = regionprops(handles.mask,'Boundingbox');
handles.b = floor(S(1).BoundingBox);

%CHANGING THIS CODE
% if size(size(handles.region)) == 3 this was code for the colours.
% %     temp = handles.skel;
% %     size(size(handles.region))
%     temp =  handles.maskedImage;
%     imgcut=distributeImage(hObject,temp,handles);
%     handles.region = imgcut;
% else
%     temp = im2double(handles.image);

temp = im2double(handles.image(:,:,handles.channelChoice));

%     mask = repmat(double(handles.mask),[1 1 3]);
temp = temp.*handles.mask;
%     temp = handles.skel;
%      imgcut=distributeImage(hObject,temp,handles);
handles.selectedRegion = temp(handles.b(2):handles.b(2)+handles.b(4), handles.b(1):handles.b(1)+handles.b(3));
% end
temp = handles.maskedImage;
imgcut = temp(handles.b(2):handles.b(2)+handles.b(4), handles.b(1):handles.b(1)+handles.b(3),1);
handles.region = imgcut;
% [skr,~] = skeleton(handles.region);
% Changed the skeleton function being used.
handles.skel = bwmorph(handles.region,'skel',inf);

% assignin('base', 'skel', handles.selectedRegion);
delete(handles.regiontemp)
x1 = handles.b(1);
x2 = handles.b(1)+handles.b(3);
y1 = handles.b(2);
y2 = handles.b(2)+handles.b(4);
hold on;
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];

[handles.ROIregionType, handles.ROIpath, handles.ROIlength,handles.ROIlines,handles.longLine,handles.bigx,handles.bigy,linesDist,linesCoord] = ALFRED_lengthpopup(handles.skel,handles.selectedRegion);
% CHECK WHICH ARE THE NAMES THAT THE TYPES WILL HAVE SO THAT WE CAN HAVE IT
% COLOR CODED LIKE: R IS AXON, B IS DISORGANISATION, C IS NEURON.
if ~strcmp(handles.ROIregionType,'None')
    switch handles.ROIregionType
        case 'Axon'
            handles.plotroi=plot(x, y, 'm--','markers',12);
        case 'MT disorganisation'
            handles.plotroi=plot(x, y, 'g--','markers',12);
            
    end
    
    
    % Save into the structure created:
    
    handles.finalROIs(handles.imageIndex).imageName = handles.fileName{handles.currentImage};
    handles.finalROIs(handles.imageIndex).ROInumber =  num2str(handles.ROIindex);
    handles.finalROIs(handles.imageIndex).type = handles.ROIregionType;
    handles.finalROIs(handles.imageIndex).length = handles.ROIlength;
    handles.finalROIs(handles.imageIndex).straightlines = handles.ROIlines;
    handles.finalROIs(handles.imageIndex).straightlineDistances = linesDist;
    handles.finalROIs(handles.imageIndex).straightlineCoord = linesCoord;
    handles.finalROIs(handles.imageIndex).box = handles.b;
    handles.finalROIs(handles.imageIndex).skel = handles.skel;
    handles.finalROIs(handles.imageIndex).perimeter = bwperim(handles.ROIpath,8);
    handles.finalROIs(handles.imageIndex).path = handles.ROIpath;
    handles.finalROIs(handles.imageIndex).channel = handles.channelChoice;
    handles.finalROIs(handles.imageIndex).virtualisation = handles.virtualisationlevel;
    handles.finalROIs(handles.imageIndex).grayscale = handles.grayscalevalue;
    handles.finalROIs(handles.imageIndex).virtImage = handles.virtualisedImage(handles.b(2):handles.b(2)+handles.b(4), handles.b(1):handles.b(1)+handles.b(3),1);
    handles.finalROIs(handles.imageIndex).regioncoordinates = table(handles.bigx,handles.bigy,'VariableNames',{'x','y'});
    % handles.finalROIs(handles.imageIndex).
    
    
    handles.ROIindex = handles.ROIindex+1;
    set(handles.undobutton,'Enable','on');
end
% set(handles.regionbutton,'Enable','off');

% cla
% imshow(handles.region)
%
% set(handles.restorebutton,'Enable','on');
% set(handles.activatemicrotubules,'Enable','on');
set(handles.calculateallbutton,'Enable','on');
set(hObject, 'Enable', 'on');
guidata(hObject, handles);
function microscopebutton_Callback(hObject, ~, handles)
% This function is going to popup ALFRED_microscopespecs and retrieve the
% data from there. Then, it saves the choice to the entire image, not just
% to any particular ROI.
try
    delete(handles.scalebarplot)
    delete(handles.scaletext);
catch
end
microscopeSpecs=ALFRED_microscopespecs;
hasChanged = str2double(string(microscopeSpecs(1)));
if strcmp(string(microscopeSpecs(8)), 'On')
    try
        numvalue = handles.metadata{1,1}{handles.currentImage}.getPixelsPhysicalSizeX(0).value();
        unit = handles.metadata{1,1}{handles.currentImage}.getPixelsPhysicalSizeX(0).unit().getSymbol();
    catch
        errordlg('Metadata not available');
        set(handles.imageSpecDisplay,'String','N/A');
        hasChanged = 3;
    end
else
    % handles.change, handles.microscope,handles.objective,handles.binning,handles.optovar,num2str(handles.numvalue),handles.unit,
    if hasChanged == 2
        numvalue = str2double(string(microscopeSpecs(6)));
        unit = string(microscopeSpecs(7));
    elseif hasChanged == 1
        microscope = string(microscopeSpecs(2));
        objective = string(microscopeSpecs(3));
        binning = string(microscopeSpecs(4));
        optovar = string(microscopeSpecs(5));
        
        [~, ~, scalingfactors] = xlsread('documents\scalingfactors.xlsx','corrected scaling factors','A2:P24');
        scalingfactors = string(scalingfactors);
        scalingfactors(ismissing(scalingfactors)) = '';
        half = scalingfactors(:,1:size(scalingfactors,2)/2);
        if ~size(find(half==microscope),1)
            half = scalingfactors(:,size(scalingfactors,2)/2:end);
        end
        [obrow,~] = find(half == objective);
        [oprow,~] = find(half == optovar);
        [~,col] = find(half == binning);
        row = intersect(obrow,oprow);
        numvalue = str2double(half(row,col));
        unit = 'µm';
    end
end
if hasChanged ~= 0
    handles.conversionFactor = double(numvalue);
    handles.physicalUnit = string(unit);
    set(handles.imageSpecDisplay,'String',num2str(double(numvalue))+string(unit));
end
hold on
y=size(handles.image,1);
x=size(handles.image,2);
sizeOfBar = 10/handles.conversionFactor;
textX = floor((x/20+sizeOfBar/2));
textY = y-150;
handles.scalebarplot=plot([floor(x/20) floor((x/20+sizeOfBar))], [y-100 y-100], '-w', 'LineWidth',5);
handles.scaletext = text(textX,textY, char(num2str(10) + handles.physicalUnit), 'HorizontalAlignment','center','Color','w');
% t.fontSize = 8;

hold off

handles.finalROIs(handles.imageIndex).conversionFactor = handles.conversionFactor;
handles.finalROIs(handles.imageIndex).physicalUnit = handles.physicalUnit;

guidata(hObject, handles);


function previousButton_Callback(hObject, ~, handles)
%Go to previous image, plot the already selected bounding boxes
handles.currentImage = handles.currentImage - 1;
handles.image = handles.imageGroup{handles.currentImage};
% handles.oldImage = handles.image;

% Inverse image: should be an independent function but goddamn matlab
% has the handles that ruin everything.
if get(handles.inverseImage,'Value')
    handles.oldImage = handles.image;
    for i = 1:size(handles.image,3)
        handles.image(:,:,i) = imcomplement(handles.image(:,:,i));
    end
    
else
    handles.image = handles.oldImage;
end

try
    imshow(handles.image(:,:,handles.channelChoice));
catch
    if size(handles.image,3)>3
        imshow(handles.image(:,:,1));
        text(0,100,'Channel 1','Color','white');
        handles.oldImage = handles.image(:,:,1);
    else
        imshow(handles.image)
%         handles.oldImage = handles.image;
    end
end


handles.channelNumber = size(handles.image,3);
channelOpt = ['All' newline num2str(linspace(1,handles.channelNumber, handles.channelNumber),'%d\n')];
set(handles.channelselectmenu,'String',channelOpt);
set(handles.textdisplay,'String',['Image Number ', num2str(handles.currentImage) ,' ', num2str(size(handles.image,1)), 'x', num2str(size(handles.image,2)) ]);
checkuint(hObject, handles)
handles.selectedRegion = handles.image;
handles.region = handles.image;
cla(handles.display,'reset');
if size(handles.image,3)>3
    imshow(handles.image(:,:,1));
    text(0,100,'Channel 1','Color','white');
else
    try
        %         imshow(handles.image(:,:,handles.channelChoice));
        handles.region = handles.image(:,:,handles.channelChoice);
        %         virtualisationbutton_Callback(@virtualisationbutton, 0, handles);
        handles = virtualising(handles);
        if handles.imagePanel
            applygrayscalethreshold_Callback(@applygrayscalethreshold, 0, handles);
            
        end
    catch
        imshow(handles.image)
    end
end
checkPrint(handles.figure1,handles)
if handles.currentImage == 1
    set(handles.previousButton,'Enable','off');
end
set(handles.undobutton,'Enable','off');
set(handles.continuebutton,'Enable','on');
guidata(hObject,handles)
function undobutton_Callback(hObject, ~, handles)
%Undo previous action (i.e. eliminate the last ROI that was saved)
delete(handles.plotroi)
handles.finalROIs(handles.imageIndex)=[];
handles.imageIndex = handles.imageIndex-1;
if handles.ROIindex>0
    handles.ROIindex = handles.ROIindex-1;
end
set(handles.undobutton,'Enable','off');
guidata(hObject,handles)
function continuebutton_Callback(hObject, ~, handles)
%This button will save any changes and move on to the next image if there
%is. If there isn't, then you can only calculate all.
handles.currentImage = handles.currentImage + 1;
handles.image = handles.imageGroup{handles.currentImage};

% Inverse image: should be an independent function but goddamn matlab
% has the handles that ruin everything.
if get(handles.inverseImage,'Value')
    handles.oldImage = handles.image;
    for i = 1:size(handles.image,3)
        handles.image(:,:,i) = imcomplement(handles.image(:,:,i));
    end
    
else
    handles.image = handles.oldImage;
end

try
    imshow(handles.image(:,:,handles.channelChoice));
catch
    if size(handles.image,3)>3
        imshow(handles.image(:,:,1));
        text(0,100,'Channel 1','Color','white');
        handles.oldImage = handles.image(:,:,1);
    else
        imshow(handles.image)
%         handles.oldImage = handles.image;
    end
end

handles.channelNumber = size(handles.image,3);
channelOpt = ['All' newline num2str(linspace(1,handles.channelNumber, handles.channelNumber),'%d\n')];
set(handles.channelselectmenu,'String',channelOpt);
set(handles.textdisplay,'String',['Image Number ', num2str(handles.currentImage) ,' ', num2str(size(handles.image,1)), 'x', num2str(size(handles.image,2)) ]);
checkuint(hObject, handles)
handles.selectedRegion = handles.image;
handles.region = handles.image;

cla(handles.display,'reset');

if size(handles.image,3)>3
    imshow(handles.image(:,:,1));
    text(0,100,'Channel 1','Color','white');
else
    try
        %         imshow(handles.image(:,:,handles.channelChoice));
        handles.region = handles.image(:,:,handles.channelChoice);
        handles = virtualising(handles);
        %virtualisationbutton_Callback(@virtualisationbutton, 0, handles);
        
        if handles.imagePanel
            applygrayscalethreshold_Callback(@applygrayscalethreshold, 0, handles);
            
        end
    catch
        imshow(handles.image)
    end
end

checkPrint(handles.figure1,handles)

if handles.currentImage == handles.totalNumberImages || handles.currentImage == length(handles.imageGroup)
    set(handles.continuebutton,'Enable','off');
end
set(handles.undobutton,'Enable','off');
set(handles.previousButton,'Enable','on');
guidata(hObject,handles)
function calculateallbutton_Callback(hObject, ~, handles)
%Calculate every region and save to a file
global DEBUG;

saveAnswer = savePopUp;
finalName = ['results/' num2str(yyyymmdd(datetime)) '-' num2str(hour(datetime)) num2str(minute(datetime)) num2str(round(second(datetime))) '_' handles.fileName{handles.currentImage} ];

if string(saveAnswer(1))=="Yes"
    
    typeROI = 'Axon';
    [~,indexAxon] = findValue(handles.finalROIs,typeROI);
    finalLengths = zeros(numel(indexAxon),1);
    
    typeROI = 'MT disorganisation';
    [~,indexMT] = findValue(handles.finalROIs,typeROI);
    finalAreasDisorg = zeros(numel(indexMT),1);
    finalAreaMT =  zeros(numel(indexMT),1);
    finalDensity = zeros(numel(indexMT),1);
    finalEccentricity = zeros(numel(indexMT),1);
    finalGeneralCurvatures = struct([]);
    finalPointCurvatures = struct([]);
    
    %Calculations on axons:
    if indexAxon > 0
        f = waitbar(0,'1','Name','Calculating axon results',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        for i = numel(indexAxon) % Iterate on all the ROIs that are axons
            if getappdata(f,'canceling')
                break
            end
            waitbar(i/numel(indexAxon),f,sprintf('ROI %1d of %1d',i,numel(indexAxon)))
            finalLengths(i) = handles.finalROIs(indexAxon(i)).length*handles.conversionFactor;
            handles.finalROIs(indexAxon(i)).unitLength = handles.finalROIs(indexAxon(i)).length*handles.conversionFactor;
        end
        delete(f)
    else
        finalLengths = 0;
    end
    checkPrint(handles.figure1,handles)
    finalName = [finalName '.mat'];
    tempROIs = handles.finalROIs; 
    save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
    % Add straightline distances/sizes to save also as another variable,
    % i'm not sure how but we don't need it now!
    
    %Calculations on MT disorganisation:
    % NOTE: curvatures are calculated using two methods: the Gaussian
    % window method (dependent on resolution but good for fine grained
    % solutions), and the Fourier method, good for course grained
    % solutions and general curvature of the image.
    imagesWithDisorganisation = [];
    
    if indexMT>0
        %             totalCurvatures=[];
        f = waitbar(0,'1','Name','Calculating MT disorganisation results',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        for i = 1:numel(indexMT)
            
            try
                waitbar(i/numel(indexMT),f,sprintf('ROI %1d of %1d',i,numel(indexMT)))
                tempV = handles.finalROIs(indexMT(i)).virtImage;
                handles.bw = im2bw(tempV,[],0.1); %#ok<IM2BW>
                handles.bw = ~handles.bw;
                %             handles.bw = imclearborder(handles.bw);
                [labeledImage, ~] = bwlabel(handles.bw);
                handles.measurements = regionprops(labeledImage,'Area','Eccentricity');
                handles.table = sortrows(struct2table(handles.measurements),'Area','descend');
                
                if isempty(intersect(imagesWithDisorganisation,handles.finalROIs(indexMT(i)).imageName))
                    imagesWithDisorganisation = [imagesWithDisorganisation ; handles.finalROIs(indexMT(i)).imageName]; %#ok<AGROW>
                end
                
                %Calculating average values from both area and eccentricity
                sumEcc=0;
                sumArea = 0;
                tot = 0;
                for j = 1:height(handles.table)
                    if handles.table{j,1}>5
                        sumArea = sumArea + handles.table{j,1};
                        sumEcc = sumEcc + handles.table{j,2};
                        tot = tot + 1;
                    end
                end
                if ~isempty(handles.finalROIs(indexMT(i)).conversionFactor) % O Nuno editou isto, não tenho a certeza que esteja a fazer o que queres
                    conversionFactor = handles.finalROIs(indexMT(i)).conversionFactor;
                else
                    conversionFactor = handles.conversionFactor;
                end
                if DEBUG; disp('before calculations'); end
                
                finalAreasDisorg(i) = sumArea/tot*conversionFactor; % Average Disorganisation area - shouldn't it be the perimeter???
                
                se = strel('square',4);
                dilated = imdilate(handles.finalROIs(indexMT(i)).path, se);
                filled = imfill(dilated);
                area = sum(sum(filled));
                
                finalEccentricity(i) = sumEcc/tot*conversionFactor; % Average Eccentricity
                handles.adjusted = handles.finalROIs(indexMT(i)).virtImage>(handles.finalROIs(indexMT(i)).grayscale/255);
                finalAreaMT(i) = area*conversionFactor; %Area of microtubules in the image
                
                
                %                 handles.imageCurvature = handles.imageCurvature.*(handles.imageCurvature<0.2);
                %                 handles.imageCurvature(handles.imageCurvature==0) = [];
                
%                 [~, ~, resultsFourier] = curvatureFourier(handles.finalROIs(indexMT(i)).regioncoordinates.x,handles.finalROIs(indexMT(i)).regioncoordinates.y);
%                 [~, ~, resultsGauss] = curvatureGaussTuner(handles.finalROIs(indexMT(i)).regioncoordinates.x,handles.finalROIs(indexMT(i)).regioncoordinates.y);
                skeleton = curvaturePreProcessing(handles.finalROIs(indexMT(i)).regioncoordinates.x,handles.finalROIs(indexMT(i)).regioncoordinates.y);
                tol=0.2;
                [N,~,~, radiusSpline, dsdt, ~] = curvatureSpline(skeleton,tol);
                resultsFourier = [];
%                 resultsSpline = [];
                finalGeneralCurvatures{i} = 1./(resultsFourier.*handles.conversionFactor);
                finalPointCurvatures{i} = 1./(radiusSpline{1}(1:1e-3:N{1}).*handles.conversionFactor);
                assignin('base','newArcLength',dsdt{1}(1:1e-3:N{1}))
                
                if DEBUG; disp('after calculations'); end
                %                 totalCurvatures = [totalCurvatures handles.imageCurvature]; %#ok<AGROW>
                
                handles.finalROIs(indexMT(i)).areasDisorg = finalAreasDisorg(i);
                handles.finalROIs(indexMT(i)).areaMT = finalAreaMT(i);
                handles.finalROIs(indexMT(i)).generalCurvatures = finalGeneralCurvatures{i};
                handles.finalROIs(indexMT(i)).individualCurvatures = finalPointCurvatures{i};
                %             handles.finalROIs(indexMT(i)).density =
                handles.finalROIs(indexMT(i)).eccentricity = finalEccentricity(i);
                % ADD FIGURE WITH THE CURVATURE OR JUST SAVE THE FILE?
            catch ME
                if DEBUG; disp(['WARNING: ',ME.identifier]); end
                continue
            end
        end
        delete(f)
    else
        finalAreasDisorg = 0;
        finalAreaMT = 0;
        finalDensity = 0;
        finalEccentricity = 0;
    end
    checkPrint(handles.figure1,handles)
    %     finalName = [finalName '.mat'];
    tempROIs = handles.finalROIs;
    save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
    
    
    
    f = waitbar(0,'1','Name','Calculating MDI results',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    for i = 1:size(imagesWithDisorganisation,1)
        try
            waitbar(i/size(imagesWithDisorganisation,1),f,sprintf('Image %1d of %1d',i,size(imagesWithDisorganisation,1)))
            [~,indexImages,~] = findValue(handles.finalROIs,imagesWithDisorganisation(i,:));
            [~,numAxon,~] = findValue(handles.finalROIs(indexImages),'Axon');
            [~,numMT,~] = findValue(handles.finalROIs(indexImages),'MT disorganisation');
            % This gives me a matrix with the overlap of the
            % disorganisation with the axons. Rows: axon, Columns: MTs
            ratios = zeros(numel(numAxon),numel(numMT));
%             'entering the cycles'
            for j = 1:numel(numAxon)
                for k = 1:numel(numMT)
%                     'calculating bbox overlap'
                    ratios(j,k)=bboxOverlapRatio(handles.finalROIs(numAxon(j)).box,handles.finalROIs(numMT(k)).box,'Min');
                    if ratios(j,k)<0.05
                        ratios(j,k)=0;
                        continue
                    end
%                     'checking skeleton overlap'
                    if ~skeletonOverlap(handles.finalROIs(numAxon(j)).box,handles.finalROIs(numMT(k)).box,handles.finalROIs(numAxon(j)).skel,handles.finalROIs(numMT(k)).skel) %if skeleton overlap is 0
                        ratios(j,k)=0;
                    end
%                     'left skeleton overlap'
                end
            end
%             'before'
            [~,indexes] = max(ratios); %gives me the max ratio per disorganisation (column)
            uniqueIndex = unique(indexes);
            for j = 1:numel(uniqueIndex)
%                 'calculating mdi'
                
                disorgIndex = find(indexes==uniqueIndex(j));
                finalAreaMT = 0;
                for k = 1:numel(disorgIndex)
                    finalAreaMT = finalAreaMT+handles.finalROIs(numMT(disorgIndex(k))).areaMT; % tirar o +1 e por + disorg(index).perimeter
                end
                % MDI = perimeterMT/length(axon)
                handles.finalROIs(numAxon(j)).MDI = finalAreaMT/handles.finalROIs(numAxon(j)).length;
            end
        catch
%             'I DIED'
            continue
        end
        
    end
    delete(f)
    checkPrint(handles.figure1,handles)
    %     finalName = [finalName '.mat'];
    tempROIs = handles.finalROIs;
    save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
    
    finalStraightLengths =  cell(numel(handles.finalROIs)-1,1);
    finalStraightPos = cell(numel(handles.finalROIs)-1,1);
    finalStraightnessIndex = cell(numel(handles.finalROIs)-1,1);
    %         handles.finalROIs(handles.imageIndex).straightlineDistances
    for i = 2:numel(handles.finalROIs)
        try
            if ~isempty(handles.finalROIs(i).conversionFactor) % O Nuno editou isto, não tenho a certeza que esteja a fazer o que queres
                conversionFactor = handles.finalROIs(i).conversionFactor;
            else
                conversionFactor = handles.conversionFactor;
            end
            finalStraightLengths{i-1} = handles.finalROIs(i).straightlineDistances.*conversionFactor;
            finalStraightPos{i-1} = handles.finalROIs(i).straightlineCoord;
            
            if strcmp(handles.finalROIs(i).type,'MT disorganisation') % THESE ARE INDEXES, i don't need the conversion factor
                finalStraightnessIndex{i-1} = sum(handles.finalROIs(i).straightlineDistances)./handles.finalROIs(i).areaMT;
            else
                finalStraightnessIndex{i-1} = sum(handles.finalROIs(i).straightlineDistances)./handles.finalROIs(i).length;
            end
        catch
            continue
        end
    end
    
%     csvname = [finalName '.csv'];
    
    % cellROIs = struct2cell(handles.finalROIs);
    i=2; % Index for what the user wants
    j=1; % index for the finalResults struct
    finalResults = struct([]);
    f = waitbar(0,'1','Name','Saving',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    checkPrint(handles.figure1,handles)
    while i<=size(saveAnswer,2)
        finalResults(j).name = saveAnswer(i);
        switch string(saveAnswer(i))
            case 'General Curvature'
                finalResults(j).value = finalGeneralCurvatures;
            case 'Individual Curvature'
                finalResults(j).value = finalPointCurvatures;
            case 'Axonal Length'
                finalResults(j).value = finalLengths;
            case 'Straight Segment Lengths'
                finalResults(j).value = finalStraightLengths;
            case 'Straight Segment Positions'
                finalResults(j).value = finalStraightPos;
            case 'Eccentricity'
                finalResults(j).value = finalEccentricity;
                %                 case 'MT Disorganisation Index'
                %                     finalResults(j).value = ;
            case 'Straightness Index'
                finalResults(j).value = finalStraightnessIndex;
            case 'MT Area'
                finalResults(j).value = finalAreaMT;
            case 'Swelling Area'
                finalResults(j).value = finalAreasDisorg;
            case 'Density'
                finalResults(j).value = finalDensity;
        end
        i = i+1;
        j = j+1;
        
    end
    delete(f)
    checkPrint(handles.figure1,handles)
    %         finalResults = struct('AxonalLengths',finalLengths,'AreaOfDisorganisation',finalAreasDisorg,'AreaMicrotubules',finalAreaMT,'Density',finalDensity,'Eccentricity',finalEccentricity,'Curvatures',finalGeneralCurvatures);
    %         c=struct2cell(finalResults);
    %         xlsname = [finalName '.csv'];
    %                 xlswrite(xlsname, vertcat(c{:}));
    %or if you want the field names:
    %         xlswrite(xlsname, [fieldnames(finalResults), c]);
    
%     struct2csv(finalResults,csvname);
    
    msgbox('The requested files have been saved.', 'Saved');
    
    set(handles.undobutton,'Enable','off');
    %catch ex
    %    disp(ex)
    %    warndlg('No conversion factor was selected. Use Microscope Specs')
    %end
    
end
% finalName = [finalName '.mat'];
tempROIs = handles.finalROIs;
save(finalName,'tempROIs'); % Saves everything into a file that can be opened later

guidata(hObject,handles)

% ------------------------------------------------------------------------------------------------------------------------------------------------
% Toggle Buttons Selection Functions
function regionselectpanel_SelectionChangedFcn(hObject, ~, handles)
handles.ROICutChoice = get(get(handles.regionselectpanel,'SelectedObject'),'String');
set(handles.regionbutton,'Enable','on');

guidata(hObject, handles);

% ------------------------------------------------------------------------------------------------------------------------------------------------
% List Menu Selection Functions
function channelselectmenu_Callback(hObject, ~, handles)

contents = cellstr(get(hObject,'String'));
value = contents{get(hObject,'Value')};
if strcmp(value,'All')
    set(handles.virtualisationbutton,'Enable','off');
    if size(handles.image,3) <= 3
        handles.region = handles.image;
    else
        imshow(handles.image(:,:,1));
        text(0,100,'Channel 1','Color','white');
        
    end
else
    handles.channelChoice = str2double(value);
    handles.region = handles.image(:,:,handles.channelChoice);
    txt = ['Channel ' num2str(handles.channelChoice)];
    text(0,100,txt,'Color','white');
    set(handles.virtualisationbutton,'Enable','on');
end

% virtualisationbutton_Callback(hObject, 0, handles);
% ALFRED_mainpanel('virtualisationbutton_Callback',handles.virtualisationbutton,[],handles);
handles.region = double(handles.region);
cla(handles.display,'reset');
%see comments on the vesselness2D file to change the parameters if necessary. True/False: are vessels bright on dark background or dark on dark
tempV=vesselness2D(handles.region(:,:,1), 1:get(handles.virtualisationscaleslider,'Value'), [1;1;1],1,true);

handles.virtualisedImage = tempV.*(tempV>0.005);
handles.virtualisationlevel = get(handles.virtualisationscaleslider,'Value');
imshow(handles.virtualisedImage);

% set(handles.calculateEnd,'Enable','on');
checkPrint(handles.figure1,handles)
handles.oldImage = handles.virtualisedImage;
set(handles.applygrayscalethreshold,'Enable','on');
% imshow(handles.region)
checkPrint(handles.figure1,handles)
guidata(hObject,handles)

% ------------------------------------------------------------------------------------------------------------------------------------------------
% Edit Text Box Functions
function virtualisationscalebox_Callback(hObject, ~, handles)
handles.virtualScaleValue = str2double(get(hObject,'String'));
if isnan(handles.virtualScaleValue)
    handles.virtualScaleValue=3;
else
    if handles.virtualScaleValue>get(handles. virtualisationscaleslider,'Max')
        handles.virtualScaleValue = get(handles. virtualisationscaleslider,'Max');
    elseif handles.virtualScaleValue<get(handles. virtualisationscaleslider,'Min')
        handles.virtualScaleValue = get(handles. virtualisationscaleslider,'Min');
    end
end
set(handles.virtualisationscalebox,'String',num2str(handles.virtualScaleValue));
if handles.virtualScaleValue ~= get(handles. virtualisationscaleslider,'Value')
    set(handles. virtualisationscaleslider,'Value',handles.virtualScaleValue)
end
guidata(hObject, handles);
function grayscalethresholdbox_Callback(hObject, ~, handles)
handles.GSThresholdValue = str2double(get(hObject,'String'));
if isnan(handles.GSThresholdValue)
    handles.GSThresholdValue=0;
else
    if handles.GSThresholdValue>get(handles.grayscalethresholdslider,'Max')
        handles.GSThresholdValue = get(handles.grayscalethresholdslider,'Max');
    elseif handles.GSThresholdValue<get(handles.grayscalethresholdslider,'Min')
        handles.GSThresholdValue = get(handles.grayscalethresholdslider,'Min');
    end
end
set(handles.grayscalethresholdbox,'String',num2str(handles.GSThresholdValue));
if handles.GSThresholdValue ~= get(handles.grayscalethresholdslider,'Value')
    set(handles.grayscalethresholdslider,'Value',handles.GSThresholdValue)
end
guidata(hObject,handles)

function imagenormalisecheck_Callback(~, ~, ~)
% Only used in checkuint function, to check whether the user whats a
% normalised image.

% ------------------------------------------------------------------------------------------------------------------------------------------------
% Slider Functions
function virtualisationscaleslider_Callback(hObject, ~, handles)
temp = get(handles.virtualisationscaleslider,'Value');
if handles.virtualScaleValue ~= temp
    set(handles.virtualisationscalebox,'String',temp)
    handles.virtualScaleValue = temp;
end
guidata(hObject, handles);
function grayscalethresholdslider_Callback(hObject, ~, handles)
temp = get(handles.grayscalethresholdslider,'Value');
if handles.GSThresholdValue ~= temp
    set(handles.grayscalethresholdbox,'String',temp)
    handles.GSThresholdValue = temp;
end
guidata(hObject, handles);

% ------------------------------------------------------------------------------------------------------------------------------------------------
% Create Functions
function virtualisationsensitivitybox_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function virtualisationscalebox_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function grayscalethresholdbox_CreateFcn(hObject, ~,  ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function grayscalemaxcontrastbox_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function grayscalemincontrastbox_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function grayscalethresholdslider_CreateFcn(hObject, ~, handles)
set(hObject,'Max',255,'Min',0,'Value',0)
set(hObject,'SliderStep',[1/255 1/255])
handles.GSThresholdValue=0;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject,handles)
function virtualisationscaleslider_CreateFcn(hObject, ~, handles)
set(hObject,'Max',5,'Min',1,'Value',3)
set(hObject,'SliderStep',[1/(5-1) 1])
handles.virtualScaleValue = 3;
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function regionselectpanel_CreateFcn(hObject, ~, ~)
set(hObject,'selectedobject','')
function channelselect_CreateFcn(hObject, ~, ~)
set(hObject,'selectedobject','')
function channelselectmenu_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, ~)
% Hint: delete(hObject) closes the figure
cancelrequest = cancelpopup;
if strcmp(cancelrequest,'Yes')
    if isequal(get(hObject, 'waitstatus'), 'waiting')
        % The GUI is still in UIWAIT, us UIRESUME
        uiresume(hObject);
    else
        % The GUI is no longer waiting, just close it
        delete(hObject);
    end
end


% --- Executes on button press in originalImageButton.
function originalImageButton_Callback(hObject, ~, handles)

% IF BUTTOn TEXT == Original: SHOW ORIGINAL, CHANGE TEXT
% ELSEIF BUTTON TEXT == Return: SHOW WHAT IT WAS, CHANGE TEXT
if strcmp(handles.originalImageButton.String, 'View Original Image')
    imshow(handles.image);
    set(handles.originalImageButton,'String','Return to image')
else
    imshow(handles.oldImage);
    set(handles.originalImageButton,'String','View Original Image')
end
checkPrint(handles.figure1,handles)
guidata(hObject, handles);


% --- Executes on button press in resetregionsbutton.
function resetregionsbutton_Callback(hObject, ~, handles)
cancelrequest = cancelpopup;
if strcmp(cancelrequest,'Yes')
    imageName =  handles.fileName{handles.currentImage};
    [~,indexAxon] = findValue(handles.finalROIs,imageName);
    for i=1:numel(indexAxon)
        handles.finalROIs(indexAxon(i)-i+1) = []; % the -i+1 is there to guarantee that we are eliminating the correct one. It guarantees the update of the indexes as we are deleting lines.
    end
    cla(handles.display,'reset');
    imshow(handles.maskedImage);
    try
        hold on
        y=size(handles.image,1);
        x=size(handles.image,2);
        sizeOfBar = handles.conversionFactor * 10;
        textX = floor((x/20+sizeOfBar/2));
        textY = y-150;
        handles.scalebarplot=plot([floor(x/20) floor((x/20+sizeOfBar))], [y-100 y-100], '-w', 'LineWidth',5);
        handles.scaletext = text(textX,textY, char(num2str(handles.conversionFactor) + handles.physicalUnit), 'HorizontalAlignment','center','Color','w');
        % t.fontSize = 8;
        
        hold off
    catch
    end
    handles.oldImage = handles.maskedImage;
    set(handles.freehand,'Enable','on');
    set(handles.rectangle,'Enable','on');
    handles.finalROIs(handles.imageIndex).imageName = handles.fileName{handles.currentImage};
    handles.finalROIs(handles.imageIndex).ROInumber =  handles.baseImageIndex;
    handles.finalROIs(handles.imageIndex).type = 'Full';
    handles.ROIindex = 1;
end
guidata(hObject, handles);


% --- Executes on button press in inverseImage.
function inverseImage_Callback(hObject, ~, handles)
% Complements the image (basically if it's black on white it becomes white
% on black).
if get(hObject,'Value')
%     'hey'
    handles.oldImage = handles.image;
    for i = 1:size(handles.image,3)
        handles.image(:,:,i) = imcomplement(handles.image(:,:,i));
    end
    
else
    handles.image = handles.oldImage;
end

try
    imshow(handles.image(:,:,handles.channelChoice));
catch
    if size(handles.image,3)>3
        imshow(handles.image(:,:,1));
        text(0,100,'Channel 1','Color','white');
        handles.oldImage = handles.image(:,:,1);
    else
        imshow(handles.image)
%         handles.oldImage = handles.image;
    end
end
guidata(hObject, handles);