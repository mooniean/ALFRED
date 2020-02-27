function varargout = ALFRED_lengthpopup(varargin)
% ALFRED_LENGTHPOPUP MATLAB code for ALFRED_lengthpopup.fig
%         Panel opening after selecting the region in the Main Image Panel, with the skeleton of the region overlapping the black and white real image.
%         Initially, the user can select from the dropdown menu which type of region is the one presented. In this version, only Microtubule Disorganisation and Axon are available, 
%         and contact the author if you'd like to have something different.
% 
%         The designation of the region will then Enable the buttons that are essential to the region at hand. Either option will enable the Start Selection.
% 
%         The user is asked to select the beginning and end of the region. This selection is done by clicking on the skeleton displayed. 
% 
%         Two yellow dots will appear on the skeleton in the closest point to the coordinates the user selected. The region will automatically appear and the user can immediately understand whether extra points
%         need to be added, with the Add Points. The amount of points to be added at one particular step can be changed in the editable text box. 
% 
%         There are two types of extra points, inner and border points, and they are available only for axons and MT disorganisation, respectively. 
%         The inner points can be used when the intended axonal path needs to pass a certain point, whether because there is a gap or because the shortest path in the algorithm is not the correct biological one. 
%         The border points will allow the user to trim the disorganisation region and remove any unwanted branches. 
% 
%         In case of an axonal region, the user can Add Lines to circumvent any imaging problems (e.g., as sometimes instead of a continuous line, the image has a dotted one).
% 
%         The Undo button will retrieve the last added points and if the user is still not convinced.
%         Reset Path will delete any selected points and the user can re-start the selection.
%         If everything is in order, the user can press Apply and then Save and Close, which will close the popup. If the user doesn't want to save anything,
%         Cancel will close the popup without saving at any point in the analysis.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ALFRED_lengthpopup

% Last Modified by GUIDE v2.5 11-Mar-2019 10:56:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ALFRED_lengthpopup_OpeningFcn, ...
    'gui_OutputFcn',  @ALFRED_lengthpopup_OutputFcn, ...
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


% --- Executes just before ALFRED_lengthpopup is made visible.
function ALFRED_lengthpopup_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
set(handles.figure1, 'Name', 'ALFRED - Region Definition');
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1,'javaframe');
jIcon=javax.swing.ImageIcon('icons\alfred_50x50.png');
jframe.setFigureIcon(jIcon); % Setting icon of the program

handles.image = varargin{1}; % Skeletonized region! not the image itself
handles.originalImage = varargin{2}; % Crop of the original image
handles.g = binaryImageGraph(handles.image);
handles.amountofpoints = 1;
handles.totalpoints = 2;
handles.numpeaks = 10;
handles.linematrix = [];

% handles.nodenumber = [1 2];
handles.regiontype = 'Axon';
handles.flag = 1;
handles.nodenumber = [];
handles.radius = 10;
handles.display=axes('Parent',handles.displaypanel);

handles.maxJump = round(max(size(handles.image))/2);
handles.temporarySkel = double(handles.image);
handles.temporarySkel(handles.temporarySkel == 0) = NaN;

hold on
imshow(handles.originalImage)
hskel = imshow(handles.temporarySkel);
set(hskel, 'AlphaData', ~isnan(handles.temporarySkel))

hold off
handles.finalX = [];
handles.finalY = [];

set(handles.radiusSlider,'Max',handles.maxJump,'Min',0,'Value',10)
set(handles.radiusSlider,'SliderStep',[1/handles.maxJump 1/handles.maxJump])
set(handles.maxJumpText,'String',handles.maxJump);


set(handles.radiusText,'String',handles.radius);
set(handles.selectedpoints,'String','');
set(handles.estimatedlengthoutput,'String','');



set(handles.addpointsbutton,'Enable','off');
set(handles.undobutton,'Enable','off');
set(handles.addlinebutton,'Enable','off');
set(handles.inputamountofpointstoadd,'Enable','off');
set(handles.startselectionbutton,'Enable','off');
set(handles.resetpathbutton,'Enable','off');
% set(handles.pointsgroup,'Enable','off');
set(handles.innerpoint,'Enable','off');
set(handles.borderpoint,'Enable','off');
set(handles.saveclosebutton,'Enable','off');
set(handles.entireimage,'Enable','off');
set(handles.applybutton,'Enable','off');
set(handles.radiusSlider,'Enable','off');
handles.graph = handles.g;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ALFRED_lengthpopup wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ALFRED_lengthpopup_OutputFcn(~, ~, handles)
% What do I want to be the outputs of this: type of region, path + length.

varargout = handles.output;
delete(handles.figure1);

function startselectionbutton_Callback(~, eventData, handles)
% Always start with two points: the endpoints!
% Run the length code. In case of error, popup to inform the user and the
% user can add more points.
if handles.flag
    handles.amountofpoints = 2;
    [handles.X, handles.Y] = ginputColour(handles.amountofpoints, [1 1 1]);
    for i=1:length(handles.X)
        [x,y,nodenumbertemp] = checkCoordinates(handles.X(i),handles.Y(i),handles.g.Nodes,handles.radius);
        handles.nodenumber = [handles.nodenumber nodenumbertemp];
        handles.finalX(end+1,1) = x;
        handles.finalY(end+1,1) = y;
    end
    handles.flag = 0;
end

hold on
handles.borderpoints = plot([handles.finalX],[handles.finalY],'y.','markers',12);
hold off

set(handles.undobutton,'Enable','on');

set(handles.applybutton,'Enable','on');
set(handles.startselectionbutton,'Enable','off');

applybutton_Callback(handles.applybutton, eventData, handles);


function addpointsbutton_Callback( ~, eventData, handles)
%Default : add only one point. Read the input on the edit box and add the
%number asked by the user. Maybe put a warning in case the input is more
%than 5. (User might prefer a line?)
%ginput has to be on this one, not on the lengthFinding function as it will
%only need the coordinates, and we save them here.
handles.amountofpoints = str2double(get(handles.inputamountofpointstoadd,'String'));
handles.totalpoints = handles.totalpoints + handles.amountofpoints;
[tempX, tempY] = ginputColour(handles.amountofpoints, [1 1 1]);
handles.X = [handles.X ; tempX];
handles.Y = [handles.Y ; tempY];

cla(handles.displaypanel,'reset')
hold on
imshow(handles.originalImage)
hskel = imshow(handles.temporarySkel);
set(hskel, 'AlphaData', ~isnan(handles.temporarySkel))
hold off

for i=1:length(tempX)
    [x,y,nodenumbertemp] = checkCoordinates(tempX(i),tempY(i),handles.graph.Nodes,handles.radius);
    handles.nodenumber = [handles.nodenumber nodenumbertemp];
    handles.finalX = [handles.finalX ; x];
    handles.finalY = [handles.finalY ; y];
end


newPoint = get(get(handles.pointsgroup,'SelectedObject'),'String');

switch newPoint
    case 'Inner Point'
        hold on
        handles.innerpoints = plot(handles.finalX(3:end),handles.finalY(3:end),'g.','markers',12);
        handles.borderpoints = plot(handles.finalX(1:2),handles.finalY(1:2),'y.','markers',12);
        hold off
    case 'Border Point'
%         try
            [~,~,firstUpdateGraph] = deleteGraphOutliers(handles.graph,handles.nodenumber(1),handles.nodenumber(end));
            
%         catch
%             [newNodes,~] = deleteGraphOutliers(handles.g,handles.nodenumber(1),handles.nodenumber(end));
%         end
%         handles.tempgraph = subgraph(handles.graph, newNodes);
%         graphx = [handles.graph.Nodes(handles.nodenumber(2),:).x handles.graph.Nodes(handles.nodenumber(end),:).x];
%         graphy = [handles.graph.Nodes(handles.nodenumber(2),:).y handles.graph.Nodes(handles.nodenumber(end),:).y];
        
%         nodeonetemp = find(table2array(handles.tempgraph.Nodes(:, 1)) == graphx(1) & table2array(handles.tempgraph.Nodes(:, 2)) == graphy(1));
%         nodetwotemp = find(table2array(handles.tempgraph.Nodes(:, 1)) == graphx(2) & table2array(handles.tempgraph.Nodes(:, 2)) == graphy(2));
        nodeonetemp = find(table2array(firstUpdateGraph.Nodes(:, 1)) == handles.finalX(2) & table2array(firstUpdateGraph.Nodes(:, 2)) == handles.finalY(2));
        nodetwotemp = find(table2array(firstUpdateGraph.Nodes(:, 1)) == handles.finalX(end) & table2array(firstUpdateGraph.Nodes(:, 2)) == handles.finalY(end));
        
        [~,setofcoordinates,secondUpdateGraph] = deleteGraphOutliers(firstUpdateGraph,nodeonetemp,nodetwotemp);
        
%         newNodes = [];
        
%         for i = 1:numel(newNodesTemp)
%             newNodes = [newNodes find(table2array(handles.graph.Nodes(:, 1)) == setofcoordinates(i).x & table2array(handles.graph.Nodes(:, 2)) == setofcoordinates(i).y)]; %#ok<AGROW>
%         end
%         handles.tempgraph = subgraph(handles.graph, newNodes);
        
        tempnodenumber = [];
        for i = 1:numel(handles.nodenumber)
            tempnodenumber = [tempnodenumber find(table2array(secondUpdateGraph.Nodes(:, 1)) == handles.finalX(i) & table2array(secondUpdateGraph.Nodes(:, 2)) == handles.finalY(i))]; %#ok<AGROW>
        end

%         [newNodes,setofcoordinates] = deleteGraphOutliers(handles.graph,handles.nodenumber(1),handles.nodenumber(end));
        handles.bigx = setofcoordinates.x;
        handles.bigy = setofcoordinates.y;
        hold on
        plot(handles.bigx,handles.bigy,'b.','LineWidth', 1);
        handles.borderpoints = plot([x handles.finalX(1) handles.finalX(2)],[y handles.finalY(1) handles.finalY(2)],'y.');
        hold off
%         handles.shortestpath = handles.path;
        handles.nodenumber = tempnodenumber;
        handles.graph = secondUpdateGraph;
end


set(handles.undobutton,'Enable','on');
set(handles.addpointsbutton,'Enable','on');
% set(handles.addlinebutton,'Enable','on');
set(handles.applybutton,'Enable','on');

% guidata(hObject,handles)
applybutton_Callback(handles.applybutton, eventData, handles);


function addlinebutton_Callback(hObject, ~, handles)
%Add a line instead of just two points. The idea is to get the coordinates
%and make a line with them, same way we do it to jump gaps.
handles.amountofpoints = 2;
[tempX, tempY] = ginputColour(handles.amountofpoints, [1 1 1]);

if ~handles.flag
    handles.totalpoints = handles.totalpoints + 2;
    handles.X = [handles.X ; tempX];
    handles.Y = [handles.Y ; tempY];
else
    handles.X = tempX;
    handles.Y = tempY;
end

for i=1:length(tempX)
    [x,y,nodenumbertemp] = checkCoordinates(tempX(i),tempY(i),handles.g.Nodes,handles.radius);
    handles.nodenumber = [handles.nodenumber  nodenumbertemp];
    handles.linematrix = [handles.linematrix  nodenumbertemp];
    %     assignin('base',handles.linematrix,'lines');
    handles.finalX = [handles.finalX ; x];
    handles.finalY = [handles.finalY ; y];
end
hold on
handles.hpoints = plot([handles.finalX],[handles.finalY],'g.','markers',12);
hold off
set(handles.undobutton,'Enable','on');
% assignin('base',handles.linematrix,'lines');
% assignin('base',handles.linematrix,'lines');
% Assume path as a straight line between these two points. Add a flag?
% Basically, check in the length finding whether there is a line between
% these. If line, jump to that one!
% Whats the most memory safe method? I can have a matrix of zeros and then
% I had a number according to the line, so that if there are different
% lines, I know that 1 meets 1 and 2 meets 2...


% try
%     [handles.totalSize,handles.path,handles.bigx,handles.bigy] = lengthFinding(hObject,handles);
% catch EX
%     EX.identifier
%     addlinebutton_Callback(handles.addpointsbutton,0,handles);
% end
% handles.finalPath=zeros(size(handles.image));
% for i = 1:length(handles.bigx)
%     handles.finalPath(handles.bigy(i),handles.bigx(i))=1;
% end
% [handles.lines,handles.longLine] = houghtest(handles.finalPath,handles.numpeaks);
%
% set(handles.estimatedlengthoutput,'String',num2str(handles.totalSize));
% set(handles.selectedpoints,'String',num2str(handles.totalpoints));
% hold on; plot(handles.bigx,handles.bigy,'r'); plot([handles.finalX],[handles.finalY],'g.'); hold off;
%
%
% set(handles.startselectionbutton,'Enable','off');
guidata(hObject,handles)

function undobutton_Callback(hObject, ~, handles)
% Undo the  last action performed.
% Checks how many points the user added last time, and remove them from the
% list
noValuesRemove = handles.amountofpoints;
while noValuesRemove>0
    noValuesRemove = noValuesRemove-1;
    handles.nodenumber(end) = [];
    handles.finalX(end) = [];
    handles.finalY(end) = [];
end
inner = false;
try
    delete(handles.innerpoints)
    inner = true;
catch
end

delete(handles.borderpoints)

if inner && (length(handles.finalX)>2)
    hold on
    handles.innerpoints = plot(handles.finalX(3:end),handles.finalY(3:end),'g.','markers',12);
    handles.borderpoints = plot(handles.finalX(1:2),handles.finalY(1:2),'y.','markers',12);
    hold off
elseif ~isempty(handles.finalX)
    hold on
    handles.borderpoints = plot([handles.finalX(1) handles.finalX(2)],[handles.finalY(1) handles.finalY(2)],'y.','markers',12);
    hold off
else
    handles.flag = 1;
    set(handles.startselectionbutton,'Enable','on');
end
set(handles.undobutton,'Enable','off');
guidata(hObject,handles)

function inputamountofpointstoadd_Callback(hObject, ~, handles)
%User input of how many points it wants to add. Default is 1.
handles.amountofpoints = str2double(get(hObject,'String'));
if isnan(handles.amountofpoints)
    handles.amountofpoints=1;
    set(handles.inputamountofpointstoadd,'String','1');
end
guidata(hObject, handles);
function typeList_Callback(hObject, ~, handles)
%Types of regions: axon, MT disorganisation.
contents = cellstr(get(hObject,'String'));
handles.regiontype = contents{get(hObject,'Value')};

% set(handles.typeList,'String',['Axon' ...
%     newline 'MT disorganisation']);

set(handles.startselectionbutton,'Enable','on');
set(handles.addpointsbutton,'Enable','on');
set(handles.inputamountofpointstoadd,'Enable','on');
% set(handles.pointsgroup,'Enable','on');

% set(handles.saveclosebutton,'Enable','on');

switch handles.regiontype
    case 'Select Region'
        set(handles.applybutton,'Enable','off');
        set(handles.entireimage,'Enable','off');
        set(handles.innerpoint,'Enable','off');
        set(handles.borderpoint,'Enable','off');
        set(handles.saveclosebutton,'Enable','off');
        set(handles.radiusSlider,'Enable','off');
    case 'Axon'
        set(handles.addlinebutton,'Enable','on');
        set(handles.innerpoint,'Enable','on');
        set(handles.borderpoint,'Enable','off');
        set(handles.entireimage,'Enable','off');
%         set(handles.saveclosebutton,'Enable','on');
        set(handles.radiusSlider,'Enable','on');
    case 'MT disorganisation'
        set(handles.entireimage,'Enable','on');
        set(handles.innerpoint,'Enable','off');
        set(handles.borderpoint,'Enable','on');
        set(handles.addlinebutton,'Enable','off');
        set(handles.radiusSlider,'Enable','off');
%         set(handles.saveclosebutton,'Enable','on');
end
guidata(hObject,handles)

function applybutton_Callback(hObject, ~, handles)
%Save the data selected, close the window.

try
    
    [handles.totalSize,handles.path,handles.bigx,handles.bigy] = lengthFinding(hObject,handles);
    if strcmp(handles.regiontype,'MT disorganisation')
        if get(handles.entireimage,'Value')
            handles.bigx=handles.graph.Nodes(:,:).x;
            handles.bigy=handles.graph.Nodes(:,:).y;
            hold on
            plot(handles.bigx,handles.bigy,'.b','LineWidth', 1);
            hold off
        else
            cla(handles.displaypanel,'reset')
            hold on
            imshow(handles.originalImage)
            hskel = imshow(handles.temporarySkel);
            set(hskel, 'AlphaData', ~isnan(handles.temporarySkel))
            hold off
            [~,setofcoordinates,cleanGraphUpdate] = deleteGraphOutliers(handles.graph,handles.nodenumber(1),handles.nodenumber(2));
            tempnodenumber = [];
            for i = 1:numel(handles.nodenumber)
                tempnodenumber = [tempnodenumber find(table2array(cleanGraphUpdate.Nodes(:, 1)) == handles.finalX(i) & table2array(cleanGraphUpdate.Nodes(:, 2)) == handles.finalY(i))]; %#ok<AGROW>
            end
            
            handles.bigx = setofcoordinates.x;
            handles.bigy = setofcoordinates.y;
            hold on
            plot(handles.bigx,handles.bigy,'.b','LineWidth', 1);
            handles.borderpoints = plot([handles.finalX(1) handles.finalX(2)],[handles.finalY(1) handles.finalY(2)],'y.');
            hold off

            handles.nodenumber = tempnodenumber;
            handles.graph = cleanGraphUpdate;
        end
    end
    
    handles.finalPath=zeros(size(handles.image));
    for i = 1:length(handles.bigx)
        handles.finalPath(handles.bigy(i),handles.bigx(i))=1;
    end
    [handles.lines,handles.longLine] = houghTest(handles.finalPath,handles.numpeaks);
    
    % Plot and show the user The size and how many points were user
    % selected
    set(handles.estimatedlengthoutput,'String',num2str(handles.totalSize));
    set(handles.selectedpoints,'String',num2str(handles.totalpoints));
    hold on; plot(handles.bigx,handles.bigy,'b.','LineWidth', 1); plot([handles.finalX],[handles.finalY],'g.'); hold off;
    
    % Average Line Values
    avg = 0;
    linesDist = [];
    linesCoord = cell(length(handles.lines),2);
    for i = 1:length(handles.lines)

        x1=handles.lines(i).point1(1);
        x2=handles.lines(i).point2(1);
        y1=handles.lines(i).point1(2);
        y2=handles.lines(i).point2(2);

        distancePoints = checkDistance(x1,x2,y1,y2);
        avg = avg + distancePoints;
        linesDist = [linesDist distancePoints]; %#ok<AGROW>
        linesCoord{i,1}= [x1 y1];
        linesCoord{i,2}= [x2 y2];
        
    end
    avg = avg/length(handles.lines);

    handles.output = {handles.regiontype, handles.finalPath, handles.totalSize,handles.lines,handles.longLine,handles.bigx,handles.bigy,linesDist,linesCoord};
    set(handles.saveclosebutton,'Enable','on');
catch EX
    warndlg('The selected points are not enough for a continuous path.')
    EX.identifier
    inner = false;
    try
        delete(handles.innerpoints)
        inner = true;
    catch
    end
    
    try
        delete(handles.borderpoints)
    catch
    end
    
    handles.flag = 1;
    
    if inner && (length(handles.finalX)>2)
        hold on
        handles.innerpoints = plot(handles.finalX(3:end),handles.finalY(3:end),'g.','markers',12);
        handles.borderpoints = plot(handles.finalX(1:2),handles.finalY(1:2),'y.','markers',12);
        hold off
    elseif ~isempty(handles.finalX)
        hold on
        handles.borderpoints = plot([handles.finalX(1) handles.finalX(2)],[handles.finalY(1) handles.finalY(2)],'y.','markers',12);
        hold off
    else
        handles.flag = 1;
        set(handles.startselectionbutton,'Enable','on');
    end
end

guidata(hObject,handles)

function cancelbutton_Callback(hObject, ~, handles)  %#ok<*DEFNU>
%Close without saving. Pop up asking if User is sure in cancelling, to
%avoid losing work.
cancelrequest = cancelpopup;
if strcmp(cancelrequest,'Yes')
%     handles.output = {'Close'};
    handles.output = {'None',0,0,0,0,0,0,0,0};
    guidata(hObject,handles)
    close(handles.figure1);
end


function saveclosebutton_Callback(~, ~, handles)
close(handles.figure1);

function inputamountofpointstoadd_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function typeList_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function figure1_CloseRequestFcn(hObject, ~, ~)
% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on button press in entireimage.
function entireimage_Callback(hObject, eventData, handles)
if get(hObject,'Value')
    cla(handles.displaypanel,'reset')
    hold on
    imshow(handles.originalImage)
    hskel = imshow(handles.temporarySkel);
    set(hskel, 'AlphaData', ~isnan(handles.temporarySkel))
    hold off
    set(handles.startselectionbutton,'Enable','off');
    set(handles.addpointsbutton,'Enable','off');
    set(handles.addlinebutton,'Enable','off');
    set(handles.applybutton,'Enable','on');
    handles.bigx=handles.graph.Nodes(:,:).x;
    handles.bigy=handles.graph.Nodes(:,:).y;
    hold on
    plot(handles.bigx,handles.bigy,'.b','LineWidth', 1);
    hold off
    set(handles.saveclosebutton,'Enable','on');
else
    set(handles.startselectionbutton,'Enable','on');
    cla(handles.displaypanel,'reset')
    hold on
    imshow(handles.originalImage)
    hskel = imshow(handles.temporarySkel);
    set(hskel, 'AlphaData', ~isnan(handles.temporarySkel))
    hold off

end
applybutton_Callback(handles.applybutton, eventData, handles);

% --- Executes on button press in resetpathbutton.
function resetpathbutton_Callback(hObject, ~, handles)
delete(handles.borderpoints) 
try 
    delete(handles.innerpoints)
catch
end

cla(handles.displaypanel,'reset')
hold on
imshow(handles.originalImage)
hskel = imshow(handles.temporarySkel);
set(hskel, 'AlphaData', ~isnan(handles.temporarySkel))
hold off

handles.nodenumber = [];
handles.finalX = [];
handles.finalY = [];
handles.flag = 1;
set(handles.startselectionbutton,'Enable','on');
guidata(hObject,handles)


% --- Executes when selected object is changed in pointsgroup.
function pointsgroup_SelectionChangedFcn(hObject, ~, handles)
% Make a selection on whether the new point to be added will be a border
% point or an inner point. The differente relies in how we deal with the
% point later. Are we ignoring what comes after or are we just jumping
% ahead?
guidata(hObject,handles)


% --- Executes on slider movement.
function radiusSlider_Callback(hObject, ~, handles)
temp = get(hObject,'Value');
if handles.radius ~= temp
    set(handles.radiusText,'String',temp)
    handles.radius = temp;
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function radiusSlider_CreateFcn(hObject, ~, handles)
set(hObject,'Max',255,'Min',0,'Value',10)
set(hObject,'SliderStep',[1/255 1/255])
handles.radius=10;

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
