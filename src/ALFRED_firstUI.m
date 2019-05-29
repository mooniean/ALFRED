function varargout = ALFRED_firstUI(varargin)
% The main panel opening when running ALFRED is the Loading Panel.
% This is the first window to popup when the software is opened and allows the user to load any type of image format, from simple jpg files to complex bioformat types such as STED stacks.
% By pressing the button "Load", the user can choose from which folder to load the images.
% A small popup will keep track of how many images have been loaded and how many are left.
% Secondly, in the case of a single image being uploaded and it being a panel, it can be divided into rectangles, simply by selecting the number of rows and columns and pressing "Update".
% The insertion of numbers can either be done by clicking the text box or by pressing tab twice.
% The white border checkbox is for when the images are separated by a defined border. "Restore" will undo the division and the image will show
% complete again.
% Each image or division will be numbered, with a number that will also be saved in the final file.
% If more than 30 images are loaded, or if a panel has more than 30 divisions, then a new page is created and the user can navigate between them by using the "Previous Page" and "Next Page" buttons.
% To open the next panel, the user simply needs to click the image intended. A new pop-up should open, the main image panel.
% At the bottom left, there is a help button in case the user needs any information regarding the panel.
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ALFRED_firstUI_OpeningFcn, ...
    'gui_OutputFcn',  @ALFRED_firstUI_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before ALFRED_firstUI is made visible.
function ALFRED_firstUI_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;

set(handles.figure1, 'Name', 'ALFRED Panel');
set(handles.updatedisplay,'Enable','off');
set(handles.back,'Enable','off');

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1,'javaframe');
jIcon=javax.swing.ImageIcon('icons/alfred_50x50.png');
jframe.setFigureIcon(jIcon);

handles.notdisplayable = imread('icons/notdisplayable.png');
handles.currentPage = 1;
handles.totalPages = 1;
handles.whiteborder = 0;
handles.widthIncr = 1;
handles.heightIncr = 1;
handles.zmax = 1;
handles.z = 1;
handles.closeWindow = 0;
handles.SubPlotAxesObj={};

set(handles.updatedisplay,'Enable','off');
set(handles.nextpagebutton,'Enable','off');
set(handles.previouspagebutton,'Enable','off');


guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ALFRED_firstUI_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

function helpbutton_Callback(~, ~, ~)
msgbox({'Welcome to A.L.F.R.E.D.' 'Load the images that you intend to analyse. If you are working with panels, open one and use the division option to separate the squares.' 'Otherwise, load all the images at once and press the first one to start.' 'If you have any doubts, please contact the developing team (see README file).'},'Help!');

function rowinput_Callback(hObject, ~, handles)
handles.heightIncr = str2double(get(hObject,'String'));
if isnan(handles.heightIncr)
    handles.heightIncr=1;
    set(handles.rowinput,'String','1');
end
guidata(hObject, handles);

function colinput_Callback(hObject, ~, handles)
handles.widthIncr = str2double(get(hObject,'String'));
if isnan(handles.widthIncr)
    handles.widthIncr=1;
    set(handles.colinput,'String','1');
end
guidata(hObject, handles);

function whiteborder_Callback(hObject, ~, handles)
handles.whiteborder = get(hObject,'Value');
guidata(hObject, handles);

function rowinput_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function colinput_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updatedisplay_Callback(hObject, ~, handles)
incr = 0;
if handles.whiteborder
    incr = 28;
end
matrix = {};
oldNames = char(handles.FileNames);
posx=0;posy=0;
width= size(handles.imagegroup{1},2)/handles.widthIncr - incr;
height= size(handles.imagegroup{1},1)/handles.heightIncr - incr;
zmax = handles.widthIncr*handles.heightIncr;
newNames= {};
for z=1:zmax
    matrix{z}=imcrop(handles.imagegroup{1},[posx posy width height]); %#ok<AGROW>
    newNames{z}=[oldNames '_' num2str(z)]; %#ok<AGROW>
    posx = posx + width + incr;
    if (mod(z,handles.widthIncr)==0) % write {handles.widthIncr+1,2*handles.widthIncr+1,3*handles.widthIncr+1,4*handles.widthIncr+1}
        posx = 0;
        posy = posy + height+incr;
    end
end
% delete(handles.display)
if size(handles.SubPlotAxesObj)>0
    for i = 1:length(handles.SubPlotAxesObj)
        delete(handles.SubPlotAxesObj{i});
    end
end
handles.SubPlotAxesObj={};
for z=1:zmax
    handles.SubPlotAxesObj{z}=subplot(handles.heightIncr,handles.widthIncr,z,'Parent',handles.displaypanel);
    p=imagesc(matrix{z});
    set(p,'ButtonDownFcn',{@copy2NewFigure,matrix,zmax,z, newNames,handles.metadata}); axis off;
end

set(handles.back,'Enable','on');

guidata(hObject, handles);

function copy2NewFigure(~,~,matrix,zmax,z,name,metadata) % (src, evt, ~)
ALFRED_mainpanel(matrix,zmax,z,name,metadata);

function back_Callback(hObject, ~, handles)
set(handles.back,'Enable','off');
%set(handles.display,'Enable','on');
for z=1:length(handles.SubPlotAxesObj)
    delete(handles.SubPlotAxesObj{z})
end
handles.display=axes('Parent',handles.displaypanel);
p=imagesc(handles.imagegroup{1});
set(p,'ButtonDownFcn',{@copy2NewFigure,handles.imagegroup{1},1,1,handles.FileNames{1}}); axis off;
guidata(hObject, handles);

function loadbutton_Callback(hObject, ~, handles) %#ok<*DEFNU>
warning('off','all')
pass = true;
handles.FileNames = {1};
while pass
    [ImageNames, PathNames] = uigetfile('*.*','Pick the images','MultiSelect','on');
    if ~iscellstr(ImageNames)
        handles.nfiles = 1;
        handles.FileNames{1} = ImageNames;
        set(handles.updatedisplay,'Enable','on');
    else
        handles.nfiles = size(ImageNames,2);
        handles.FileNames = ImageNames;
    end
    
    if handles.nfiles < 1
        errordlg('Please select at least one image.','Insufficient Number of Files');
        continue;
    end
    pass = false;
end

f = waitbar(0,'Loading...','Name','Loading images...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

handles.datagroup = {};
handles.metadata = {};
handles.imagegroup = {};
step = 1;
for i=1:handles.nfiles
    file = fullfile(PathNames, handles.FileNames{i});
    if getappdata(f,'canceling')
        break
    end
    try
        
        data = bfopen(file);
        waitbar(i/handles.nfiles,f,sprintf('Image %1d of %1d',i,handles.nfiles))
        
        handles.metadata{i} = data(:,4);
        handles.datagroup{i} = data(:,1);
        currentDataGroup = handles.datagroup{i};
        for j = 1:size(handles.datagroup{i})
            image = concatCells(currentDataGroup{j});
            switch size(image,3)
                case 1
                    image = cat(3,image,image,image);
                case 2
                    image = cat(3,image,image(:,:,1));
            end
            handles.imagegroup{step} = image;
            step = step + 1;
        end
    catch
        image = imread(file);
        handles.metadata{1} = 0;
        handles.datagroup{1} = 0;
        handles.imagegroup{step} = image;
        step = step + 1;
    end
    %     image = imread(file);
    %     if size(image,3)<3 %careful with two channels
    %         image = cat(3,image,image,image);
    %     end
    %     handles.imagegroup{i} = image;
    
end
delete(f)

handles.nperpage = 30;

if handles.nfiles > handles.nperpage
    handles.totalPages = ceil(handles.nfiles/handles.nperpage);
    set(handles.nextpagebutton,'Enable','on');
elseif handles.nfiles == 1 && size(handles.datagroup,1)>1
    handles.totalPages = ceil(size(handles.datagroup,1)/handles.nperpage);
    set(handles.nextpagebutton,'Enable','on');
else
    handles.nperpage = handles.nfiles;
end


if size(handles.SubPlotAxesObj)>0
    for i = 1:length(handles.SubPlotAxesObj)
        delete(handles.SubPlotAxesObj{i})
    end
end
handles.display=axes('Parent',handles.displaypanel);
if handles.nfiles == 1 && (length(handles.imagegroup) == 1)
    p=imagesc(handles.imagegroup{1});
    set(p,'ButtonDownFcn',{@copy2NewFigure,handles.imagegroup,1,1,handles.FileNames{1},handles.metadata{1}}); axis off;
    set(handles.filenamedisplay, 'String', handles.FileNames{1});
    set(handles.updatedisplay,'Enable','on');
else
    handles.SubPlotAxesObj={};
    if size(handles.datagroup,1)>1
        ImageNames = handles.FileNames;
        for j = 1:size(handles.datagroup,1)
            handles.FileNames{j}=ImageNames;
        end
    end
    for z=1:handles.nperpage
        subPlotRows=floor(sqrt(handles.nperpage));
        subPlotCols=ceil(handles.nperpage/subPlotRows);
        handles.SubPlotAxesObj{z}=subplot(subPlotRows,subPlotCols,z,'Parent',handles.displaypanel);
        if size(handles.imagegroup{z},3) <= 3
            p=imagesc(handles.imagegroup{z});
        else
            p=imshow(handles.imagegroup{z}(:,:,1));
            text(0,100,'Channel 1','Color','white');
        end
        set(p,'ButtonDownFcn',{@copy2NewFigure,handles.imagegroup,handles.nfiles,z, handles.FileNames,handles.metadata}); axis off; title(num2str(z));
    end
end

guidata(hObject, handles);

function nextpagebutton_Callback(hObject, ~, handles)
step = handles.currentPage * handles.nperpage;
beginning = step + 1;
ending = beginning + (handles.nperpage-1);
if handles.nfiles>1
    if ending > handles.nfiles
        ending = handles.nfiles;
    end
else
    if ending > size(handles.datagroup{1},1)
        ending = size(handles.datagroup{1},1);
    end
end

for i = 1:length(handles.SubPlotAxesObj)
    delete(handles.SubPlotAxesObj{i})
end

handles.SubPlotAxesObj={};

for z = beginning:ending
    plotIndex = z-step;
    subPlotRows=floor(sqrt(handles.nperpage));
    subPlotCols=ceil(handles.nperpage/subPlotRows);
    handles.SubPlotAxesObj{plotIndex}=subplot(subPlotRows,subPlotCols,plotIndex,'Parent',handles.displaypanel);
    if size(handles.imagegroup{z},3) <= 3
        p=imagesc(handles.imagegroup{z});
    else
        p=imshow(handles.imagegroup{z}(:,:,1));
        text(0,100,'Channel 1','Color','white');
    end
    set(p,'ButtonDownFcn',{@copy2NewFigure,handles.imagegroup,handles.nfiles,z, handles.FileNames,handles.metadata}); axis off; title(num2str(z));
end

set(handles.previouspagebutton,'Enable','on');
handles.currentPage = handles.currentPage + 1;
if handles.currentPage == handles.totalPages
    set(handles.nextpagebutton,'Enable','off');
end
guidata(hObject,handles);

function previouspagebutton_Callback(hObject, ~, handles)
step = handles.currentPage * handles.nperpage;
ending = step;
beginning = ending - (handles.nperpage-1);

for i = 1:length(handles.SubPlotAxesObj)
    delete(handles.SubPlotAxesObj{i})
end
plotIndex = 1;
handles.SubPlotAxesObj={};
for z = beginning:ending
    subPlotRows=floor(sqrt(handles.nperpage));
    subPlotCols=ceil(handles.nperpage/subPlotRows);
    handles.SubPlotAxesObj{plotIndex}=subplot(subPlotRows,subPlotCols,plotIndex,'Parent',handles.displaypanel);
    if size(handles.imagegroup{z},3) <= 3
        p=imagesc(handles.imagegroup{z});
    else
        p=imshow(handles.imagegroup{z}(:,:,1));
        text(0,100,'Channel 1','Color','white');
    end
    set(p,'ButtonDownFcn',{@copy2NewFigure,handles.imagegroup,handles.nfiles,z, handles.FileNames,handles.metadata}); axis off; title(num2str(z));
    plotIndex = plotIndex + 1;
end

set(handles.nextpagebutton,'Enable','on');
handles.currentPage = handles.currentPage - 1;
if handles.currentPage == 1
    set(handles.previouspagebutton,'Enable','off');
end
guidata(hObject,handles);
