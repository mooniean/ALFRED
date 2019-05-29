function varargout = ALFRED_microscopespecs(varargin)
%         Before ALFRED calculates all of the values mentioned, it needs a value for the conversion from pixels to another unit of length. 
%         There are three options that a user can choose from in order to select the conversion.
% 
%         First, the user can select the metadata. Certain types of files, such as images froma STED microscope, will include metadata with each image. 
%         This usually includes the pixel/unit conversion needed. As such, with a simple checkbox, that value will be accessed and the user doesn't need to do anything further.
% 
% 
%         Second, if the images don't provide the correct needed conversion, then the user can select from a predefined scaling factor spreadsheet. 
%         Currently, this can only be altered manually previous to running of the program, but it's an easy solution. 
%         On the folder ndocuments, there is an Excel spreadsheet named "scalingfactors.xlsx". Inside, there is a table with the conversion values for the selected microscope.
%         When the user selects the microscope, the Objective will automatically update with the available options from the sheet, as well as Optovar and Binning options after each has been sequentially selected. 
%         Finally, ALFRED will get the designated value from the provided table, as soon as the user presses Apply.
% 
%         Third, the user can insert manually the value for conversion, where the default unit will be micrometers (µm), but can be changed with the drop-down menu. 
%         If the user wishes to obtain the values in pixels, just insert manually 1 and no conversion will be calculated.
% 
%         After each selection, the user can simply press the button "Apply" and continue the analysis. On the other hand, if the user doesn't want any changes, the "Cancel" button
%         can be pressed.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ALFRED_microscopespecs_OpeningFcn, ...
    'gui_OutputFcn',  @ALFRED_microscopespecs_OutputFcn, ...
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

% --- Executes just before ALFRED_microscopespecs is made visible.
function ALFRED_microscopespecs_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;
set(handles.figure1, 'Name', 'ALFRED - Microscope Specifications');
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jframe=get(handles.figure1,'javaframe');
jIcon=javax.swing.ImageIcon('icons\alfred_50x50.png');
jframe.setFigureIcon(jIcon); % Setting icon of the program

handles.microscope = 'Microscope';
handles.metadata = 'Off';
handles.optovar = 'Optovar';
handles.binning = 'Binning';
handles.objective = 'Objective';
handles.unit = 'Unit';
handles.numvalue = 0;
handles.change = 0;
handles.result ={handles.change, handles.microscope,handles.objective,handles.binning,handles.optovar,num2str(handles.numvalue),handles.unit,handles.metadata};

set(handles.metadatacheck,'Value',0);
guidata(hObject, handles);
% UIWAIT makes ALFRED_microscopespecs wait for user response (see UIRESUME)
uiwait(handles.figure1);

function varargout = ALFRED_microscopespecs_OutputFcn(~, ~, handles)
varargout{1} = handles.result;
delete(handles.figure1);

%Functions executed when things are selected from the lists
function microscopeChoice_Callback(hObject, ~, handles) %#ok<*DEFNU>
contents = cellstr(get(hObject,'String'));
handles.microscope = contents{get(hObject,'Value')};
switch handles.microscope
    case 'AxioCam 506m'
        listboxItems = {'Objective','10','40 air', '40 sil oil', '60','100'};
    case 'AxioCam MR5m'
        listboxItems = {'Objective','10','40 air', '60','100'};
    case 'Microscope'
        listboxItems = {'Objective'};
        set(handles.optovarChoice,'Value', 1); %The position of the selected value has to change to 1, otherwise it can become out of bounds
        set(handles.binningChoice,'Value', 1);
        set(handles.objectiveChoice,'Value', 1);
        set(handles.optovarChoice,'String', {'Optovar'} );
        set(handles.binningChoice,'String', {'Binning'} );
        handles.optovar = 'Optovar';
        handles.binning = 'Binning';
        handles.objective = 'Objective';
        handles.change = 0;
end
set(handles.objectiveChoice,'String', listboxItems );
guidata(hObject,handles)
function objectiveChoice_Callback(hObject,~, handles)
contents = cellstr(get(hObject,'String'));
handles.objective = contents{get(hObject,'Value')};
% if strcmp(handles.objective, '40 sil oil')
%     if strcmp(handles.optovar, '1')
%         listboxItems = {'Binning','1x1'};
%     else
%         listboxItems = {'Binning','1x1','2x2'};
%     end
% else
%     listboxItems = {'Binning','1x1','2x2', '3x3', '4x4','5x5'};
% end
% set(handles.binningChoice,'String', listboxItems );
listboxItems = {'Optovar','1','1.25', '1.6', '2'};
set(handles.optovarChoice,'String', listboxItems );
guidata(hObject,handles)

function optovarChoice_Callback(hObject, ~, handles)
contents = cellstr(get(hObject,'String'));
handles.optovar = contents{get(hObject,'Value')};
if strcmp(handles.objective, '40 sil oil')
    if strcmp(handles.optovar, '1')
        listboxItems = {'Binning','1x1'};
    else
        listboxItems = {'Binning','1x1','2x2'};
    end
else
    listboxItems = {'Binning','1x1','2x2', '3x3', '4x4','5x5'};
end
set(handles.binningChoice,'String', listboxItems );
guidata(hObject,handles)
function binningChoice_Callback(hObject, ~, handles)
contents = cellstr(get(hObject,'String'));
handles.binning = contents{get(hObject,'Value')};
guidata(hObject,handles)


function value_Callback(hObject, ~, handles)
handles.numvalue = str2double(get(hObject,'String'));
if isnan(handles.numvalue)
    handles.numvalue=0;
    set(handles.value,'String','');
else
    if get(handles.unitList,'Value') == 1
        set(handles.unitList,'Value',3);
        handles.unit = 'um';
    end
end
guidata(hObject, handles);
function unitList_Callback(hObject, ~, handles)
contents = cellstr(get(hObject,'String'));
handles.unit = contents{get(hObject,'Value')};
guidata(hObject,handles)


%Functions executed on button press.
function applybutton_Callback(hObject, ~, handles)
if get(handles.microscopeChoice,'Value') ~= 1 && get(handles.optovarChoice,'Value') ~= 1 && get(handles.binningChoice,'Value') ~= 1 && get(handles.objectiveChoice,'Value') ~= 1
    handles.change = 1;
elseif get(handles.unitList,'Value') ~= 1 && handles.numvalue > 0
    handles.change = 2;
end
handles.result = {handles.change, handles.microscope,handles.objective,handles.binning,handles.optovar,num2str(handles.numvalue),handles.unit,handles.metadata};

guidata(hObject, handles)
close(handles.figure1)

function cancelbutton_Callback(hObject, ~, handles)
handles.microscope = 'Microscope';
handles.optovar = 'Optovar';
handles.binning = 'Binning';
handles.objective = 'Objective';
handles.value = 0;
handles.unit = 'Unit';
handles.change = 0;
handles.result = {handles.change, handles.microscope,handles.objective,handles.binning,handles.optovar,num2str(handles.numvalue),handles.unit,handles.metadata};

guidata(hObject,handles)

close(handles.figure1);
function microscopeChoice_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function objectiveChoice_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function optovarChoice_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function binningChoice_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function value_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function unitList_CreateFcn(hObject, ~, ~)
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
function metadatacheck_Callback(hObject, ~, handles)
if get(hObject,'Value')
    set(handles.microscopeChoice,'Value', 1);
    set(handles.microscopeChoice,'String',{'Microscope'});
    set(handles.objectiveChoice,'String',{'Objective'});
    set(handles.optovarChoice,'String', {'Optovar'} );
    set(handles.binningChoice,'String', {'Binning'} );
    set(handles.value, 'Enable', 'off');
    set(handles.optovarChoice,'Value', 1); %The position of the selected value has to change to 1, otherwise it can become out of bounds
    set(handles.binningChoice,'Value', 1);
    set(handles.objectiveChoice,'Value', 1);
    handles.metadata = 'On';
    handles.microscope = 'Microscope';
    handles.optovar = 'Optovar';
    handles.binning = 'Binning';
    handles.objective = 'Objective';
    handles.change = 3;
else
    set(handles.microscopeChoice,'Value', 1);
    set(handles.microscopeChoice,'String',{'Microscope','AxioCam 506m','AxioCam MR5m'});
    set(handles.value, 'Enable', 'on');
    handles.metadata = 'Off';
end
guidata(hObject,handles)