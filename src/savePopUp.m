function varargout = savePopUp(varargin)
% SAVEPOPUP MATLAB code for savePopUp.fig
% Prompts the user to save or cancel.    

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @savePopUp_OpeningFcn, ...
                   'gui_OutputFcn',  @savePopUp_OutputFcn, ...
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

% --- Executes just before savePopUp is made visible.
function savePopUp_OpeningFcn(hObject, ~, handles, varargin)
% Update handles structure
handles.output = hObject;
guidata(hObject, handles);
% UIWAIT makes savePopUp wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = savePopUp_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
delete(handles.figure1);

function save_Callback(hObject, ~, handles)
handles.output = {'Yes'};

if get(handles.generalCurvature,'Value') == true
    handles.output{end+1} = 'General Curvature';
end
if get(handles.individCurvature,'Value') == true
    handles.output{end+1} = 'Individual Curvature';
end
if get(handles.axonLength,'Value') == true
    handles.output{end+1} = 'Axonal Length';
end
if get(handles.straightSegLengths,'Value') == true
    handles.output{end+1} = 'Straight Segment Lengths';
end
if get(handles.straightSegPos,'Value') == true
    handles.output{end+1} = 'Straight Segment Positions';
end
if get(handles.eccentricity,'Value') == true
    handles.output{end+1} = 'Eccentricity';
end
if get(handles.MDI,'Value') == true
    handles.output{end+1} = 'MT Disorganisation Index';
end
if get(handles.SI,'Value') == true
    handles.output{end+1} = 'Straightness Index';
end
if get(handles.mtArea,'Value') == true
    handles.output{end+1} = 'MT Area';
end
if get(handles.swellingArea,'Value') == true
    handles.output{end+1} = 'Swelling Area';
end
if get(handles.density,'Value') == true
    handles.output{end+1} = 'Density';
end

guidata(hObject,handles)
uiresume(handles.figure1);




function cancel_Callback(hObject, ~, handles) %#ok<*DEFNU>
handles.output = {'No','[]'};

guidata(hObject,handles)
uiresume(handles.figure1);

function generalCurvature_Callback(~, ~, ~)
function individCurvature_Callback(~, ~, ~)
function axonLength_Callback(~, ~, ~)
function straightSegLengths_Callback(~, ~, ~)
function MDI_Callback(~, ~, ~)
function SI_Callback(~, ~, ~)
function straightSegPos_Callback(~, ~, ~)
function eccentricity_Callback(~, ~, ~)
function mtArea_Callback(~, ~, ~)
function swellingArea_Callback(~, ~, ~)
function density_Callback(~, ~, ~)

function figure1_CloseRequestFcn(hObject, ~, ~)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% --- Executes on button press in selectAll.
function selectAll_Callback(hObject, ~, handles)
% hObject    handle to selectAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selectAll
if get(hObject,'Value')
   set(handles.generalCurvature,'Value',1);
   set(handles.individCurvature,'Value',1);
   set(handles.axonLength,'Value', 1);
   set(handles.straightSegLengths,'Value',1);
   set(handles.straightSegPos,'Value', 1);
   set(handles.eccentricity,'Value', 1);
   set(handles.MDI,'Value',1);
   set(handles.SI,'Value', 1);
   set(handles.mtArea,'Value', 1);
   set(handles.swellingArea,'Value',1);
   set(handles.density,'Value', 1);
end
guidata(hObject,handles)
