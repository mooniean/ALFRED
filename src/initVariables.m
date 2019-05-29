function initVariables(hObject, handles)
% Initialisation of the variables in the le.


% Restore all the buttons to initial state
handles.region = handles.image;
handles.selectedregion = handles.image;

set(handles.normalize,'Value',1);
set(handles.checkAvDisArea,'Value',1);
set(handles.checkTotalDisArea,'Value',1);
set(handles.checkRegions,'Value',1);
set(handles.checkJunctions,'Value',1);
set(handles.checkEcc,'Value',1);
set(handles.checkTotalArea,'Value',1);
% set(handles.checkAvCurv,'Value',1);
% set(handles.checkRadius,'Value',1);
% set(handles.checkStraight,'Value',1);
set(handles.blackbackground,'Value',1);

set(handles.meanarea,'String','');
set(handles.totalarea,'String','');
set(handles.numregions,'String','');
set(handles.numjunctions,'String','');
set(handles.avEccentricity,'String','');
set(handles.totWhiteArea,'String','');
% set(handles.avCurvature,'String','');
% set(handles.straightOut,'String','');
% set(handles.avRadius,'String','');
set(handles.textdisplay,'String',['Image Number ', num2str(handles.current) ,' ', num2str(size(handles.image,1)), 'x', num2str(size(handles.image,2)) ]);

set(handles.savebutton,'Enable','off');
set(handles.calculateEnd,'Enable','off');
set(handles.activatemicrotubules,'Enable','off');
set(handles.regionbutton,'Enable','on');
set(handles.restorebutton,'Enable','off');
set(handles.applycontrast,'Enable','off');
set(handles.grayscale,'Enable','off');

handles.scaleValue = 3;
handles.senseValue = 1;
handles.conMinValue = 0.01;
handles.conMaxValue = 0.99;
handles.grayScaleValue = 0;
set(handles.contrastmin,'Value',handles.conMinValue);
set(handles.valConMin,'String',get(handles.contrastmin,'Value'))
set(handles.contrastmax,'Value',handles.conMaxValue);
set(handles.valConMax,'String',get(handles.contrastmax,'Value'))
set(handles.sensitivity,'Value',handles.senseValue);
set(handles.valSensitivity,'String',get(handles.sensitivity,'Value'))
set(handles.scale,'Value',handles.scaleValue);
set(handles.valScale,'String',get(handles.scale,'Value'))
set(handles.grayscaleSlider,'Value',handles.grayScaleValue);
set(handles.grayValue,'String',get(handles.grayscaleSlider,'Value'))

imshow(handles.region)
guidata(hObject, handles);

end