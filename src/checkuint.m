function checkuint(hObject, handles)
% Verifies if the image is of type uint and converts it if not, i.e.,
% normalizes pixels to be in the range 0 - 255.


if max(max(max(handles.image)))>255 && get(handles.imagenormalisecheck,'Value')
    handles.image = uint8(handles.image / 256);
end
guidata(hObject, handles)
end