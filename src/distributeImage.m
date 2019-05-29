function imgcut = distributeImage(hObject,temp,handles)
% Distributes the region selected per all three channels.

imgcut(:,:,1) = temp(handles.b(2):handles.b(2)+handles.b(4), handles.b(1):handles.b(1)+handles.b(3),1);
imgcut(:,:,2) = temp(handles.b(2):handles.b(2)+handles.b(4), handles.b(1):handles.b(1)+handles.b(3),2);
imgcut(:,:,3) = temp(handles.b(2):handles.b(2)+handles.b(4), handles.b(1):handles.b(1)+handles.b(3),3);
guidata(hObject, handles);
end