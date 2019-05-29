function [handles] = virtualising(handles)
% Function that will execute the virtualisation of the received binary
% image. Updates handles,


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
set(handles.regionbutton,'Enable','off');
end