function checkPrint(hObject,handles)
% Checks which regions have been selected and designated in an image, particularly when changing back and forth between images using the main panel.

name = handles.fileName{handles.currentImage};
[existInStruct,index] = findValue(handles.finalROIs,name);
if existInStruct == 0
    %It's the first time we're opening this image. Indexes need to be
    %updated, and we need to save this image.
    handles.imageIndex = handles.imageIndex + 1;
    handles.ROIindex = 1;
    handles.finalROIs(handles.imageIndex).imageName = handles.fileName{handles.currentImage};
    handles.finalROIs(handles.imageIndex).ROInumber =  handles.baseImageIndex;
    handles.finalROIs(handles.imageIndex).type = 'Full';
    handles.finalROIs(handles.imageIndex).conversionFactor = handles.conversionFactor;
    try
        handles.finalROIs(handles.imageIndex).metadata = handles.metadata{handles.currentImage};
    catch
        handles.finalROIs(handles.imageIndex).metadata = 0;
    end
elseif existInStruct == 1
    handles.ROIindex = str2double(handles.finalROIs(index).ROInumber)+1;
elseif existInStruct > 1
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
handles.oldImage = handles.region;
guidata(hObject,handles)

end