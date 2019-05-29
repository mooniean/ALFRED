%Extra lines of code so I can clean up the different files

%% Lengthpopup


% plotImageGraph(handles.g)
% try
%     [handles.totalSize,handles.path,handles.bigx,handles.bigy] = lengthFinding(hObject,handles);
% catch EX
%     EX.identifier
%     addlinebutton_Callback(handles.addlinebutton,0,handles);
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


% [handles.totalSize,handles.path,handles.bigx,handles.bigy] = lengthFinding(hObject,handles);
% handles.finalPath=zeros(size(handles.image));
% for i = 1:length(handles.bigx)
%     handles.finalPath(handles.bigy(i),handles.bigx(i))=1;
% end
% [handles.lines,handles.longLine] = houghtest(handles.finalPath,handles.numpeaks);
% set(handles.estimatedlengthoutput,'String',num2str(handles.totalSize));
% set(handles.selectedpoints,'String',num2str(handles.totalpoints));
% hold on; plot(handles.bigx,handles.bigy,'r'); plot([handles.finalX],[handles.finalY],'g.'); hold off;
%             [handles.totalSize,handles.path,handles.bigx,handles.bigy] = lengthFinding(hObject,handles);

%         x1=handles.g.Nodes(handles.lines(i).point1,:).x;
%         x2=handles.g.Nodes(handles.lines(i).point2,:).x;
%         y1=handles.g.Nodes(handles.lines(i).point1,:).y;
%         y2=handles.g.Nodes(handles.lines(i).point2,:).y;
%     longestLineDistance = checkDistance(handles.longLine(1,1),handles.longLine(1,2),handles.longLine(2,1),handles.longLine(2,2));
%     fprintf('Filename \t Length of Axon \t Average Line \t Largest Line \n %s \t %f \t %.4f \t %.4f \n',handles.fileName{handles.currentImage},handles.ROIlength,avg,longestLineDistance);
      
    
    
    
    % handles.output = {handles.regiontype, handles.path, handles.totalSize,handles.lines,handles.longLine,handles.bigx,handles.bigy,linesDist,linesCoord};
% guidata(hObject,handles)



%%