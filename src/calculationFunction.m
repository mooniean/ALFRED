%% Calculation Function
function [handles] = calculationFunction(handles)
%Calculate every region and save to a file
global DEBUG;

saveAnswer = savePopUp;
finalName = ['results/' num2str(yyyymmdd(datetime)) '-' num2str(hour(datetime)) num2str(minute(datetime)) num2str(round(second(datetime))) '_' handles.fileName{handles.currentImage} ];

if string(saveAnswer(1))=="Yes"
    
    
    % Dividing types of ROIs between Axons and Disorganisation
    typeROI = 'Axon';
    [~,indexAxon] = findValue(handles.finalROIs,typeROI);
    finalLengths = zeros(numel(indexAxon),1);
    
    typeROI = 'MT disorganisation';
    [~,indexMT] = findValue(handles.finalROIs,typeROI);
    
    % Variable initiation
    finalAreasDisorg = zeros(numel(indexMT),1);
    finalAreaMT =  zeros(numel(indexMT),1);
    finalDensity = zeros(numel(indexMT),1);
    finalEccentricity = zeros(numel(indexMT),1);
    finalGeneralCurvatures = struct([]);
    finalPointCurvatures = struct([]);
    finalArcLength = struct([]);
    
    % Retroadding the conversionFactor to all of the ROIs before
    i = 1;
    conversionFactor = handles.finalROIs(i).conversionFactor;
    while isempty(conversionFactor)
        conversionFactor = handles.finalROIs(i).conversionFactor;
        i=i+1;
    end
    
    %Calculations on axons:
    if indexAxon > 0
        f = waitbar(0,'1','Name','Calculating axon results',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        for i = numel(indexAxon) % Iterate on all the ROIs that are axons
            if getappdata(f,'canceling')
                break
            end
            waitbar(i/numel(indexAxon),f,sprintf('ROI %1d of %1d',i,numel(indexAxon)))
            finalLengths(i) = handles.finalROIs(indexAxon(i)).length*handles.conversionFactor;
            handles.finalROIs(indexAxon(i)).unitLength = handles.finalROIs(indexAxon(i)).length*handles.conversionFactor;
        end
        delete(f)
    else
        finalLengths = 0;
    end
    
    checkPrint(handles.figure1,handles)
    finalName = [finalName '.mat'];
    tempROIs = handles.finalROIs;
    save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
    
    % Add straightline distances/sizes to save also as another variable,
    % i'm not sure how but we don't need it now!
    
    %Calculations on MT disorganisation:
    % NOTE: curvatures are calculated using two methods: the smoothing filter
    % method (dependent on resolution but good for fine grained
    % solutions), and the Fourier method, good for course grained
    % solutions and general curvature of the image.
    imagesWithDisorganisation = [];
    
    if indexMT>0
        %             totalCurvatures=[];
        f = waitbar(0,'1','Name','Calculating MT disorganisation results',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        for i = 1:numel(indexMT)
            
            try
                waitbar(i/numel(indexMT),f,sprintf('ROI %1d of %1d',i,numel(indexMT)))
                tempV = handles.finalROIs(indexMT(i)).virtImage;
                handles.bw = im2bw(tempV,[],0.1); %#ok<IM2BW>
                handles.bw = ~handles.bw;
                %             handles.bw = imclearborder(handles.bw);
                [labeledImage, ~] = bwlabel(handles.bw);
                handles.measurements = regionprops(labeledImage,'Area','Eccentricity');
                handles.table = sortrows(struct2table(handles.measurements),'Area','descend');
                
                if isempty(intersect(imagesWithDisorganisation,handles.finalROIs(indexMT(i)).imageName))
                    imagesWithDisorganisation = [imagesWithDisorganisation ; handles.finalROIs(indexMT(i)).imageName]; %#ok<AGROW>
                end
                
                %Calculating average values from both area and eccentricity
                sumEcc=0;
                sumArea = 0;
                tot = 0;
                for j = 1:height(handles.table)
                    if handles.table{j,1}>5
                        sumArea = sumArea + handles.table{j,1};
                        sumEcc = sumEcc + handles.table{j,2};
                        tot = tot + 1;
                    end
                end
                if ~isempty(handles.finalROIs(indexMT(i)).conversionFactor) % O Nuno editou isto, não tenho a certeza que esteja a fazer o que queres
                    conversionFactor = handles.finalROIs(indexMT(i)).conversionFactor;
                else
                    conversionFactor = handles.conversionFactor;
                end
                if DEBUG; disp('before calculations'); end
                
                finalAreasDisorg(i) = sumArea/tot*conversionFactor; % Average Disorganisation area - shouldn't it be the perimeter???
                
                se = strel('square',4);
                dilated = imdilate(handles.finalROIs(indexMT(i)).path, se);
                filled = imfill(dilated);
                area = sum(sum(filled));
                
                finalEccentricity(i) = sumEcc/tot*conversionFactor; % Average Eccentricity
                handles.adjusted = handles.finalROIs(indexMT(i)).virtImage>(handles.finalROIs(indexMT(i)).grayscale/255);
                finalAreaMT(i) = area*conversionFactor; %Area of microtubules in the image
                
                
                %                 handles.imageCurvature = handles.imageCurvature.*(handles.imageCurvature<0.2);
                %                 handles.imageCurvature(handles.imageCurvature==0) = [];
                
                %                 [~, ~, resultsFourier] = curvatureFourier(handles.finalROIs(indexMT(i)).regioncoordinates.x,handles.finalROIs(indexMT(i)).regioncoordinates.y);
                %                 [~, ~, resultsGauss] = curvatureGaussTuner(handles.finalROIs(indexMT(i)).regioncoordinates.x,handles.finalROIs(indexMT(i)).regioncoordinates.y);
                
                skeleton = curvaturePreProcessing(handles.finalROIs(indexMT(i)).regioncoordinates.x,handles.finalROIs(indexMT(i)).regioncoordinates.y);
                tol=0.2;
                [N,~,~, radiusSpline, dsdt, ~] = curvatureSpline(skeleton,tol);
                
                %resultsFourier = [];
                %                 resultsSpline = [];
                %finalGeneralCurvatures{i} = 1./(resultsFourier.*handles.conversionFactor);
                finalPointCurvatures{i} = 1./(radiusSpline{1}(1:1e-3:N{1}).*handles.conversionFactor);
                finalArcLength{i} =dsdt{1}(1:1e-3:N{1});
%                 assignin('base','newArcLength',dsdt{1}(1:1e-3:N{1}))
                
                if DEBUG; disp('after calculations'); end
                %                 totalCurvatures = [totalCurvatures handles.imageCurvature]; %#ok<AGROW>
                
                handles.finalROIs(indexMT(i)).areasDisorg = finalAreasDisorg(i);
                handles.finalROIs(indexMT(i)).areaMT = finalAreaMT(i);
                %                 handles.finalROIs(indexMT(i)).generalCurvatures = finalGeneralCurvatures{i};
                handles.finalROIs(indexMT(i)).individualCurvatures = finalPointCurvatures{i};
                handles.finalROIs(indexMT(i)).arcLength = finalArcLength{i};
                
                %             handles.finalROIs(indexMT(i)).density =
                handles.finalROIs(indexMT(i)).eccentricity = finalEccentricity(i);
                % ADD FIGURE WITH THE CURVATURE OR JUST SAVE THE FILE?
            catch ME
                if DEBUG; disp(['WARNING: ',ME.identifier]); pause(); end
                delete(f)
                continue
            end
        end
        if (isgraphics(f))
            delete(f)
        end
    else
        finalAreasDisorg = 0;
        finalAreaMT = 0;
        finalDensity = 0;
        finalEccentricity = 0;
    end
    checkPrint(handles.figure1,handles)
    %     finalName = [finalName '.mat'];
    tempROIs = handles.finalROIs;
    save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
    save('results/arcLength.mat','finalArcLength');
    
    
    
    f = waitbar(0,'1','Name','Calculating MDI results',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    for i = 1:size(imagesWithDisorganisation,1)
        try
            waitbar(i/size(imagesWithDisorganisation,1),f,sprintf('Image %1d of %1d',i,size(imagesWithDisorganisation,1)))
            [~,indexImages,~] = findValue(handles.finalROIs,imagesWithDisorganisation(i,:));
            [~,numAxon,~] = findValue(handles.finalROIs(indexImages),'Axon');
            [~,numMT,~] = findValue(handles.finalROIs(indexImages),'MT disorganisation');
            % This gives me a matrix with the overlap of the
            % disorganisation with the axons. Rows: axon, Columns: MTs
            ratios = zeros(numel(numAxon),numel(numMT));
            % %             'entering the cycles'
            for j = 1:numel(numAxon)
                for k = 1:numel(numMT)
                    % %                     'calculating bbox overlap'
                    ratios(j,k)=bboxOverlapRatio(handles.finalROIs(numAxon(j)).box,handles.finalROIs(numMT(k)).box,'Min');
                    if ratios(j,k)<0.05
                        ratios(j,k)=0;
                        continue
                    end
                    % %                     'checking skeleton overlap'
                    if ~skeletonOverlap(handles.finalROIs(numAxon(j)).box,handles.finalROIs(numMT(k)).box,handles.finalROIs(numAxon(j)).skel,handles.finalROIs(numMT(k)).skel) %if skeleton overlap is 0
                        ratios(j,k)=0;
                    end
                    % %                     'left skeleton overlap'
                end
            end
            %             'before'
            [~,indexes] = max(ratios); %gives me the max ratio per disorganisation (column)
            uniqueIndex = unique(indexes);
            for j = 1:numel(uniqueIndex)
                %                 'calculating mdi'
                
                disorgIndex = find(indexes==uniqueIndex(j));
                finalAreaMT = 0;
                for k = 1:numel(disorgIndex)
                    finalAreaMT = finalAreaMT+handles.finalROIs(numMT(disorgIndex(k))).areaMT; % tirar o +1 e por + disorg(index).perimeter
                end
                % MDI = perimeterMT/length(axon)
                handles.finalROIs(numAxon(j)).MDI = finalAreaMT/handles.finalROIs(numAxon(j)).length;
            end

        catch
            %             'I DIED'
            delete(f)
            continue
        end
        
    end
    if (isgraphics(f))
        delete(f)
    end
    checkPrint(handles.figure1,handles)
    %     finalName = [finalName '.mat'];
    tempROIs = handles.finalROIs;
    save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
    
    finalStraightLengths =  cell(numel(handles.finalROIs)-1,1);
    finalStraightPos = cell(numel(handles.finalROIs)-1,1);
    finalStraightnessIndex = cell(numel(handles.finalROIs)-1,1);
    %         handles.finalROIs(handles.imageIndex).straightlineDistances
    for i = 2:numel(handles.finalROIs)
        try
            if ~isempty(handles.finalROIs(i).conversionFactor) % O Nuno editou isto, não tenho a certeza que esteja a fazer o que queres
                conversionFactor = handles.finalROIs(i).conversionFactor;
            else
                conversionFactor = handles.conversionFactor;
            end
            finalStraightLengths{i-1} = handles.finalROIs(i).straightlineDistances.*conversionFactor;
            finalStraightPos{i-1} = handles.finalROIs(i).straightlineCoord;
            
            straightLines = finalROIs(i).straightlines;
            sizeStraightLines = length(straightLines);
            distances=[];
            segmentNumber = 1;
            distanceIndex = 1;
            k=1;
            while k<=sizeStraightLines
                try
                    if (straightLines(k).rho == straightLines(k+1).rho) && (straightLines(k).theta == straightLines(k+1).theta)
                        distances = [distances finalROIs(i).straightlineDistances(k)+finalROIs(i).straightlineDistances(k+1)];
                        k=k+1;
                    else
                        distances = [distances finalROIs(i).straightlineDistances(k)];
                    end
                catch
                    distances = [distances finalROIs(i).straightlineDistances(k)];
                end
                k=k+1;
            end
            
            if strcmp(handles.finalROIs(i).type,'MT disorganisation') % THESE ARE INDEXES, i don't need the conversion factor
                %                 finalStraightnessIndex{i-1} = sum(handles.finalROIs(i).straightlineDistances)./handles.finalROIs(i).areaMT;
                
                
                straightnessRatio = sum(distances)/sum(sum(finalROIs(indexMT(i)).skel))*100;
                %                 sum(sum(finalROIs(indexMT(i)).skel))
                % straightnessRatio = sum(tempROIs(roiIndex).straightlineDistances)/sum(sum(tempROIs(roiIndex).skel))*100
            else
                %                 finalStraightnessIndex{i-1} = sum(handles.finalROIs(i).straightlineDistances)./handles.finalROIs(i).length;
                totallength =tempROIs(i).length;
                straightnessRatio = sum(distances)/totallength*100;
                
            end
            finalStraightnessIndex{i-1} = straightnessRatio;
            numberSegments = length(distances);
        catch
            continue
        end
    end
    
    %     csvname = [finalName '.csv'];
    
    % cellROIs = struct2cell(handles.finalROIs);
    i=2; % Index for what the user wants
    j=1; % index for the finalResults struct
    finalResults = struct([]);
    f = waitbar(0,'1','Name','Saving',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    checkPrint(handles.figure1,handles)
    while i<=size(saveAnswer,2)
        finalResults(j).name = saveAnswer(i);
        switch string(saveAnswer(i))
            case 'General Curvature'
                finalResults(j).value = finalGeneralCurvatures;
            case 'Individual Curvature'
                finalResults(j).value = finalPointCurvatures;
            case 'Axonal Length'
                finalResults(j).value = finalLengths;
            case 'Straight Segment Lengths'
                finalResults(j).value = finalStraightLengths;
            case 'Straight Segment Positions'
                finalResults(j).value = finalStraightPos;
            case 'Eccentricity'
                finalResults(j).value = finalEccentricity;
                %                 case 'MT Disorganisation Index'
                %                     finalResults(j).value = ;
            case 'Straightness Index'
                finalResults(j).value = finalStraightnessIndex;
            case 'MT Area'
                finalResults(j).value = finalAreaMT;
            case 'Swelling Area'
                finalResults(j).value = finalAreasDisorg;
            case 'Density'
                finalResults(j).value = finalDensity;
        end
        i = i+1;
        j = j+1;
        
    end
    delete(f)
    checkPrint(handles.figure1,handles)
    %         finalResults = struct('AxonalLengths',finalLengths,'AreaOfDisorganisation',finalAreasDisorg,'AreaMicrotubules',finalAreaMT,'Density',finalDensity,'Eccentricity',finalEccentricity,'Curvatures',finalGeneralCurvatures);
    %         c=struct2cell(finalResults);
    %         xlsname = [finalName '.csv'];
    %                 xlswrite(xlsname, vertcat(c{:}));
    %or if you want the field names:
    %         xlswrite(xlsname, [fieldnames(finalResults), c]);
    
    %     struct2csv(finalResults,csvname);
    
    msgbox('The requested files have been saved.', 'Saved');
    
    set(handles.undobutton,'Enable','off');
    %catch ex
    %    disp(ex)
    %    warndlg('No conversion factor was selected. Use Microscope Specs')
    %end
    
end
% finalName = [finalName '.mat'];
tempROIs = handles.finalROIs;
save(finalName,'tempROIs'); % Saves everything into a file that can be opened later
end
% guidata(hObject,handles)