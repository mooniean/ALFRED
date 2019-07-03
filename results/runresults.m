

finalName = [num2str(yyyymmdd(datetime)) '-' num2str(hour(datetime)) num2str(minute(datetime)) num2str(round(second(datetime))) '_' tempROIs(1).imageName ];


typeROI = 'Axon';
[~,indexAxon] = findValue(tempROIs,typeROI);
finalLengths = zeros(numel(indexAxon),1);

typeROI = 'MT disorganisation';
[~,indexMT] = findValue(tempROIs,typeROI);
finalAreasDisorg = zeros(numel(indexMT),1);
finalAreaMT =  zeros(numel(indexMT),1);
finalDensity = zeros(numel(indexMT),1);
finalEccentricity = zeros(numel(indexMT),1);
finalGeneralCurvatures = struct([]);
finalPointCurvatures = struct([]);
i = 1;
conversionFactor = tempROIs(i).conversionFactor;
while isempty(conversionFactor)
    conversionFactor = tempROIs(i).conversionFactor;
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
        
        finalLengths(i) = tempROIs(indexAxon(i)).length*conversionFactor;
        
        tempROIs(indexAxon(i)).unitLength = tempROIs(indexAxon(i)).length*conversionFactor;
    end
    delete(f)
else
    finalLengths = 0;
end
% checkPrint(handles.figure1,handles)
finalName = [finalName '.mat'];
finalROIs = tempROIs; 
save(finalName,'finalROIs'); % Saves everything into a file that can be opened later
% Add straightline distances/sizes to save also as another variable,
% i'm not sure how but we don't need it now!

%Calculations on MT disorganisation:
% NOTE: curvatures are calculated using two methods: the Gaussian
% window method (dependent on resolution but good for fine grained
% solutions), and the Fourier method, good for course grained
% solutions and general curvature of the image.
imagesWithDisorganisation = [];
if indexMT>0
    %             totalCurvatures=[];
    f = waitbar(0,'1','Name','Calculating MT disorganisation results',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    for i = 1:numel(indexMT)
        waitbar(i/numel(indexMT),f,sprintf('ROI %1d of %1d',i,numel(indexMT)))
        tempV = tempROIs(indexMT(i)).virtImage;
        bw = im2bw(tempV,[],0.1); %#ok<IM2BW>
        bw = ~bw;
        %             handles.bw = imclearborder(handles.bw);
        [labeledImage, ~] = bwlabel(bw);
        measurements = regionprops(labeledImage,'Area','Eccentricity');
        table = sortrows(struct2table(measurements),'Area','descend');
        
        if isempty(intersect(imagesWithDisorganisation,tempROIs(indexMT(i)).imageName))
            imagesWithDisorganisation = [imagesWithDisorganisation ; tempROIs(indexMT(i)).imageName]; %#ok<AGROW>
        end
        
        %Calculating average values from both area and eccentricity
        sumEcc=0;
        sumArea = 0;
        tot = 0;
        for j = 1:height(table)
            if table{j,1}>5
                sumArea = sumArea + table{j,1};
                sumEcc = sumEcc + table{j,2};
                tot = tot + 1;
            end
        end
        if ~isempty(tempROIs(indexMT(i)).conversionFactor) % O Nuno editou isto, não tenho a certeza que esteja a fazer o que queres
            conversionFactor = tempROIs(indexMT(i)).conversionFactor;

        end
        finalAreasDisorg(i) = sumArea/tot*conversionFactor; % Average Disorganisation area - shouldn't it be the perimeter???
        finalEccentricity(i) = sumEcc/tot*conversionFactor; % Average Eccentricity
        adjusted = tempROIs(indexMT(i)).virtImage>(tempROIs(indexMT(i)).grayscale/255);
        finalAreaMT(i) = sum(sum(adjusted))*conversionFactor; %Area of microtubules in the image
        
        
%% Curvature Calculation: Uncomment this part if needed     
        
        [~, ~, radiiRot] = curvatureFourier(tempROIs(indexMT(i)).regioncoordinates.x,tempROIs(indexMT(i)).regioncoordinates.y);
        
%         [~, ~, results] = curvatureGauss(tempROIs(indexMT(i)).regioncoordinates.x,tempROIs(indexMT(i)).regioncoordinates.y);
        
        
%         finalPointCurvatures{i} = 1./(results.*conversionFactor);
        
        finalGeneralCurvatures{i} = 1./(radiiRot.*conversionFactor);

        tempROIs(indexMT(i)).generalCurvatures = finalGeneralCurvatures{i};
%         tempROIs(indexMT(i)).individualCurvatures = finalPointCurvatures{i};
%% End of Curvature Calculation
        
        tempROIs(indexMT(i)).areasDisorg = finalAreasDisorg(i); %#ok<*SAGROW>
        tempROIs(indexMT(i)).areaMT = finalAreaMT(i);
        
        %             tempROIs(indexMT(i)).density =
        tempROIs(indexMT(i)).eccentricity = finalEccentricity(i);
        % ADD FIGURE WITH THE CURVATURE OR JUST SAVE THE FILE?
    end
    delete(f)
else
    finalAreasDisorg = 0;
    finalAreaMT = 0;
    finalDensity = 0;
    finalEccentricity = 0;
end
% finalName = [finalName '.mat'];
finalROIs = tempROIs; 
save(finalName,'finalROIs'); % Saves everything into a file that can be opened later



f = waitbar(0,'1','Name','Calculating MDI results',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for i = 1:size(imagesWithDisorganisation,1)
    waitbar(i/size(imagesWithDisorganisation,1),f,sprintf('Image %1d of %1d',i,size(imagesWithDisorganisation,1)))
    [~,indexImages] = findValue(tempROIs,imagesWithDisorganisation(i,:));
    [~,numAxon] = findValue(tempROIs(indexImages),'Axon');
    [~,numMT] = findValue(tempROIs(indexImages),'MT disorganisation');
    % This gives me a matrix with the overlap of the
    % disorganisation with the axons. Rows: axon, Columns: MTs
    ratios = zeros(numel(numAxon),numel(numMT));
    for j = 1:numel(numAxon)
        
        for k = 1:numel(numMT)
            ratios(j,k)=bboxOverlapRatio(tempROIs(numAxon(j)).box,tempROIs(numMT(k)).box,'Min');
            if ratios(j,k)<0.05
                ratios(j,k)=0;
                continue
            end
            if ~skeletonOverlap(tempROIs(numAxon(j)).box,tempROIs(numMT(k)).box,tempROIs(numAxon(j)).skel,tempROIs(numMT(k)).skel) %if skeleton overlap is 0
                ratios(j,k)=0;
            end
        end
    end
    [~,indexes] = max(ratios); %gives me the max ratio per disorganisation (column)
    uniqueIndex = unique(indexes);
    for j = 1:numel(uniqueIndex)
        
        disorgIndex = find(indexes==uniqueIndex(j));
        perimeterMT = 0;
        for k = 1:sum(disorgIndex)
            perimeterMT = perimeterMT+tempROIs(disorgIndex(k)).perimeter; % tirar o +1 e por + disorg(index).perimeter
        end
        % MDI = perimeterMT/length(axon)
        tempROIs(uniqueIndex(j)).MDI = perimeterMT/tempROIs(uniqueIndex(j)).length;
    end
    
    
end
delete(f)
% checkPrint(handles.figure1,handles)
finalName = [finalName '.mat'];
finalROIs = tempROIs;
save(finalName,'finalROIs'); % Saves everything into a file that can be opened later

finalStraightLengths =  cell(numel(tempROIs)-1,1);
finalStraightPos = cell(numel(tempROIs)-1,1);
finalStraightnessIndex = cell(numel(tempROIs)-1,1);
%         tempROIs(handles.imageIndex).straightlineDistances
for i = 2:numel(tempROIs)
    if ~isempty(tempROIs(i).conversionFactor) % O Nuno editou isto, não tenho a certeza que esteja a fazer o que queres
        conversionFactor = tempROIs(i).conversionFactor;
    else
        conversionFactor = 1;
    end
    finalStraightLengths{i-1} = tempROIs(i).straightlineDistances.*conversionFactor;
    finalStraightPos{i-1} = tempROIs(i).straightlineCoord;
    
    if strcmp(tempROIs(i).type,'MT disorganisation') % THESE ARE INDEXES, i don't need the conversion factor
        finalStraightnessIndex{i-1} = sum(tempROIs(i).straightlineDistances)./tempROIs(i).areaMT;
    else
        finalStraightnessIndex{i-1} = sum(tempROIs(i).straightlineDistances)./tempROIs(i).length;
    end
end


%% Saving To an Excel File. 
% csvname = [finalName '.csv'];
% 
% 
% i=2; % Index for what the user wants
% j=1; % index for the finalResults struct
% finalResults = struct([]);
% f = waitbar(0,'1','Name','Saving',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% % checkPrint(handles.figure1,handles)
% while i<=size(saveAnswer,2)
%     finalResults(j).name = saveAnswer(i);
%     switch string(saveAnswer(i))
%         case 'General Curvature'
%             finalResults(j).value = finalGeneralCurvatures;
%         case 'Individual Curvature'
%             finalResults(j).value = finalPointCurvatures;
%         case 'Axonal Length'
%             finalResults(j).value = finalLengths;
%         case 'Straight Segment Lengths'
%             finalResults(j).value = finalStraightLengths;
%         case 'Straight Segment Positions'
%             finalResults(j).value = finalStraightPos;
%         case 'Eccentricity'
%             finalResults.(j).value = finalEccentricity;
%             %                 case 'MT Disorganisation Index'
%             %                     finalResults(j).value = ;
%         case 'Straightness Index'
%             finalResults(j).value = finalStraightnessIndex;
%         case 'MT Area'
%             finalResults(j).value = finalAreaMT;
%         case 'Swelling Area'
%             finalResults(j).value = finalAreasDisorg;
%         case 'Density'
%             finalResults(j).value = finalDensity;
%     end
%     i = i+1;
%     j = j+1;
%     
% end
% delete(f)
% struct2csv(finalResults,csvname);
% 
% msgbox('The requested files have been saved.', 'Saved');
% 
% finalName = [finalName '.mat'];
% finalROIs = tempROIs;
% save(finalName,'finalROIs'); % Saves everything into a file that can be opened later


%% Run straightlines

roiIndex = 3;
straightLines = tempROIs(roiIndex).straightlines;
sizeStraightLines = length(straightLines);
distances=[];
segmentNumber = 1;
distanceIndex = 1;
i=1;
while i<=sizeStraightLines
    try
    if (straightLines(i).rho == straightLines(i+1).rho) && (straightLines(i).theta == straightLines(i+1).theta)
        distances = [distances tempROIs(roiIndex).straightlineDistances(i)+tempROIs(roiIndex).straightlineDistances(i+1)];
        i=i+1;
    else
        distances = [distances tempROIs(roiIndex).straightlineDistances(i)];
    end
    catch
        distances = [distances tempROIs(roiIndex).straightlineDistances(i)];
    end
    i=i+1;
end

straightnessRatio = sum(distances)/tempROIs(roiIndex).length*100;
numberSegments = length(distances);