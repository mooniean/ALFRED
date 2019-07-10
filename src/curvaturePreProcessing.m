function skeleton=curvaturePreProcessing(coordX,coordY)

%build a binary image out of the coordinates
offset=20; %zero-pad image so I have like a frame around it
skeleton = zeros(max(coordX)+offset,max(coordY)+offset);
for i = 1:length(coordX)
    skeleton(coordX(i)+offset/2,coordY(i)+offset/2) = 1;
end

%find branchpoints
brPoints=bwmorph(skeleton,'branchpoints');

%eliminate branchpoints and their 24 nearest neighbours from the skeleton
[rowBr,colBr] = find(brPoints);
for k=1:length(rowBr)
    for i=-2:2
        for j=-2:2
            skeleton(rowBr(k)+i,colBr(k)+j)=0;
        end
    end
end

%clean groups of up to 10 pixels that are isolated in the skeleton
skeleton=bwareaopen(skeleton,10);