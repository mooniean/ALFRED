function state = skeletonOverlap(boxA,boxB,skelA,skelB)
% Input: Bounding Boxes boxA and boxB (coordinates of the top left corner of the boxes, width and height), image skeletons skelA and skelB.
% Output: Boolean State (true or false)
% This function evaluates whether skelA and skelB overlap with each other. 
% The skeletons are comprised inside boxA and boxB, respectively, in the original image. 
% The function will create a virtual overlap matrix with the size of the boxes and there is a
% coordinate transformation from "original image" to "new image". 
% As soon as (or if) the state variable becomes true, the function ends and returns that value.

state=false;

% 'inside skeleton overlap'
% newImage = zeros(max(boxA(3)+boxA(1),boxB(3)+boxB(1)),max(boxA(4)+boxA(2),boxB(4)+boxB(2)));
newImage = zeros(max(boxA(4)+boxA(2),boxB(4)+boxB(2)),max(boxA(3)+boxA(1),boxB(3)+boxB(1)));

% coordx=boxA(1)-min(boxA(3),boxB(3));
% coordy=boxA(2)-min(boxA(4),boxB(4));
% coordx=boxA(1)-min(boxA(1),boxB(1));
% coordy=boxA(2)-min(boxA(2),boxB(2));

coordx=boxA(2)-min(boxA(2),boxB(2));
coordy=boxA(1)-min(boxA(1),boxB(1));
for i=1:size(skelA,1)
    if max(skelA(i,:))<1
        continue
    end
    for j=1:size(skelA,2)
        if skelA(i,j)==1
            newImage(coordx+i,coordy+j)=newImage(coordx+i,coordy+j)+1;
        end
    end
end


% 
% coordx=boxB(1)-min(boxA(3),boxB(3));
% coordy=boxB(2)-min(boxA(4),boxB(4));

% coordx=boxB(1)-min(boxA(1),boxB(1));
% coordy=boxB(2)-min(boxA(2),boxB(2));

coordx=boxB(2)-min(boxA(2),boxB(2));
coordy=boxB(1)-min(boxA(1),boxB(1));

for i=1:size(skelB,1)
    if max(skelB(i,:))<1
        continue
    end
    for j=1:size(skelB,2)
        if skelB(i,j)==1
            if newImage(coordx+i,coordy+j)==1
                state = true;
                break
            end
        end
    end
end


end