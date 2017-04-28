function p = Permutationtest_jaccard(imgList, numGroup1, numGroup2, thresh, sameSubjectsFlag, numPerm, fOut, maskOutFlag, ROIList)
%
%   This function is a way to compare two groups of subjects using each
%   individual subject's data image (ie eGlob.nii, eLoc.nii, Ki.nii).
%   This function utilizes permutationTest.m
%
%   Parameters: imgList    - the flist of data images for each subject. The
%                            first images in the list are in group 1.
%               numGroup1  - the number in Group 1
%               numGroup2  - the number in Group 2
%               thresh     - Threshold for data images. Only needed if
%                            data images are not binary.
%                            Default: 0
%               sameSubjectsFlag - only use if comparing
%                                  different scans of same group
%                                  Deafult: 0
%               numPerm     - Number of permutations for the jaccard matrix.
%                             Default: 100000
%               fOut        - folder output location
%               maskOutFlag - create a mask of only voxels where each image
%                             has a value and apply mask to each brain
%                             Default: 1
%               ROIList     - Region Of Interest Image flist (one for each group)
%                             Used to compare only certain regions between
%                             groups. Not proven to work.
%                             * Not used, might be taken out
%                   
%   Output:         One file is saved with the name
%                   'Jaccard_[imgList].mat'. It stores the following:
%                       JI_all  : Statistical Matrix
%                       p       : p Value
%                       numPerm : number of permuations
%                   pVals.mat - contains p for each threshold passed in and
%                               actualThresholds 
%                   group1_[thresh].nii overlap of images at each threshold
%                               for each group.
%
%
%   Robert Lyday - 4/30/12
%   Updated for distribution on 6/20/12

%**************************************************************************
% Deal with function parameters
%**************************************************************************
if nargin < 9 || ROIList == 0
    overlapImgs = 0;
else
    overlapImgs = wfu_read_flist(ROIList);
end
if nargin < 8
    maskOutFlag = 1;
end
if nargin < 6
    numPerm = 100000;
end
if nargin < 5
    sameSubjectsFlag = 0;
end
if nargin < 4
    thresh = 0;
    fprintf('No threshold given, assuming all matricies are binary!\n');
end
if nargin < 3
    error('Need at least 3 arguments, imgList, numGroup1, numGroup2\n');
end

%**************************************************************************
% Preprocess images
%**************************************************************************
imgs = wfu_read_flist(imgList);
if nargin < 7
    fOut = wfu_fileparts(imgList);
end
[~, fn] = wfu_fileparts(imgList);
%********************************
% ROI Overlap
%********************************
if overlapImgs
    % load image
    temp = load_untouch_nii(deblank(overlapImgs(1,:)));
    
    % create mask for everywhere not equal to 1
    overlaps{1} = temp.img;
    % save mask of everywhere equal to 1
    fileName = strcat(fOut,'/group1.nii');
    
    temp.img = 1 - overlaps{1};
    save_untouch_nii(temp,fileName);
    temp = load_untouch_nii(deblank(overlapImgs(2,:)));
    
    overlaps{2} = temp.img;
    fileName = strcat(fOut,'/group2.nii');
    temp.img = 1 - overlaps{2};
    save_untouch_nii(temp,fileName);
        overlapFinal = (overlaps{1} + overlaps{2})/2;
else
    overlapFinal = 1;
end

% The intersectMask is used as a standard for thresholding
% this image is for images with voxel sizes of 4x4x5 mm. If another size is
% needed you will need to reslice this image.
intersectMask = load_untouch_nii('intersectMask.nii'); 
intersectMaskNum = sum(intersectMask.img(:) > 0);

%********************************
% Mask out
%********************************
if maskOutFlag
    for i = 1:size(imgs,1)
        tempImgFn = deblank(imgs(i,:));
        tempImg = load_untouch_nii(tempImgFn);
        if i == 1
            totalImg = tempImg;
            totalImg.img = totalImg.img > 0;
        else
            totalImg.img = totalImg.img + (tempImg.img > 0);
        end
    end
    % Save overlap image just to have it
    mkdir(fOut);
    tempFN = strcat(fOut,'/overlapped.nii');
    save_untouch_nii(totalImg, tempFN);
    % Remove anything > 0 that is not == max 
    % (only keep what exists in every image);
    tempImg1 = totalImg.img ~= max(totalImg.img(:));
    tempImg2 = totalImg.img > 0;
    totalImg.img = (tempImg1 + tempImg2) == 2;
    maskOutImg = 1 - totalImg.img;
    intersectMask.img = intersectMask.img - totalImg.img;
    intersectMask.img = intersectMask.img == 1;
    intersectMaskNum = sum(intersectMask.img(:));
    totalImg.img = maskOutImg;
    maskName = strcat(fOut, '/mask.nii');
    save_untouch_nii(totalImg, maskName);
else
    maskOutImg = 1;
end


%********************************
% Threshold
%********************************
pCounter = 0;

for ii = 1:length(thresh)
    threshPerc = thresh(ii);
    if threshPerc > 1
        threshPerc = threshPerc/100;
    end
    for i = 1:size(imgs,1)
        tempImgFn = deblank(imgs(i,:));
        tempImg = load_untouch_nii(tempImgFn);
        
        if threshPerc
            
            % Make image a vector
           
            tempImg.img = tempImg.img .* maskOutImg; % Apply mask
            tempImg.img = tempImg.img .* overlapFinal; % Apply overlap
            imgVect = tempImg.img(:);

            % Sort Vector
            imgVect = sort(imgVect,'descend');
            
            % Find location of thresh percentage
            threshLoc = intersectMaskNum * threshPerc;% ** OLD THRESHOLD METHOD -> round(length(imgVect)*threshPerc);
            % Find Value of vector at threshLoc
            threshVal = imgVect(round(threshLoc));
            % Binarize img
            tempImg.img = tempImg.img > threshVal;
            
     
            threshValues(i) = sum(tempImg.img(:))/sum(imgVect>0);
            
        end
       
        tempFN = strcat(fOut, '/thresholded_', num2str(i),'.nii');
        save_untouch_nii(tempImg, tempFN);
        tempMask = maskOutImg .* intersectMask.img;
        
        temp = tempImg.img(:);
        
        
        clus{i}=temp;
        
        
        
    end
    % Save Thresholds

    threshFn = strcat(fOut, '/thresh',num2str(threshPerc*100), '.mat');
    save(threshFn, 'threshValues');
    actualThresholds(ii) = mean(threshValues);
    
    fprintf('\n\n Files preprocessed!\n');
    fprintf('Expected threshold %.2f \n Actual Threshold %.2f \n', threshPerc, mean(threshValues));

%**************************************************************************
% Create Jaccard
%**************************************************************************
    fprintf('Creating Jaccard\n');
    for i = 1:length(clus)
        for j = 1:i
            N11 = clus{i} + clus{j};
            N11 = N11 == 2;
            N11 = sum(N11);
            d1 = sum(clus{i});
            d2 = sum(clus{j});
            
            
            JI_all(i,j) = (N11 ./ (d1 + d2 + - N11));
            JI_all(j,i) = JI_all(i,j);
            
            % Ensure that the Jaccard Index is between 0 and 1
            if JI_all(i,j) > 1
                fprintf('Jaccard Index Too big at %d, %d\n',i,j);
            elseif JI_all(i,j) < 0
                fprintf('Jaccard Index too small at %d, %d\n', i, j);
            end
            
        end
    end
    
    JfileName = strcat(fOut,'/Jaccard_',fn,'.mat');
    save(JfileName,'JI_all');
%     if threshPerc == .2
%         Origname = strcat(fOut,'/Jaccard_20.mat');
%         save(Origname,'JI_all');
%     end
    Origname = strcat(fOut,'/Jaccard_',num2str(threshPerc *100),'.mat');
    save(Origname,'JI_all');
    
%**************************************************************************
% Permutation
%**************************************************************************  
    fprintf('Permutation...\n');
    p(ii) = permutationTest(JfileName, numGroup1, numGroup2, sameSubjectsFlag, numPerm,1);
    
    % If average threshold value == 1, increase counter
    if mean(threshValues) == 1 
        pCounter = pCounter + 1;
    end
    fprintf('p Value: %.2f\n',p(ii));
    
    % *******************************
    %   Create Group Overlap Images
    %*********************************
    
    % Group 1
    totalImg = 0;
    for i = 1:numGroup1
        tempImgFn = deblank(imgs(i,:));
        tempImg = load_untouch_nii(tempImgFn);
        tempImg.img = tempImg.img .* maskOutImg;
        tempImg.img = tempImg.img .* overlapFinal; % Not really used
        imgVect = tempImg.img(:);
        imgVect = sort(imgVect,'descend');
        threshLoc = round(intersectMaskNum * threshPerc);
        threshVal = imgVect(threshLoc);
        tempImg.img = tempImg.img > threshVal;

        totalImg = totalImg + tempImg.img;
    end
    tempImg.img = totalImg/numGroup1;
    sigFN = strcat(fOut,'/group1_', num2str(round(threshPerc*100)),'.nii');
    save_untouch_nii(tempImg, sigFN);
    
    % Group 2
    totalImg = 0;
    for i = numGroup1+1:size(imgs,1)
        tempImgFn = deblank(imgs(i,:));
        tempImg = load_untouch_nii(tempImgFn);
        tempImg.img = tempImg.img .* maskOutImg;
        tempImg.img = tempImg.img .* overlapFinal;
        imgVect = tempImg.img(:);
        imgVect = sort(imgVect,'descend');
        threshLoc = round(intersectMaskNum * threshPerc);
        threshVal = imgVect(threshLoc);
        tempImg.img = tempImg.img > threshVal;

        totalImg = totalImg + tempImg.img;
    end
    tempImg.img = totalImg/numGroup2;
    
    sigFN = strcat(fOut,'/group2_', num2str(round(threshPerc*100)),'.nii');
    save_untouch_nii(tempImg, sigFN);
    
    
    
    if pCounter > 1
        break; % end loop if you get 1 twice
    end
end

fn = strcat(fOut , '/pVals.mat');
save(fn, 'p', 'actualThresholds');


