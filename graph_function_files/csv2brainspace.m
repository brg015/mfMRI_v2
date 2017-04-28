function csv2brainspace(csvFlist, fOutFlist, atlasImg)
% 
% Parameters:
%     csvFlist:   an flist of csvs that you want to map to brain space
%     fOutFlist:  an flist of image output filenames that correspond to the csv flist
%     atlasImg:   the nii file that shows how each of the 449 nodes map to the brain
%     

% this is essentially all that the LCBN's flist reader does
CSVs = textread(csvFlist,'%s','commentstyle','matlab');
CSVs = char(CSVs);

fOuts = textread(fOutFlist,'%s','commentstyle','matlab');
fOuts = char(fOuts);

% find the locations that each of your 449 regions make up on the voxel
% level.
atlas = load_untouch_nii(atlasImg);
for i = 1:max(atlas.img(:))
    atlasRegions{i} = find(atlas.img(:) == i);
end


for i = 1:size(CSVs,1)
    data = load(deblank(CSVs(i,:)));
    tempImg = atlas; % make a copy of each image
    % for each value in the csv, put in the appropriate locations based on
    % atlas location.
    for ii = 1:length(data)
        tempImg.img(atlasRegions{ii}) = data(i);
    end
    save_untouch_nii(tempImg, deblank(fOuts(i,:)));
    
end
        