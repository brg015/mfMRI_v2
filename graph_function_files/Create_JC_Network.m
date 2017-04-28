function fJacc = Create_JC_Network(group1Flist, group2Flist, folderOut, unionMask, thresh, ROI, voxDim)
if nargin < 7
    voxDim = 1;
end

if nargin < 6 || isempty(ROI)
    ROI.img = 1;
else
    ROI = load_untouch_nii(ROI);
end

if nargin < 5
    thresh = 0;
end

if nargin < 4
    unionMask = 1;
end

if nargin < 3
    error('Need at least 2 inputs!');
end



wfu_mkdir(folderOut);
if voxDim == 1
    vMask = load_untouch_nii('/lcbn2/lcbn/shared_images/intersectMask_444.nii');
elseif voxDim == 2
    vMask = load_untouch_nii('/lcbn2/lcbn/shared_images/intersectMask_445.nii');
elseif voxDim == 3
    vMask = load_untouch_nii('/victer4/ABC/SPM_Stats/intersectMask_1.5.nii');
elseif voxDim ==4
    vMask = load_untouch_nii('/lookahead1/1_Jakob_final/jakob_atlas_flipped_final_4mm_binary.nii');
end
vMaskPointer = find(vMask.img);



% create single flist
gFN = [folderOut, '/combined.flist'];
fid = fopen(gFN, 'w');
dataFiles = wfu_read_flist(group1Flist);
for i = 1:size(dataFiles,1)
    fprintf(fid, '%s\n', dataFiles(i,:));
end
dataFiles = wfu_read_flist(group2Flist);
for i = 1:size(dataFiles,1)
    fprintf(fid, '%s\n', dataFiles(i,:));
end
fclose(fid);

dataFiles = wfu_read_flist(gFN);


numberSubjects = size(dataFiles,1);

Origname = [folderOut, '/MI_all_JC.mat'];
newMask = vMask;
newMask.img = newMask.img * 0;
for i = 1:numberSubjects
    img1 = load_untouch_nii(deblank(dataFiles(i,:)));
    newMask.img = newMask.img + (img1.img > thresh);
end
newMask.img = newMask.img > 0;
if unionMask
    vMask = newMask;
    vMaskPointer = find(vMask.img>0);
end


for i = 1:numberSubjects
    img1 = load_untouch_nii(deblank(dataFiles(i,:)));
    img1.img = img1.img .* (img1.img > thresh) .* ROI.img;
    fprintf('\t%d of %d\n', i, numberSubjects)
    for j = 1:i
        fprintf('\t\t%d of %d\n', j, i);
        if j==i
            img2 = img1;
        else
            img2 = load_untouch_nii(deblank(dataFiles(j,:)));
            img2.img = img2.img .* (img2.img > thresh);
        end
        vals1 = img1.img(vMaskPointer);
        vals2 = img2.img(vMaskPointer);
        % Calculate Czekanowski Index
        temp = 1 - (sum(abs((double(vals1(:)) - double(vals2(:)))))/sum(vals1(:) + vals2(:)));
        % Jaccardize it
        temp = temp/(2- temp);
        MI_all(i,j) = temp;
        MI_all(j,i) = MI_all(i,j);
        
    end
end
JI_all = MI_all;
JI_orig = JI_all;
save(Origname,'JI_all','JI_orig');
fJacc = Origname;