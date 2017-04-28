function fJacc = Create_JC_Network_450(group1Flist, group2Flist, folderOut, thresh)

% This code is used on ROI based networks

if nargin < 4
    thresh = 0;
end

if nargin < 3
    error('Need at least 3 inputs!');
end



wfu_mkdir(folderOut);


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

for i = 1:numberSubjects
    fprintf('\t%d of %d\n', i, numberSubjects)
    data1 = load(deblank(dataFiles(i,:)));
    data1 = data1 .* (data1 > thresh);
    for j = 1:i
        fprintf('\t\t%d of %d\n', j, i);
        if j==i
            img2 = img1;
        else
            data2 = load(deblank(dataFiles(i,:)));
            data1 = data2 .* (data2 > thresh);
        end
        
        % Calculate Czekanowski Index
        temp = 1 - (sum(abs((double(data1(:)) - double(data2(:)))))/sum(data1(:) + data2(:)));
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