function  trans = transpose_mat(csvfile,outfolder)
% Transpose matrix
file = csvread(csvfile);
name = strsplit(csvfile,'/');
name = name{end};
name = name(1:end-4);

name1 = strcat(name,'_untransposed.csv');
name = strcat(name,'.csv');
% newfile = [outfolder, name];
M = [];
M = transpose(file);
csvwrite(fullfile(outfolder,name1),file);
csvwrite(fullfile(outfolder,name),M);


    