function an = run_allnodes(Flist, Outfolder,tranposemat1,tranposemat2,hasNaN,sameSubjectFlag)

% runs multiple select permutations on ALL individual nodes, shotgun
% approach
%Flist = a xls file for list of csvfiles to run
% transposemat1 = 1 if the matrix 1 needs to be transposed
% transposemat2 = 1 if the matrix 2 needs to be transposed
% hasNaN = 1 if need to clear NaN
% SameSubjectFlag = 1 if the two sets are from same subject group

[numeric, strings]= xlsread(Flist);
count = 1;
MAT1 = [];
for i = 1: size(strings,1)
    
    epi = csvread(strings{i,1});
    sem = csvread(strings{i,2});
    
    
    name = strsplit(strings{i,1},'\');
    name = name{end};
    name = name(1:end-8);

    

    %     [nodes] = csvread(Nodemat);
    if(tranposemat1 == 1)
        epi = transpose(epi);
    end
    if(tranposemat2 == 1)
        sem = transpose(sem);
    end
    for j = 1:size(epi,2)
       
        fprintf('num: %s %i',name , j);
        pval = run_permutation(epi(:,j),sem(:,j),sameSubjectFlag,Outfolder,hasNaN,name);
        MAT{j,1} = pval;      
    end     
    fileName = strcat(Outfolder,'/',name,'pval.csv');
    csvwrite(fileName,MAT);
    MAT = [];
    MAT1(:,i) = MAT;
end
an = MAT1;
end


