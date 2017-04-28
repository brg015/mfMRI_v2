function an = run_multi(Flist, Nodemat, Outfolder,tranposemat1,tranposemat2,hasNaN,sameSubjectFlag,runwhole)

% runs multiple select permutations on particular nodes
% NOTE: changed the name(1:end-8) value to appropiate size of the file name
% you want to have
%Flist = a xls file to path of matrices csv files(not compatiable with mat). Put in two column to compare the
%two matrices
%Nodemat = a xls file with list of nodes to be used in a permutation test,
%each all nodes in a column would be used in a permutation test, 0 if you
%want to runwhole
%Outfolder = path for results/graphs to be saved
% transposemat1 = 1 if the matrix 1 needs to be transposed
% transposemat2 = 1 if the matrix 2 needs to be transposed
% hasNaN = 1 if need to clear NaN values
% SameSubjectFlag = 1 if the two sets are from same subject group,0
% otherwise
% runwhole = 0 for particular nodes, 1 to run with all nodes thus making
% nodemat necessary
[numeric, strings]= xlsread(Flist);
count = 1;
for i = 1: size(strings,1)
    
    epi = csvread(strings{i,1});
    sem = csvread(strings{i,2});
    
    
    name = strsplit(strings{i,1},'\');
    name = name{end};
    % end 
    name = name(1:end-8);
    EpiMat = [];
    SemMat = [];
    MAT{1,1} = 'name';
    MAT{1,2} = 'pvalue';
    
    if(runwhole ==0)
    [nodes,txt] = xlsread(Nodemat);
    %     [nodes] = csvread(Nodemat);
    end
    if(tranposemat1 == 1)
        epi = transpose(epi);
    end
    if(tranposemat2 == 1)
        sem = transpose(sem);
    end
    if(runwhole ==0)
    for j = 1:size(nodes,2)
        for k = 1:size(nodes,1)
      
            if nodes(k,j) == 0
                sizemat = k-1;
                break
            else
             
                    EpiMat = cat(2,EpiMat,epi(:,nodes(k,j)));
                    SemMat = cat(2,SemMat,sem(:,nodes(k,j)));
                
            end
        end
         name1 = strcat(name,txt{1,j});
        %name1=name;
        pval = run_permutation(EpiMat,SemMat,sameSubjectFlag,Outfolder,hasNaN,name1);
        MAT{count+1,2} = pval;
        MAT{count+1,1} = name1;
        EpiMat = [];
        SemMat = [];
        count= count+1;
        
    end     
    end
    if(runwhole ==1)
        pval = run_permutation(sem,epi,sameSubjectFlag,Outfolder,hasNaN,name);
        MAT{count+1,2} = pval;
        MAT{count+1,1} = name;
        count =count +1;
    end
end
an = MAT;

