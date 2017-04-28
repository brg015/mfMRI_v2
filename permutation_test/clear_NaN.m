function clearmat = clear_NaN(mat1,folderOut)
% Method to clear out any NaN from the adjacency matrices in order to run
% further calculations on them

k=1;
    removednode = [folderOut, '/removed_nodes.mat'];
    for c= 1:size(mat1,2)
        c= c- k+1;
        col = mat1(:,c)
        d = col(1,1);
        e =mat2str(d);
        e = e(1);
        if e == 'N'
            M(k,1) = c;
            mat1(:,c) = [];
            k = k+1;
            
        end
    end
     save(removednode, 'M')
     clearmat = mat1;