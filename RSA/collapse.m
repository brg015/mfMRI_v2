function [indx] = collapse(conds)
%%    
    if size(conds,1) == 1
        indx = conds;
    else
        indx = sum(conds);
    end
