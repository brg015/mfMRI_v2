function [behav] = make_behav_mat(SPM,SL,conds,b_filter)

    if ~isempty(SL.b_filter)
        behav = zeros(1,length(SPM.Vbeta));

        for i = 1:length(SPM.Vbeta)
            if wild_strfind(SPM.Vbeta(i).descrip,SL.b_filter) && SL.conds(i) ~= 0
                behav(i) = 1;
            end
        end
    else
        behav = ones(1,length(SPM.Vbeta));
    end
    
%     if t_mat_save ==1 
%         save(strcat('b_behav_',b_filter),'behav');
%     end
    
