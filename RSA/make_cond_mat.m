function [conds] = make_cond_mat(SPM,SL)
        
    conds = zeros(length(SL.interest),length(SPM.Vbeta));
    
    for ncond = 1:length(SL.interest)

        convect = zeros(1,length(SPM.Vbeta));
        
        for b = 1:length(SPM.Vbeta)
            bool = wild_strfind(SPM.Vbeta(b).descrip,SL.interest{ncond});
            if bool == 1
               convect(b) =  str2num(SPM.Vbeta(b).descrip(strfind(SPM.Vbeta(b).descrip,SL.id_start)+length(SL.id_start):strfind(SPM.Vbeta(b).descrip,SL.id_end)-1));
            end
        end
        
        new = convect;
        new(new == 0) = [];
        
        for c = 1:length(new)
            for i = 1:length(SPM.Vbeta)
                if ~isempty(strfind(SPM.Vbeta(i).descrip,strcat(SL.id_start,num2str(new(c)),SL.id_end)))
                    convect(i) = new(c);
                    if length(find(convect == new(c))) == 2
                        break
                    end
                end
            end
        end
        
%         if t_mat_save == 1
%             name = interest{ncond};
%             name(strfind(name,'*')) = [];
%             save(strcat('E1_',name),'convect');
%         end
        conds(ncond,:) = convect;
        
    end