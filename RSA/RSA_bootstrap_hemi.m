function [out,noise_out]=RSA_bootstrap_hemi(subject,R,tmp,vox,out,ii,out_name,noise_out)     
global SL;
out_array={};  % Actual values
out_index={};  % String match
%=========================================================================%
% Correlation Calculations
%=========================================================================%
[out_index,out_array]=RSA_model_Identity1_hemi(R,tmp,ii,out_name);
% if strcmp(SL.design.calc{ii},'Spear') || strcmp(SL.design.calc{ii},'Euclid') ... 
%         || strcmp(SL.design.calc{ii},'Kendall')
%     v1=SL.design.matrix{ii}(~isnan(SL.design.matrix{ii})); % Design vector
%     v2=R(~isnan(SL.design.matrix{ii}));                    % Data vector
%     switch SL.design.calc{ii}
%         case 'Spear',   RHO=corr([v1,v2],'type','Spearman'); 
%                         outval(SL.analysis.multi_boot+1)=RHO(1,2);
%         case 'Euclid',  outval(SL.analysis.multi_boot+1)=1/sqrt(sum((v1-v2-mean(v2)).^2));
%         case 'Kendall', outval(SL.analysis.multi_boot+1)=ktaua([v1,v2]);
%     end
%     out_array{1}=outval(SL.analysis.multi_boot+1);
%     out_index={'Key'};
% end
% 
% switch SL.design.calc{ii}  
%     case 'AtoB'
%        [out_index,out_array]=RSA_model_AtoB(R,ii);      
%     case 'Univariate'
%        [out_index,out_array]=RSA_model_univariate(tmp,ii);
%     case 'Anova1'
%        [out_index,out_array]=RSA_model_Anova1(R,ii);
%     case 'mvpa' % IN WORK
%         RSA_model_mvpa();
%     case 'ER_special'
%         out=RSA_model_ER_special(R,vox,tmp,out);
%     case 'Identity1'
%         [out_index,out_array]=RSA_model_Identity1(R,tmp,ii,out_name);
% end

%=========================================================================%
% Write to searchlight volumes
%=========================================================================%
% Update 12/10 for 4D volumes
% if SL.region.noSL==0
    % Do Stuff
%     try
        % Note* if out is set, and out_index is never updated, then this
        % loop is never hit and thus 'out' is fine
        for jj=1:length(out_index)
            if isnumeric(out_index{jj})
                out{out_index{jj}}.mat(vox,:)=out_array{jj};
%             elseif strcmp(out_index{jj},'Key')
%                 out{ii}.mat(vox,:)=out_array{jj};
            else
%                 try
                    out{strcmp(out_index{jj},out_name)}.mat(vox,:)=out_array{jj};
%                 catch err
%                     continue;
%                 end
            end
        end
%     catch err
%         display('Error writing volumes in RSA_bootstrap');
%         keyboard;
%     end
%=========================================================================%
% Output Matlab files
%=========================================================================%
% elseif SL.region.noSL==1
%     [~,ROI,~]=fileparts(SL.region.mask);
%     WriteTo=fullfile(SL.dir.outpath,'ROI',subject,SL.design.calc{ii});
%     if ~exist(WriteTo,'dir'), mkdir(WriteTo); end
%     save(fullfile(WriteTo,[ROI '_' SL.design.save_str{ii} '.mat']),...
%         'out_array','out_index');
% end
% 
% if isfield(SL.design,'noise')
%     if SL.design.noise(ii)==1
%         noise_out{ii}.mat(vox,:)=v2;
%     end
% end

% % Compute bootstrap randomization vectors - dependent upon the vector,
% % which is dependent upon the design matrix.
% if SL.analysis.multi_boot>0
%     SL.analysis.boot=zeros(SL.analysis.multi_boot,N);
%     % May be able to initialize this earlier - nope, vector dependent
%     % length, so would occasionally fail.
%     for jj=1:SL.analysis.multi_boot
%         SL.analysis.boot(jj,:)=randperm(N); % [permute X order] matrix
%     end
%     % Randomize the design matrix half the time and the data half the time
%     for jj=1:SL.analysis.multi_boot/2
%         switch SL.design.calc{ii}
%             case 'Spear'
%                 RHO=corr([v1(SL.analysis.boot(jj,:)),v2],'type','Spearman'); 
%                 outval(jj)=RHO(1,2);
%             case 'Euclid'
%                 outval(jj)=...
%                     1/sqrt(sum((v1(SL.analysis.boot(jj,:)-v2-mean(v2)).^2)));
%         end
%     end
%     for jj=SL.analysis.multi_boot/2+1:SL.analysis.multi_boot
%         switch SL.design.calc{ii}
%             case 'Spear'
%                 RHO=corr([v1,v2(SL.analysis.boot(jj,:))],'type','Spearman'); 
%                 outval(jj)=RHO(1,2);
%             case 'Euclid'
%                 outval(jj)=...
%                     1/sqrt(sum((v1-v2(SL.analysis.boot(jj,:)-mean(v2)).^2)));
%         end
%     end
%     % Calculate values and save to output
%     sd=std(outval(1:SL.analysis.multi_boot));
%     m=mean(outval(1:SL.analysis.multi_boot));
%     anti_pval=sum(v <= m)/SL.analysis.multi_boot;
%     pos_pval=sum(v >= m)/SL.analysis.multi_boot;
%     out{2}(vox,ii)=anti_pval;
%     out{3}(vox,ii)=pos_pval;
%     out{4}(vox,ii)=m;
%     out{5}(vox,ii)=sd;
% end

end