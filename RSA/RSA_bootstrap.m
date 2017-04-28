function [out]=RSA_bootstrap(subject,R,vox,out,ii,out_name)     
global SL;
out_array={};  % Actual values
out_index={};  % String match
%=========================================================================%
% Correlation Calculations
%=========================================================================%
if strcmp(SL.design.calc{ii},'Spear') || strcmp(SL.design.calc{ii},'Euclid') ... 
        || strcmp(SL.design.calc{ii},'Kendall') || strcmp(SL.design.calc{ii},'Mean')
    v1=SL.design.matrix{ii}(~isnan(SL.design.matrix{ii})); % Design vector
    v2=R(~isnan(SL.design.matrix{ii}));                    % Data vector
    switch SL.design.calc{ii}
        case 'Spear',   RHO=corr([v1,v2],'type','Spearman'); 
                        outval(SL.analysis.multi_boot+1)=RHO(1,2);
        case 'Euclid',  outval(SL.analysis.multi_boot+1)=1/sqrt(sum((v1-v2-mean(v2)).^2));
        case 'Kendall', outval(SL.analysis.multi_boot+1)=ktaua([v1,v2]);
        case 'Mean',    outval(SL.analysis.multi_boot+1)=mean(v2(v1==1));
    end
    out_array{1}=outval(SL.analysis.multi_boot+1);
    out_index={'Key'};
end
% 
switch SL.design.calc{ii}  
    case 'AtoB'
       [out_index,out_array]=RSA_model_AtoB(R,ii);      
    case 'Anova1'
       [out_index,out_array]=RSA_model_Anova1_v2(R,ii);
    case 'Identity1'
        % Updated to v2 on 10/29 to streamline processing
        [out_index,out_array]=RSA_model_Identity1_v2(R,ii);
end

%=========================================================================%
% Write to searchlight volumes
%=========================================================================%
% Update 12/10 for 4D volumes
if SL.region.noSL==0
    % Do Stuff
    try
        % Note* if out is set, and out_index is never updated, then this
        % loop is never hit and thus 'out' is fine
        for jj=1:length(out_index)
            if isnumeric(out_index{jj})
                out{out_index{jj}}.mat(vox,:)=out_array{jj};
            elseif strcmp(out_index{jj},'Key')
                out{ii}.mat(vox,:)=out_array{jj};
            else
                out{strcmp(out_index{jj},out_name)}.mat(vox,:)=out_array{jj};
            end
        end
    catch err
        display('Error writing volumes in RSA_bootstrap');
        keyboard;
    end
%=========================================================================%
% Output Matlab files
%=========================================================================%
elseif SL.region.noSL==1
    [~,ROI,~]=fileparts(SL.region.mask);
    WriteTo=fullfile(SL.dir.outpath,'ROI',subject,SL.design.calc{ii});
    if ~exist(WriteTo,'dir'), mkdir(WriteTo); end
    save(fullfile(WriteTo,[ROI '_' SL.design.save_str{ii} '.mat']),...
        'out_array','out_index');
end

end