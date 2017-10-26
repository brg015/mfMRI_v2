function [out]=RSA_bootstrap(subject,R,vox,out,ii,out_name,tmp)    
% subject (str)      -> subjectID
% R (matrix NxN)     -> brain matrix from cross-cor
% ii                 -> model index
% vox (#int)         -> voxel index (for SL models)
% out {cell(volume)} -> brain maps (for SL models)
% out_name           -> file name (for SL models)
% tmp (matrix)       -> brain matrix from volume
%
% NOTE: SL save directly to the output maps, thus why they are passed
% into the function (see Write to Searchlight Volumes)
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
    out_index{1}=[SL.design.save_str{ii} '_key'];
end

% out_array{} has all the values
% out_index{} has the file name to save to
switch SL.design.calc{ii}  
    case 'AtoB'
       [out_index,out_array]=RSA_model_AtoB(R,ii);      
    case 'Anova1'
       [out_index,out_array]=RSA_model_Anova1_v2(R,ii);
    case 'Identity1'
        % Updated to v2 on 10/29 to streamline processing
        [out_index,out_array]=RSA_model_Identity1_v2(R,ii);
    % Added 6/23/17
    case 'MVPA' % SVM only now
        [out_index,out_array]=RSA_model_MVPA_SVM(tmp(:,SL.design.trial_include{ii}),ii); 
    % Added 6/23/17
    case 'distance'
        [out_index,out_array]=RSA_model_distance(tmp(:,SL.design.trial_include{ii}),ii); 
    % Added 10/26/17
    case 'MRegression'
        [out_index,out_array]=RSA_model_regression(R,ii);
end

%=========================================================================%
% Write to searchlight volumes
%=========================================================================%
% Update 12/10 for 4D volumes
% Updated 5/17/17 for custom models
if SL.region.noSL==0
    % Do Stuff
    try
        % Note* if out is set, and out_index is never updated, then this
        % loop is never hit and thus 'out' is fine
        for jj=1:length(out_index), out{strcmp(out_index{jj},out_name)}.mat(vox,:)=out_array{jj}; end
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