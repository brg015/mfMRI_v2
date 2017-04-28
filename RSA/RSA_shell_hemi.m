function RSA_shell_hemi()
%=========================================================================%
%% Code Setup
%=========================================================================%
% Updates 11/12/2014
% 1) Added many presets to be internal, for example, all of SL.analysis and
% all of SL.preprocess. Plus a few others, basically, simplified the Master
% script by making more internal assumptions
global SL FIRm;
keyboard;
% Setting fields if they do not exist
bugger=0;
%=========================================================================%
% Defunct options
%=========================================================================%
if ~isfield(SL,'preprocess'),       SL.preprocess.on=0; end
if ~isfield(SL,'noise'),            SL.noise.on=0; end
%=========================================================================%
% Default options
%=========================================================================%
SL.region.noSL=SL.region.use_mask;  % Set equal as default
if ~isfield(SL,'ID'),             SL.ID.check=0; end
if ~isfield(SL,'regress'),        SL.regress.on=0; end
if SL.regress.on==1
    if ~isfield(SL.regress,'type');       SL.regress.type='normal'; end
    if ~isfield(SL.regress,'cond');       SL.regress.cond=[]; end
end
if ~isfield(SL.design,'mumford'), SL.design.mumford=0; end
if ~isfield(SL,'err'),            SL.err=0; end
if ~isfield(SL,'fisher'),         SL.fisher=0; end % Fisher transform R

if ~isfield(SL,'run'),            SL.run.include=0; end

if ~isfield(SL.region,'noModel'), SL.region.noModel=0; end
if ~isfield(SL.region,'noSL'),    SL.region.noSL=0; end

if ~isfield(SL.dir,'QA'),         SL.dir.QA='QA'; end
if ~isfield(SL.dir,'check'),      SL.dir.check=1; end
if ~isfield(SL.dir,'overwrite'),  SL.dir.overwrite=0; end

% If no design is set, then mark is as null
if ~isfield(SL.design,'interactions'); SL.design.interactions{1}=[]; end
if ~isfield(SL.design,'model');        SL.design.model=[]; end
if ~isfield(SL.design,'matrix');       SL.design.matrix{1}=[]; end
if ~isfield(SL.design,'save_str');     SL.design.save_str{1}=''; end
if ~isfield(SL.design,'custom');       SL.design.custom=0; end
if ~isfield(SL.design,'calc');         SL.design.calc{1}=''; end

% Prevents the setupt of output maps in advance
if ~isfield(SL,'skip_maps');           SL.skip_maps=0; end
%=========================================================================%
% Special Analysis Strings
%=========================================================================%
% Temporary override variable to specially alter functions for me
if ~isfield(SL,'geib'),           SL.geib=0; end

% ER Match special analyses
if ~isfield(SL,'ERmatch');        SL.ERmatch=0; end
if SL.ERmatch==1, SL.skip_maps=1; end

if ~isfield(SL,'Hemi');           SL.Hemi=0; end
if ~isfield(FIRm,'on');           FIRm.on=0; end
if ~isfield(SL,'GBU');            
    SL.GBU=0; 
    NModels=length(SL.design.save_str);
else
    NModels=7;
end
% Current functionality 11/12
% -> Runs check and analysis at the same time
%=========================================================================%
%% SL.analysis
%=========================================================================%
%   .sd_trim (num)         => Number of standard deviations to trim 
%   .voxel_per (num)       => Minimum number of voxels to include in SL
%   .multi_boot (num)      => Bootstrap iterations (Set ==0 for now)
%                          => This is currently not implemented
%   .sparse_sampling (num) => Allows for sampling of partial data
%   .on (num)              => Never used, just easy way to init.
if ~isfield(SL,'analysis'),                 SL.analysis.on=1; end 
if ~isfield(SL.analysis,'sd_trim'),         SL.analysis.sd_trim=3; end
if ~isfield(SL.analysis,'voxel_per'),       SL.analysis.voxel_per=6; end
if ~isfield(SL.analysis,'multi_boot'),      SL.analysis.multi_boot=0; end
if ~isfield(SL.analysis,'sparse_sampling'), SL.analysis.sparse_sampling=1; end
%=========================================================================%
SL.err=0;     % No error currently exists
TrialLimit=0; % Says min. number of trials per condition

if ~exist(SL.dir.outpath,'dir'), mkdir_tree(SL.dir.outpath); end

% Setup RSA models (custom)
if isfield(SL,'model'), SL=RSA_model_design(SL); end
% Setup mask list if needed
if isfield(SL.region,'mask'), mask_list=SL.region.mask; else mask_list=[]; end
%=========================================================================%
%% Code Start
%=========================================================================%
fprintf(strcat('VECTORIZED SEARCHLIGHT\t',datestr(clock),'\n'));
for cursub = 1:length(SL.dir.subjects)
    try
    %=====================================================================%
    % Setup design matrix space
    %=====================================================================%
    % Set flag and subject specific parameters
    flag=0;
    SL.design.ID_idx=[];
    SL.design.ID_descrip={};
    SL.design.ID_file={};
    SL.design.Box=[];
    SL.files=[];
    
    % CD into subject directory
    cd(strcat(SL.dir.stpath,SL.dir.subjects{cursub}));
    % Alert user of loading & do so
    fprintf(strcat('....Loading\t',SL.dir.subjects{cursub},'''s SPM.mat\n'));
    if exist('SPM.mat','file')
        load('SPM.mat');
    else
        if SL.design.mumford==0,
            fprintf('\tSubject has no SPM.mat\n'); continue;
        else
            SPM=[];
        end
    end
    
    %=====================================================================%
    %% Find and load beta files
    %=====================================================================%
    SPM=RSA_find_betas(SPM,strcat(SL.dir.stpath,SL.dir.subjects{cursub}),...
        SL.dir.subjects{cursub});
    SL.V=spm_vol(SL.design.ID_file{1});
    fprintf(strcat('....Loading\t',SL.dir.subjects{cursub},'''s beta images\n'));
    % Dont load files as of yet
    if SL.geib==0
        [SL.files] = load_st_betas(SL);
    end
    %=====================================================================%
    %% Quick data check:
    %=====================================================================%
    display(['....Checking ' SL.dir.subjects{cursub} ' trial counts']);
    for ii=length(SL.design.Box)+1:length(SL.design.cond_str), SL.design.Box(ii)=0; end
    for ii=1:length(SL.design.cond_str)
        t_str=SL.design.cond_str{ii}{1};
        for jj=2:length(SL.design.cond_str{ii})
            t_str=[t_str '_' SL.design.cond_str{ii}{jj}];
        end
        display(['........Condition: ' t_str ' - ' num2str(SL.design.Box(ii))]);
    end
    % Removed conditions that have no trials (impl. for SWAT 6/2/15)
    % c0=find(SL.design.Box==0);
    % SL.design.cond_str(c0)=[];
    
    
    if ~isempty(find(SL.design.Box<TrialLimit))
        display('....Lacking trials, skipping'); continue; 
    end
    %=====================================================================%
    %% Arrange design matrices
    %=====================================================================% 
    % Define
    % SL.design.matrix => design matrix
    % SL.design.power  => matrix intensity
    keyboard;
    for ii=1:NModels
        SL.design.matrix{ii}=[];   % Initialize model
        SL.design.power{ii}=[];
        SL=RSA_design_matrix(SL,ii); 
        % Identity 1 model fix
        if strcmp(SL.design.calc{ii},'Identity1')
            if ~isfield(SL.design.anova{ii},'row')
                SL.design.anova{ii}.row=0;
            end
        end
    end

    if flag==1, 
        fprintf([SL.dir.subjects{cursub} ' FAILED\n']);
        if bugger==1, keyboard; else continue; end
    end
    
    %=====================================================================%
    %% Setup searchlight
    %=====================================================================%
    % Create searchlight indices
    fprintf(strcat('....Creating\t',SL.dir.subjects{cursub},'''s Search Light indices\n'));

    %=====================================================================%
    %% Remove Outliers
    %======================================================================%
    SL.analysis.good_vox=find(~isnan(nanmean(SL.files,2))==1);
    if ~isempty(SL.analysis.sd_trim)
        for ii=SL.analysis.good_vox'
            [SL.files(ii,:),~] = trimts(SL.files(ii,:),SL.analysis.sd_trim);
        end
    end
    % For safety, check for bad volumes as well - this has ended up being
    % specific to Hemi study
    keyboard;
    I1=mean(isinf(SL.files));
    I2=nanstd(SL.files);
    I3=mean(isnan(SL.files));
    I=find(I2==0);
    if ~isempty(I), 
        display('....Subject has bad volumes - skipping');
        continue;
    end
    %=====================================================================%
    %% Special Analysis Setup 
    %=====================================================================%
    keyboard;
    % Hemi Only
    if SL.Hemi==1,
        % Dead
    elseif SL.Hemi==2
       % If non-memory mode simply pulls information from behavorial files
       % for subsequent analysis use
       SL_hemi_mode(cursub);
       if ~strcmp(SL.Hemi_set.model,'Memory')
           % List sim
           R=zeros(length(SL.design.list));
           for ii=1:72
               I1=find(SL.design.list==ii);
               R(I1,I1)=1;
           end
           SL.design.matrix{2}=SL.design.matrix{2}.*R;
           clear R I1;
           
           % Special for semantic similarity
           dir1=fullfile(SL.root_dir,'fmri_reanalysis','hemi');
           load(fullfile(dir1,'SemanticCosineSim.mat')); % Loads R
           % Need the mapping for it as well
           A=excel_reader(fullfile(dir1,'hemi_stimlist_basic.csv'));
           Wlist=A{1}.col;
           % Now to determine the overlap for sorting

           c=1; Kl=[]; Im=[]; 
           for ii=1:length(SL.design.word_match)
               I=find(strcmp(SL.design.word_match{ii},Wlist));
               if ~isempty(I)
                   Im(ii)=I;
               else
                   Im(ii)=1;
                   Kl(c)=ii; c=c+1;
               end  
           end
           Rm=R(Im,Im); % Resorted R
           if ~isempty(Kl)
               Rm(Kl,Kl)=NaN;
           end
           R=Rm;
           SL.design.matrix{1}=SL.design.matrix{1}.*R;
       end     
    end

    %=====================================================================%
    %% ROI analysis
    %=====================================================================%
    % Hemi Analysis Area...
    % 1) Behave cross-check
    data{1}.col{cursub}=SL.dir.subjects{cursub};
    for ii=1:length(SL.design.anova)
        if cursub==1, data{ii+1}.header=SL.design.save_str{ii}; end
        data{ii+1}.col(cursub)=sum(sum(SL.design.anova{ii}.f{1}));
    end

    if cursub==length(SL.dir.subjects)
        write_struct(data,'D:\Data\fmri_reanalysis\hemi\young\trial_counts_anova.csv')
    end

    % 2) Iseries analysis
    % for each Iseries get data;
    %=====================================================================%
    %% Rankcorr analyses?
    %=====================================================================%
    if SL.geib==1

        display('Warning: special analysis settings need checked');
        keyboard;
        
        Nb=cumsum(SL.design.Box); N=Nb(end);
    sdir=fullfile(SL.dir.outpath,SL.dir.subjects{cursub});
    sdir_save=fullfile('D:\Data\fmri_reanalysis\hemi\Analysis\SL_5vox_sd2_10N_sfn_Iseries_rankcorr_sfn_v4',SL.dir.subjects{cursub});
    if ~exist(sdir_save,'dir'), mkdir(sdir_save); end
    f1=fullfile(sdir,'CR_Identity1_Iseries.img');
    f1t=fullfile(sdir,'CR_Identity1_Zmean.img');
    f2=fullfile(sdir,'FA_Identity1_Iseries.img');
    f2t=fullfile(sdir,'FA_Identity1_Zmean.img');
    f3=fullfile(sdir,'MI_Identity1_Iseries.img');
    f3t=fullfile(sdir,'MI_Identity1_Zmean.img');
    f4=fullfile(sdir,'HT_Identity1_Iseries.img');
    f4t=fullfile(sdir,'HT_Identity1_Zmean.img');
    
    vN1=SL.design.lang_norm(Nb(4)+1:Nb(5));
    vF1=SL.design.lang_freq(Nb(4)+1:Nb(5));
    vN2=SL.design.lang_norm(Nb(5)+1:Nb(6));
    vF2=SL.design.lang_freq(Nb(5)+1:Nb(6));
    vN3=SL.design.lang_norm(Nb(6)+1:Nb(7));
    vF3=SL.design.lang_freq(Nb(6)+1:Nb(7));
    vN4=SL.design.lang_norm(Nb(7)+1:Nb(8));
    vF4=SL.design.lang_freq(Nb(7)+1:Nb(8));
    
    f1s=fullfile(sdir_save,'CR_Identity1_Iseries_corr_rank.img');
    f1s2=fullfile(sdir_save,'CR_Identity1_Iseries_corr_freq.img');
    f2s=fullfile(sdir_save,'FA_Identity1_Iseries_corr_rank.img');
    f2s2=fullfile(sdir_save,'FA_Identity1_Iseries_corr_freq.img');
    f3s=fullfile(sdir_save,'MI_Identity1_Iseries_corr_rank.img');
    f3s2=fullfile(sdir_save,'MI_Identity1_Iseries_corr_freq.img');
    f4s=fullfile(sdir_save,'HT_Identity1_Iseries_corr_rank.img');
    f4s2=fullfile(sdir_save,'HT_Identity1_Iseries_corr_freq.img');
    
    f5s=fullfile(sdir_save,'OLD_Identity1_Iseries_corr_rank.img');
    f5s2=fullfile(sdir_save,'OLD_Identity1_Iseries_corr_freq.img');
    f6s=fullfile(sdir_save,'NEW_Identity1_Iseries_corr_rank.img');
    f6s2=fullfile(sdir_save,'NEW_Identity1_Iseries_corr_freq.img');
    f7s=fullfile(sdir_save,'ALL_Identity1_Iseries_corr_rank.img');
    f7s2=fullfile(sdir_save,'ALL_Identity1_Iseries_corr_freq.img');
    
    % CRs
    l1=load_nii(f1); d=size(l1.img);
    D1=reshape(l1.img,[prod(d(1:3)),size(l1.img,4)]);
    i1=isnan(vN1); i2=isnan(vF1); i3=isnan(sum(D1'));
    xN=corr(D1(~i3,~i1)',vN1(~i1)','type','spearman');
    xF=corr(D1(~i3,~i2)',vF1(~i2)','type','spearman');
    T=load_nii(f1t);   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f1s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f1s2);
    clear IN3img INimg T;
    
    % FAs
    l1=load_nii(f2); d=size(l1.img);
    D2=reshape(l1.img,[prod(d(1:3)),size(l1.img,4)]);
    i1=isnan(vN2); i2=isnan(vF2); i3=isnan(sum(D1'));
    xN=corr(D2(~i3,~i1)',vN2(~i1)','type','spearman');
    xF=corr(D2(~i3,~i2)',vF2(~i2)','type','spearman');
    T=load_nii(f2t);   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f2s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f2s2);
    clear IN3img INimg T;
    
    % MIs
    l1=load_nii(f3); d=size(l1.img);
    D3=reshape(l1.img,[prod(d(1:3)),size(l1.img,4)]);
    i1=isnan(vN3); i2=isnan(vF3); i3=isnan(sum(D1'));
    xN=corr(D3(~i3,~i1)',vN3(~i1)','type','spearman');
    xF=corr(D3(~i3,~i2)',vF3(~i2)','type','spearman');
    T=load_nii(f3t);   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f3s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f3s2);
    clear IN3img INimg T;
    
    % HTS
    l1=load_nii(f4); d=size(l1.img);
    D4=reshape(l1.img,[prod(d(1:3)),size(l1.img,4)]);
    i1=isnan(vN4); i2=isnan(vF4); i3=isnan(sum(D1'));
    xN=corr(D4(~i3,~i1)',vN4(~i1)','type','spearman');
    xF=corr(D4(~i3,~i2)',vF4(~i2)','type','spearman');
    T=load_nii(f4t);   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f4s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f4s2);
    clear IN3img INimg T;
    
    % OLD
    D=[D3 D4];
    vN=[vN3 vN4];
    vF=[vF3 vF4];
    i1=isnan(vN); i2=isnan(vF); i3=isnan(sum(D'));
    xN=corr(D(~i3,~i1)',vN(~i1)','type','spearman');
    xF=corr(D(~i3,~i2)',vF(~i2)','type','spearman');
    T=load_nii(f1t); % Any will do   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f5s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f5s2);
    clear IN3img INimg T;
    
    % NEW
    D=[D1 D2];
    vN=[vN1 vN2];
    vF=[vF1 vF2];
    i1=isnan(vN); i2=isnan(vF); i3=isnan(sum(D'));
    xN=corr(D(~i3,~i1)',vN(~i1)','type','spearman');
    xF=corr(D(~i3,~i2)',vF(~i2)','type','spearman');
    T=load_nii(f1t); % Any will do   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f6s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f6s2);
    clear IN3img INimg T;
    
     % ALL
    D=[D1 D2 D3 D4];
    vN=[vN1 vN2 vN3 vN4];
    vF=[vF1 vF2 vF3 vF4];
    i1=isnan(vN); i2=isnan(vF); i3=isnan(sum(D'));
    xN=corr(D(~i3,~i1)',vN(~i1)','type','spearman');
    xF=corr(D(~i3,~i2)',vF(~i2)','type','spearman');
    T=load_nii(f1t); % Any will do   
    INimg=NaN(length(i3),1); INimg(~i3)=xN;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f7s);
    clear IN3img INimg;
    INimg=NaN(length(i3),1); INimg(~i3)=xF;
    IN3img=reshape(INimg,[d(1:3)]);
    T.img=IN3img; save_nii(T,f7s2);
    clear IN3img INimg T;
    continue;
    end
%     
%     % Beware of NaNs
%     keyboard;
    %=====================================================================%
    fprintf(strcat('....Running searchlight on\t',SL.dir.subjects{cursub},'\n'));  
    if SL.region.noSL==1
    %=====================================================================%
    %% ROI analysis
    %=====================================================================%
        ROI_tmp=nan(length(mask_list),length(SL.design.matrix{1}));
        for jj=1:length(mask_list)
            SL.region.mask=mask_list{jj};
            X=spm_vol(SL.region.mask);
            if sum(sum(SL.V.mat-X.mat))~=0,
               sdisp('OH NO! Your mask doesn''t fit',1); 
               %error(['Everything is now broken: it is suggested you ' ...
                %   'reslice your mask and try again']);
            end
            display(['........Mask: ' mask_list{jj}]);
            % Some mask (especially in GT) lack full coverage. This
            % functions to skip those mask.
            [SL.LOC] = vect_sl_locations(SPM,SL);  
            if isempty(SL.LOC(1).voi),
                fprintf('..No coverage, skipping\n'); continue;
            end
            
            [~,mname,~]=fileparts(SL.region.mask);
            ROI_save=fullfile(SL.dir.outpath,'ROI',SL.dir.subjects{cursub},'data',[mname,'.mat']); 
            if ~exist(ROI_save,'file')
                SL.analysis.loop=randperm(length(SL.LOC));
                if FIRm.on==0, tic; searchlight(SL.dir.subjects{cursub},cursub); toc, end
            else
                continue;
            end
            %=============================================================%
            % FIR Pattern Conn Micro
            %=============================================================%
            if FIRm.on==1
                if FIRm.ROI==1
                    SL.LOC(1).box=[SL.LOC(:).voi]; % Box is all vois
                    SL.LOC(2:end)=[];
                    tmp = SL.files(SL.LOC(1).box,:);
                    tmp(isnan(nanmean(tmp,2)),:)=[]; % Remove NaN voxels
                    ROI_tmp(jj,:)=mean(tmp,1);   
                elseif FIRm.ROI==0
                    ROI_tmp=SL.files; break;
                end
            end
        end
        % This analysis is for all regions...
        if FIRm.on==1
            ER_FIR_models(ROI_tmp,cursub);
            clear ROI_tmp;
        end
    %=====================================================================%
    %% Searchlight analysis
    %=====================================================================%
    else % Natural searchlight zone
        [SL.LOC] = vect_sl_locations(SPM,SL);  
        SL.analysis.loop=randperm(length(SL.LOC));
        tic; searchlight_hemi_2(SL.dir.subjects{cursub},cursub); toc
    end
    
    if SL.Hemi==1, 
        SL.design.save_str={}; 
        SL.design.matrix={};
        SL.design.anova=[];
    end % Reset
    
    catch err
        smart_err(err); keyboard
        display([SL.dir.subjects{cursub} ': failed']);
    end
    % Exit if ERMatch analysis is run
    if SL.ERmatch==1, break; end
end % Subject loop

return;

