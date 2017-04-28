function RSA_shell()
%=========================================================================%
%% Code Setup
%=========================================================================%
% Updates 11/12/2014
% 1) Added many presets to be internal, for example, all of SL.analysis and
% all of SL.preprocess. Plus a few others, basically, simplified the Master
% script by making more internal assumptions
global SL FIRm;

% Setting fields if they do not exist and display warnings
bugger=0;
sdisp('Warning: Regression turned off',1); pause(2);
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
if ~isfield(SL.analysis,'sd_trim'),         SL.analysis.sd_trim=2; end
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
    [SL.files] = load_st_betas(SL);
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
    % Between Subject Similarity Analysis
    if SL.ERmatch==1, 
%         if exist(fullfile(SL.dir.outpath,'N21_dat.mat'),'file')
            X=load(fullfile('D:\Data\Wing\ERMatch_serv1\Analysis\STwa\SL_ER_vols_5vox_sd2_10N\','N21_dat.mat'));
            SL.files=X.SL.files;
            SL.design=X.SL.design;
            clear X; 
%         else
%             if cursub==1, S={}; end
%             S=ER_ss_analysis(S,cursub); 
%             if cursub<length(SL.dir.subjects), continue; end
%             save(fullfile(SL.dir.outpath,'N21_dat.mat'),'SL','-v7.3'); 
%         end
    end
    % Hemi Only
    if SL.Hemi==1,
        SL.design.calc{1}='Anova1';
        SL.design.save_str{1}='SemanticSim_runfix';
        SL.design.anova{1}.type='Anova1';
        SL.design.anova{1}.cond=[1 2];
        SL.design.anova{1}.names={'On' 'Off'};
        v=[]; for ii=1:length(SL.design.Box)
            eval(['dv' num2str(ii) '=ones('...
                num2str(SL.design.Box(ii)) ',' num2str(SL.design.Box(ii)) ');']);
            v=[v [',dv' num2str(ii)]];
        end
        v(1)=[];
        eval(['X=blkdiag(' v ');']);
        if length(X)~=length(SL.design.ID_idx),
            fprintf('Bad Subject Broke\n');
            continue;
        end

        SL.design.anova{1}.f{1}=X;
        Y=X; 
        Y(X==0)=1; 
        Y(X==1)=0;
        SL.design.anova{1}.f{2}=Y;
        
        % NOTE: This has been updated to remove bw run comparisons... the
        % section of code provided below is what allowed this
        D=cumsum(SL.design.Box);
        for aa=1:length(SL.design.Box)
            for bb=1:length(SL.design.Box)
                if aa~=bb
                   if aa==1, I1=1:D(aa); else I1=D(aa-1)+1:D(aa); end
                   if bb==1, I2=1:D(bb); else I2=D(bb-1)+1:D(bb); end
                   in=~eye(length(I1),length(I2));
                   SL.design.anova{1}.f{2}(I1,I2)=in;
                end
            end
        end
        clear X Y;
        
        indx=1;
        SL.design.matrix{indx}=ones(sum(SL.design.Box),sum(SL.design.Box));
        SL.design.matrix{indx}=tril(SL.design.matrix{indx},-1);
        nan_idx=tril(ones(length(SL.design.matrix{indx})),-1)==0;
        SL.design.matrix{indx}(nan_idx)=nan;
        if strcmp(SL.design.calc{indx},'Anova1') || strcmp(SL.design.calc{indx},'Anova2') ...
            || strcmp(SL.design.calc{indx},'Identity1');
            if isfield(SL.design.anova{indx},'f')
                for jj=1:length(SL.design.anova{indx}.f)
                    SL.design.anova{indx}.f{jj}(nan_idx)=0; % Anova map, so ==0 
                    if std(std(SL.design.anova{indx}.f{jj}))==0
                        SL.design.anova{indx}.f{jj}=[];
                        display('Warning: Anova Matrices have no variability, removing');
                    end
                end
                SL.design.anova{indx}.f=SL.design.anova{indx}.f(~cellfun('isempty',SL.design.anova{indx}.f));
            end
        end
    elseif SL.Hemi==2
       SL_hemi_mode(cursub);
    end
    
    % GBU setup
    if SL.GBU==1
        % Setup spearman models
        SL.design.calc{8}='Spear';
        SL.design.save_str{8}='SpearModel';
        SL.design.custom(8)=0;
        SL.design.matrix{8}=zeros(length(SL.design.ID_idx),length(SL.design.ID_idx));
        D=cumsum(SL.design.Box);
        for aa=1:length(SL.design.Box)
            for bb=1:length(SL.design.Box)
                iv=8-abs(aa-bb);
                % Gotta get cords
                if aa==1, r=1:D(1); else r=D(aa-1)+1:D(aa); end
                if bb==1, c=1:D(1); else c=D(bb-1)+1:D(bb); end
                SL.design.matrix{8}(r,c)=iv;
                clear iv r c;
            end
        end
        clear D;
    end
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
               error(['Everything is now broken: it is suggested you ' ...
                  'reslice your mask and try again']);
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
        % Needs forced to overcome saved N21_dat.mat file...
        if SL.ERmatch==1, SL.design.SSL=3; end
        [SL.LOC] = vect_sl_locations(SPM,SL);  
        SL.analysis.loop=randperm(length(SL.LOC));
        tic; searchlight(SL.dir.subjects{cursub},cursub); toc
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

if FIRm.on==1
    % Setup save strings...
    switch FIRm.mem
        case 1, s1='cnt';
        case 2, s1='HL';
    end
    switch FIRm.model
        case 'CSA'
            Sdir=fullfile(SL.dir.outpath,'Connectivity');
        case 'TSA'
            Sdir=fullfile(SL.dir.outpath,'Temporal');
    end
    M=mean(FIRm.Region_E1,2);
    S=std(FIRm.Region_E1');
    T=M'./(S/sqrt(21)); 
    save_asA=fullfile(Sdir,['ER_' s1 '.nii']);
    [~]=ROI_constructor(mask_list,T,[],'pure',save_asA,3);
    
    M=mean(FIRm.Region_E2,2);
    S=std(FIRm.Region_E2');
    T=M'./(S/sqrt(21)); 
    save_asB=fullfile(Sdir,['ME_' s1 '.nii']);
    [~]=ROI_constructor(mask_list,T,[],'pure',save_asB,3);
end
