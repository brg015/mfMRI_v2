function RSA_shell_v2()
%=========================================================================%
%% Code Setup
%=========================================================================%
% Updates 11/12/2014
% 1) Added many presets to be internal, for example, all of SL.analysis and
% all of SL.preprocess. Plus a few others, basically, simplified the Master
% script by making more internal assumptions
% Updates 4/9/2017
% 1) Removed FIRm

global SL;
% Setting fields if they do not exist and display warnings
bugger=0;
sdisp('Warning: Regression turned off',1); pause(2);
sdisp('Warning: Beta version',1); pause(2);
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
% SL.run.include == 1 -> only includes comparisons between runs
% SL.run.include == 2 -> only includes comparisons within runs

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

NModels=length(SL.design.save_str);

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
if ~isfield(SL.analysis,'sd_trim'),         SL.analysis.sd_trim=2.5; end
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
    %% Setup design matrix space - clear all subject variables
    %=====================================================================%
    % Set flag and subject specific parameters
    flag=0;
    SL.design.ID_idx=[];
    SL.design.ID_descrip={};
    SL.design.ID_file={};
    SL.design.Box=[];
    SL.files=[];
    SL.design.f=[];
    SL.design.matrix=[];
    SL.run.matrix=[];
    SL.design.run=[];
    %=====================================================================%
    %% Load in subject SPM.mat file
    %=====================================================================%
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
    %% Remove Outliers & Setup searchlight
    %=====================================================================%
    fprintf(strcat('....Creating\t',SL.dir.subjects{cursub},'''s Search Light indices\n'));
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
    %% Arrange design matrices
    %=====================================================================% 
    % Define
    % SL.design.matrix => design matrix
    % SL.design.power  => matrix intensity
    for ii=1:NModels
        SL.design.matrix{ii}=[];   % Initialize model
        SL.design.power{ii}=[];
        RSA_design_matrix_v3(ii); 
    end
    if flag==1, 
        fprintf([SL.dir.subjects{cursub} ' FAILED\n']);
        if bugger==1,  else continue; end
    end
    %=====================================================================%
    fprintf(strcat('....Running searchlight on\t',SL.dir.subjects{cursub},'\n'));  
    if SL.region.use_mask==2
        %=================================================================%
        % Voxel Networks
        %=================================================================%
        V = spm_vol(SL.design.ID_file{1});        
        a1=cumsum(SL.design.Box);
        keyboard;
        for ii=1:length(SL.design.Box)
            b1=SL.files(SL.analysis.good_vox,a1(ii)-a1(1)+1:a1(ii));
            R=single(corr(b1')); data.R=R; data.V=V;
            data.vox=SL.analysis.good_vox; data.R=int16(data.R*1000);
            save(fullfile(SL.dir.outpath,[SL.dir.subjects{ii} '_' SL.desgn.cond_names{ii} '.mat']),'data','-v7.3');
            clear data R;
        end
        clear a1 b1;
       continue; 
    end
    if SL.region.noSL==1
    %=====================================================================%
    %% ROI analysis
    %=====================================================================%
        % ROI_tmp=nan(length(mask_list),length(SL.design.matrix{1}));

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
                tic; searchlight(SL.dir.subjects{cursub},cursub); toc
            else
                continue;
            end
            SL.dir.check=0; % Only run this for first iteration
            
        end
    %=====================================================================%
    %% Searchlight analysis
    %=====================================================================%
    else % Natural searchlight zone
        [SL.LOC] = vect_sl_locations(SPM,SL);  
        SL.analysis.loop=randperm(length(SL.LOC));
        tic; searchlight(SL.dir.subjects{cursub},cursub); toc
    end
    
    catch err
        smart_err(err); 
        display([SL.dir.subjects{cursub} ': failed']);
    end
end % Subject loop

