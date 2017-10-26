% Master_RSA.m mfMRI_v1
% Matlab 2015a 
% Winter 2016
% Karl (Searchlight) & BRG updates
%
% Description
%   Searchlight script. Performs searchlight style analysis of any
%   arbitrary function. The searchlight volume is a 'cube' by default.
%   Searchlight volumes can also be defined as ROIs e.g. multivariate
%   analyses can be run discretely. Finally, the script is designed to
%   function with Master_Network_Maker.m. It does so by saving beta-series
%   with respect to an ROI. 
%
% clear existing variables and declare SL
clc; clear all;
global SL;

% Update log
% 3/28/2016  -> added comments and fixed custom models
% 6/23/2017  -> code for distance computations and SVM classifier
%            * still needs better documentation
% 10/26/2017 -> added code for multiple RDMs
%            * QA doesn't work with this yet             
%=========================================================================%
%% Add paths
%=========================================================================%
% mfMRI_v2, function_files, and spm8 need to be added
% xjview is not essential, but my preference
OP='serv1';
switch OP
    case 'local', root_dir='Z:\';
    case 'serv1', root_dir='D:\Data\';
end
wrk_dir=fullfile(root_dir,'\Wing\ERMatch_serv1\');
code_dir{1}=fullfile(wrk_dir,'\MultipleRDMs\mfMRI_v2-master\');
code_dir{2}=fullfile(root_dir,'Geib\Scripts\Public\SPM\spm8\');
code_dir{3}=fullfile(root_dir,'Geib\Scripts\Public\SPM\xjview\');
code_dir{4}=fullfile(root_dir,'Geib\Scripts\Public\function_files');
for ii=1:length(code_dir)
    addpath(genpath(code_dir{ii}));
end
%=========================================================================%
%% SL.dir
%=========================================================================%
% Variables:
%   .stpath (string)     => Location of single trials to use
%   .subjects (cell)     => Each element is a subject ID (folder)
%   .outpath (string)    => Where to save the analysis
%   .overwrite (logical) => Determines if files should be overwritten, if
%                           is set to 0, will skip over any existing
%                           output (still testing...)
%   .QA (string)         => Quality assurance directory
SL.dir.stpath = fullfile(wrk_dir,'STwa\');
SL.dir.subjects={'13100' '13205' '13220' '13448' '13491' ...
    '13510' '13523' '13549' '13552' '13562' ...
    '13617' '13655' '13658' '13683' '13693' ...
    '13720' '13743' '13779' '13793' '13807' ...
    '13845'};
SL.dir.outpath =  fullfile(wrk_dir,'\testing\');
SL.dir.overwrite=1;

SL.dir.QA='QA'; 
SL.dir.QAcheck=0;

%*Expected directory of 1st SPM.mat file is...
%  SPM=fullfile(SL.dir.stpath,SL.dir.subjects{1},'SPM.mat');
%=========================================================================%
%% SL.design & SL.ID & SL.run
%=========================================================================%
% Varaibles
%   .mumford (binary)       => if ==1 indicates that a mumford style
%                              analysis was run. Meaning that beta_xxxx.img
%                              files do not exist and that file names must
%                              be processed as opposed to SPM.mat files
%   .mumford_dir (string)   => subject sub-directory indicating where the
%                              brain volumes have been saved
%   .mumford_str (string)   => wildcard string indicating the file type of
%                              the saved volumes
SL.design.mumford=0;
SL.design.mumford_dir='betas';
SL.design.mumford_str='img';
%*Expected directory of volumes if mumford==1
% Files=fullfile(SL.dir.stpath,SL.dir.subjects{1},SL.design.mumford_dir,...
%  ['*' SL.design.mumford_str]);
% There is direct link between the file description in the SPM.mat file and
% the beta volumes. When Mumford is turned on, it skips the SPM.mat
% reference and directly analyzes the file names
%
% Variables
%   .SSL (num)              => Size of the searchlight (must be odd)
SL.design.SSL=5; % size of search light: 3x3x3 USE ODD NUMBERS
%   .cond_str (cell)        => Conditions to include in RDM
%       => 1* Order determines order of sorted matrix
%       => 2* Order determined by ID match, if exist (should)
%       => Requires nicely defined (descriptive) Vbeta fields, works by
%       performing an & conjunction on any strings within a cell within a
%       Vbeta.descrip string. If all strings found, then included in said
%       condition. Overlaps must never occur. A single '|' statement is
%       permitted
SL.design.cond_str={...
    {'Phase1','Run','Viv1' '|' 'Viv2'},...
    {'Phase1','Run','Viv3' '|' 'Viv4'},...
    {'Phase2','Run','Viv1' '|' 'Viv2'},...
    {'Phase2','Run','Viv3' '|' 'Viv4'}};
%   .check (num)      => Determines if IDs should be checked, if this ==0
%                        then IDs are not checked
%   .match (cell x2)  => Which strings to look between for the ID, the
%                        string being examined is in Vbeta
%   .overlap{x}(array)=> Determines number of overlaps to look for. An
%                        overlap array consists of zeros with a 1 and 2
%                        used to indicate which conditions (deterimined by
%                        SL.design.cond_str) should have aligned IDs 
%   .setcheck (array) => Determines which type of overlap to look for, one
%                        condition must exist for each overlap defined
%                           0 : 1 == 2
%                           1 : maximum overlap between 1 and 2
SL.ID.check=1;
SL.ID.match={'ID' '*bf'}; 
SL.ID.overlap{1}=[1 0 2 0]; % Overlap set1
SL.ID.overlap{2}=[0 1 0 2]; % Overlap set2
SL.ID.setcheck=[0 0];
% Variables
%   SL.run.include (bin)  => determines if within run items should be
%                            excluded from the analysis, if ==1, then this
%                            is so
%   SL.run.match {cell x2}=> Same as SL.ID.match, it determines the strings
%                            to look between in order to determine what run
%                            each beta belongs to
%   * in the below formulation, it suggests we should not include items
%   within runs, and that the run number is found in vbeta
%   '...Run1_Tri60...' e.g. run #1 in this case
SL.run.include=1;
SL.run.match={'Run' 'Tri'};
%=========================================================================%
%% SL.region
%=========================================================================%
%   .use_mask (num)   => Determines if a mask should be used 
%   .mask (cell)      => Exact location of the mask to use - if MVPA is
%                        being done, then this must be set to a valid mask.
%                        As it is used in the initialization of the
%                        toolbox.
%   .noSL (num)       => If set ==1 then searchlight is turned OFF. This is
%                        true of MVPA or RSA analyses
%   .noModel (logical)=> ==0 if you want to run no models, this is
%                        typically true if you're just trying to save beta
%                        series or such
SL.region.use_mask=0;
SL.region.noSL=0;
SL.region.noModel=0;

mask_dir=fullfile(root_dir,'\Geib\WFU_ROIs_resliced_f_negX\');
X=dir([mask_dir '*.nii']);
for ii=1:length(X)
    SL.region.mask{ii}=fullfile(mask_dir,X(ii).name);
end
%*Typically use_mask and no_SL are set equal. For example, while you could
% run a searchlight within a mask, this seems a bit odd
%=========================================================================%
%% SL
%=========================================================================%
%   .fisher (logical)       => controls fisher transformation correlation
%                              values, by defualt this is set ==0
%   .analysis.sd_trim (num) => number of standard deviations to windsorize
%                              the time series by
%=========================================================================%
%% SL.design
%=========================================================================%
%*Not needed for beta series analysis
%   .calc{#} (string)  => Multiple types of calculations are possible, the
%                         most common is 'Identity1'
%       Identity1 => examines the difference between identical items and
%                    mismatched items
%           .interactions{X}=[a b; c d]   => conditions to include the model in a row column format.
%           .save_str{#} (string)         => name of model for output
%           .custom(#) (logical)          => 0
%           .Identity1{#}.regress.on (logical)=> 0 (in testing atm)
%           .Identity1{#}.type (string)       => Identity1
%           .Identity1{#}.names (cell)        => Names of on/off diag
%           .Identity1{#}.diag=[a,b]          => defines the on diag cells
%           .Identity1{#}.row (logical)       => zscore rows (2=both,1=row,0=col)
%           .Identity1{#}.FourD (logical)     => save 4D output
%       Spear, Euclid, Mean or Kendall => for custom models only
%           .interactions{X}              => same as above
%           .save_str{#}                  => same as above
%           .custom(#)=1                  => says to include custom model
%       MRegression => for custom models only
%           .interactions{X}              => same as above
%           .save_str{#}                  => same as above
%           .custom(#)=1                  => says to include custom model
%           .ortho(#)=1                   => says to orthogonalize to first model    
%   .model{X} (file string)               => indicates data file to load
%       *loaded data file should be a .mat with 2 variables R &
%       stim_ID_num. R is a cross correlation matrix. stim_ID_num maps the
%       rows/columns of the matrix to SL.ID.match s.t. numbers pulled from
%       the Vbeta descriptions can match matrix R
%       *only needs set if this is an Spear, Euclid or Kendall analysis
%           design.interactions{X} => cells to compute average over
%           design.save_str{X}     => name of the output volume
%           custom(X)              => ==1 if design is a custom matrix
%       *if MRegression then this should be a cell list of models
%=========================================================================%
%% Examples of all three models
%=========================================================================%
model=1;
% SL.design.calc{model}='Identity1';
% SL.design.interactions{model}=[3 1; 3 2; 4 1; 4 2];
% SL.design.save_str{model}=['ERS'];
% SL.design.custom(model)=0;
% SL.design.Identity1{model}.regress.on=0;
% SL.design.Identity1{model}.type='Identity1';
% SL.design.Identity1{model}.names={'On' 'Off'};
% SL.design.Identity1{model}.diag=[3 1; 4 2]; %index to cross-phase ID diagonal containing block cells 
% SL.design.Identity1{model}.row=1;
% SL.design.Identity1{model}.FourD=1;
% model=model+1;
% 
% SL.design.calc{model}='Anova1';
% SL.design.interactions{model}=[1 1; 2 1; 2 2];
% SL.design.save_str{model}='Encoding';
% SL.design.custom(model)=0;
% model=model+1;
% 
SL.design.calc{model}='Spear';
SL.design.model{model}=['D:\Data\Wing\ERMatch_serv1\RSA_models\features\model.mat'];
SL.design.interactions{model}=combn(1:4,2);
SL.design.save_str{model}='Enc_Spear';
SL.design.custom(model)=1;
model=model+1; 

SL.design.calc{model}='MRegression';
SL.design.model{model}={'D:\Data\Wing\ERMatch_serv1\RSA_models\VGG16_Places_mat_files\model_layer_nan_3.mat',...
    'D:\Data\Wing\ERMatch_serv1\RSA_models\VGG16_Places_mat_files\model_layer_nan_11.mat',...
    'D:\Data\Wing\ERMatch_serv1\RSA_models\VGG16_Places_mat_files\model_layer_nan_21.mat',...
    'D:\Data\Wing\ERMatch_serv1\RSA_models\VGG16_Places_mat_files\model_layer_nan_26.mat'};
SL.design.interactions{model}=combn(1:4,2);
SL.design.save_str{model}='Enc_VGG_Spear';
SL.design.ortho(model)=1; % This doesn't work
SL.design.custom(model)=1;
model=model+1; 

% Beta phase models...
SL.design.calc{model}='distance';
SL.design.distance{model}='euclidean';
SL.design.save_str{model}='euclidean';
SL.design.class{model}=[1 2];
SL.design.custom(model)=0;
model=model+1;

SL.design.calc{model}='MVPA';
SL.design.save_str{model}='MVPA';
SL.design.class{model}=[1 2];
SL.design.custom(model)=0;
model=model+1;

RSA_shell_v2;
%=========================================================================%
%% Output Structure
%=========================================================================%
% ROI (beta-series) analysis
%  Data is saved to SL.dir.outpath. Each subject is given a single folder
%  in which the folder 'data' is present. The data folder will have one
%  file for each ROI. Files will named [SL.dir.QA '_' SL.region.mask{#}].
%  
%  Data files have three variables 'R', 'tmp' and 'design'
%   R      => cross-correlation matrix from the given region
%   tmp    => data matrix from ROI i.e. (voxel X beta)
%   design => Information on the general design structure
%
% QA directory (SL.dir.QA)
%  This folder will be located in SL.dir.outpath and named SL.dir.QA. A
%  single folder will exist for each model run. These folders will have
%  contrast maps for each condition. The folder will also have .csv files
%  detailing which beta files have been included in each condition.























