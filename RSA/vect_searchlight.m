% Vect_searchlight
% Master Level Script
% Karl - BRG MVPA edits 2014 (Winter)
%
clc; clear all;
% Description: Overlap exists between the mvpa settings and SL settings.
% This will not be updated as mvpa is being included as a simple addition
% to SL possiblities. See variables for more in-depth description.
%=========================================================================%
% General Inputs:
%=========================================================================%
SL.stpath = 'J:\ERMatch_Sol\Analysis\STwa\';
%'\\ccn-cabeza34.win.duke.edu\C$\fmristudies\ERMatch\Analysis\SingleTrial_ER';
SL.outpath =  'J:\ERMatch_Sol\Analysis\STwa\SL_vols_5vox\2ndOrdERS\Viv4_MVPA\';
%'\\ccn-cabeza34.win.duke.edu\C$\fmristudies\ERMatch\Analysis\SingleTrial_ER\SL_vols\sub';
SL.subjects =  {'13549' '13552' '13562' '13617' '13655' '13658' '13683'};
%=========================================================================%
% Searchlight Input Variables:
%=========================================================================%

SL.interest = {'Viv4'}; % CAN ONLY USE ONE CONDITION RIGHT NOW!
%{'FCDM1'}
SL.name = 'full';  %'cnd1sac1';
%filter for first encoding trial (not encoding 2 or retrieval) uses
%extracted ID to find pair
SL.b_filter = ''; % binary behav
%'RATING4'
SL.id_start = 'ID';
SL.id_end   = '*bf'; 
SL.hm = 0; %1 if you want hits, 0 for misses - binary from behav
SL.model = 'shuffle_incr'; %'shuffle' if you want random (shuffle) condition, 'match' for id match, 'shuffle_incr' calc similarity for every id with every other. 'match_zdistrib' for z-scored ID sim values, 
SL.SSL = 5; % size of search light: 3x3x3 USE ODD NUMBERS

if ~exist(SL.outpath,'dir')
    mkdir(SL.outpath);
end
if SL.hm == 1
    SL.name = strcat(SL.name,'_hit');
elseif SL.hm == 0 && isempty(SL.b_filter)
    SL.hm = 1;
    SL.name = strcat(SL.name,'_all');
else
    SL.name = strcat(SL.name,'_miss');
end
addpath(fileparts(which('vect_searchlight')));

%=========================================================================%
% MVPA variables Input Variables:
%=========================================================================%
% SL.model='mvpa'; 
SL.mvpa.on=1;
SL.mvpa.algorithm='defualt'; 

%=========================================================================%
%% Code Start
%=========================================================================%
fprintf(strcat('VECTORIZED SEARCHLIGHT\t',datestr(clock),'\n'));
for cursub = 1:length(SL.subjects)
    
    cd(strcat(SL.stpath,'\',SL.subjects{cursub}));
    fprintf(strcat('....Loading\t',SL.subjects{cursub},'''s SPM.mat\n'));
    load SPM;
    
    % Loading in voxel information from all beta files
    SL.V=spm_vol(SPM.Vbeta(1).fname);
    fprintf(strcat('....Loading\t',SL.subjects{cursub},'''s beta images\n'));
    [SL.files] = load_st_betas(SPM);
keyboard
    % Returns SL.conds - which should be an index of any Vbeta area that 
    % compares with SL.interest
    fprintf(strcat('....Extracting\t',SL.subjects{cursub},'''s condition vector\n'));
    [SL.conds] = make_cond_mat(SPM,SL);
    
    fprintf(strcat('....Extracting\t',SL.subjects{cursub},'''s behavior\n'));
    [SL.behav] = make_behav_mat(SPM,SL);
    
    SL.indx = collapse(SL.conds);    
    SL.btod = find(SL.indx == 0);
    
    SL.files(:,SL.btod) = [];
    SL.indx(:, SL.btod) = [];
    SL.behav(:,SL.btod) = [];
    
    fprintf(strcat('....Creating\t',SL.subjects{cursub},'''s Search Light indices\n'));
    [SL.LOC] = vect_sl_locations(SPM,SL);
    
    fprintf(strcat('....Running searchlight on\t',SL.subjects{cursub},'\n'));
    [SL.out] = searchlight(SPM,SL,SL.subjects{cursub});

    clear files conds behav btod out mean_out
end % Subject Loop