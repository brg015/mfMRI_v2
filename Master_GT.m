% Master_GT.m mfMRI_v1
% Matlab 2015a 
% BRG Winter 2016
%=========================================================================%
%% Initialize
%=========================================================================%
% Note: Many of the variables here are directly related to path names,
% really all this script is doing is setting up paths. At some point the
% structure will be more clearly written out (maybe)
clc; close all; clear all;
global SL;
global ntwk;
serv='serv1';
switch serv
    case 'local', root_dir='X:\';
    case 'serv1', root_dir='D:\Data\';
end
SL.root_dir=root_dir;
wrk_dir=fullfile(root_dir,'\fmri_reanalysis\hybridwords\');
% Specify where your function files are
code_dir{1}=fullfile(root_dir,'Geib\Scripts\Matlab_Scripts\mfMRI_v1');
code_dir{2}=fullfile(root_dir,'Geib\SPM\spm8\');
code_dir{3}=fullfile(root_dir,'Geib\SPM\xjview\');
code_dir{4}=fullfile(root_dir,'Geib\Scripts\Matlab_Scripts\function_files');
 for ii=1:length(code_dir), addpath(genpath(code_dir{ii})); end
%=========================================================================%
%% SL.dir
%=========================================================================%
%   .subjects (cell)  => list of subjects to include, I think this might
%                        only be used to define adjacency files, but it
%                        could be used as an index in other areas too
%   .save (string)    => output directory
SL.dir.subjects={'32158jdb' '32225jhs'};
SL.dir.save=fullfile(wrk_dir,'testing','GT_measures');
%=========================================================================%
% ntwk.analysis (cell array) => determines GT measures to compute
%   Eig           => eigenvariate centrality
%   C_Coef        => clustering coeffecient
%   Degree        => degree
%   wDegree       => weighted degree
%   Path_L        => path length
%   Global_Effc   => global efficiency
%   wMod1         => signed louvain modularity
%   Mod1          => louvain modularity
%   Local_Effc    => local efficiency
%   BW_cent       => betweeness centrality
%   PR_cent       => page range centrality
% * These are all caled line 28-69 in GT_compare_measures.m

% ntwk.atlas (str)           => determines the network used
% ntwk.contrast.
%   pos (cell array)         => positive network names
%   neg (cell array)         => negative network names
%   nam (cell array)         => name of contrast
%   pos_R (cell array)       => set of positive adjacencies to include
%   neg_R (cell array)       => set of negative adjacencies to include
% ntwk.FSR (logic)           => Computes FSR 

ntwk.analysis={'Global_Effc'};
ntwk.atlas='AAL_negX';
ntwk.FSR=1;

ntwk.contrast.pos={'BetaHit'};
ntwk.contrast.neg={'BetaMiss'};
ntwk.contrast.nam={'BetaHvM'};
for ii=1:length(SL.dir.subjects)
    ntwk.contrast.pos_R{1}{ii}=fullfile(wrk_dir,'testing','ROI','Group',...
        'AAL_negX','beta','RetHit',[SL.dir.subjects{ii} '.mat']);
     ntwk.contrast.neg_R{1}{ii}=fullfile(wrk_dir,'testing','ROI','Group',...
        'AAL_negX','beta','RetMiss',[SL.dir.subjects{ii} '.mat']);
end
%*note that multiple sets of analyses can be run at once...
%=========================================================================%
% SL.transform     => beta exponent transform of adjacency matrix 0 = off
% SL.write.gephi   => generate output for gephi 0 = off
% SL.write.compute => flag indicating that GT measures need computed
SL.transform=0;
SL.write.gephi=0;
SL.write.compute=1;
%=========================================================================%
% ONLY pertinent of gephi/subnetwork architecture is computed...
ntwk.gephi.thresh=2.58; % N=17
ntwk.gephi.sign='T';
ntwk.gephi.outZ=3.5;
ntwk.gephi.save='D:\Data\fmri_reanalysis\hybridwords\gephi\alpha01_z35_0815_manuscript_beta_permbin\';
% ntwk.gephi.save=fullfile(root_dir,'fmri_reanalysis\hybridwords\gephi\alpha01_z35_0815_manuscript_ppi\');
ntwk.gephi.node={39};
ntwk.gephi.extra.file='D:\Data\fmri_reanalysis\hybridwords\PPI\gephi_network_values.csv';
ntwk.gephi.extra.col=[10 11 12 13];

GT_shell_v3;
%=========================================================================%
%% Output
%=========================================================================%
% Folders are created for the positive, negative, and difference contrast
%  Positive/Negative folders have
%   .mat files for the output data
%   _avg.csv files for average
%   _avg.mat files for average
%   _avg.img files mapped back to brain space
%  nam (difference folders) have
%   _ranksum_[p or z][.csv,.img,.mat] with ranksum stats
%   _tstats_full with means, SEs, T-values for the conditions
%   _Ttest_[p or tval][.csv,.img,.mat] with all t-test stats
%   _permutation[.csv,.mat,.img] with permutation stats
      


