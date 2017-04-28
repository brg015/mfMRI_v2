% Master_Network_Maker.m mfMRI_v1
% Matlab 2015a 
% BRG Winter 2016
% 
clc; clear all; close all;
global SL;
%=========================================================================%
%% Add paths
%=========================================================================%
OP='serv1';
switch OP
    case 'local', SL.dir.root='Z:\';
    case 'serv1', SL.dir.root='D:\Data\';
end
wrk_dir=fullfile(SL.dir.root,'\fmri_reanalysis\hybridwords\');
code_dir{1}=fullfile(SL.dir.root,'Geib\Scripts\Matlab_Scripts\mfMRI_v1');
code_dir{2}=fullfile(SL.dir.root,'Geib\Scripts\Matlab_Scripts\function_files');
for ii=1:length(code_dir), addpath(genpath(code_dir{ii})); end
%=========================================================================%
%% SL.dir
%=========================================================================%
%   .subjects (cell)  => should be the same as from Master_RSA
%   .outpath (string) => is SL.dir.outpath from Master_RSA + ROI
SL.dir.subjects={'32158jdb' '32225jhs'}; 
SL.dir.outpath=fullfile(wrk_dir,'testing','ROI');
%=========================================================================%
%% SL.ROI
%=========================================================================%
%   .prefix (string) => [SL.dir.QA '_'], if a QA dir wasn't set, then this
%                       can be left blank e.g. ''
% HYBRID is a neg rotation*
SL.ROI.set='AAL_negX';
SL.ROI.prefix='Int_';
%=========================================================================%
%% SL.analysis.ROI
%=========================================================================%
%   .save (cell str)   => Description of the type of networks you're
%                         interested in creating
%   .Ibox (matrix)     => binary matrix. Each row is representative of a
%                         single network, 1 indicates the inclusion of a
%                         condition, 0 indicates its not included.
%                         Conditions are defined within the RSA script.
%   .overwrite (bin)   => Option to overwrite existing data
SL.analysis.ROI.save={'RetHit' 'RetMiss'};
SL.analysis.ROI.Ibox=[0 0 1 0;...
    0 0 0 1];
SL.analysis.ROI.overwrite=1;
%=========================================================================%
ROI_shell_v2;

%=========================================================================%
%% Output Structure
%=========================================================================%
% Files are saved to [SL.dir.outpath 'Group' SL.ROI.set]
%  Subfolders exist under beta for SL.analysis.ROI.save conditions
%   Each of those folders has a file for every subject
%    These files contain variables R and R1
%     R => original connectivity matrix
%     R1=> connectivity matrix in which negatives and diagonal have been 0


























