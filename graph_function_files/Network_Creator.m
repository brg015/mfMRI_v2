% Function file
% BR Geib 2012
% Designed from beta_series_correlation_nomars by Dennis Thompson
function  GT_Network_Creator...
    (SPM_loc,ROI_loc, Events,trimsd,sav_dir,correlation_method)
% beta_series(SPM_loc, ROI_loc, Events,trimsd,sav_dir)
%
% takes the beta series from one roi (the seed) and correlates it to
% all the voxels in the brain and saves the results as an image
% SPM_loc - a string that points to a single SPM mat file location
% results will be written in to the same location as the SPM mat file
%
% ROI_loc - a string that points to a single nifti format ROI file location
%
% Events - A cell array of strings that will be used to identify the Events in the beta series
% This allows for a multipass string isolation
% For example Events = {'GreenProbe_Correct', 'Sn(1)','bf(1)'}
% First locate the set of events that match 'GreenProbe_Correct',
% Second locate the set of events that match Sn(1)
% Third locate the set of events that match bf(1)
% The final set of event will be the those in the intersection of the
% three sets
% Events can be a single string
%
% trimsd is the is number of standard deviations to use if
% you wish the raw data to be Windsorized
% set to 0 if you wish the raw data to be used
%
% sav_dir just dictates where the output should be written to

if ~exist('trimsd','var'), trimsd = 0; end
if ~iscell(Events), Events = {Events}; end

if ~exist(sav_dir,'dir'), mkdir(sav_dir); end

% load SPM mat file
load(SPM_loc);

% locate beta_images related to Event
P = location_of_beta_images_from_event_discription(Events,SPM);
% get beta values (in Vector form) 
vbetas = spm_get_data(P,SPM.xVol.XYZ);
% get ROI
roi = spm_get_data(ROI_loc,SPM.xVol.XYZ);
idx = find(roi);
% extract mean of ROI from each beta
for n = 1:size(vbetas,1)
    mean_roi(n) = mean(vbetas(n,idx));
end


if trimsd > 0,
    mean_roi = trimts(mean_roi, trimsd, []);
end
% for each voxel in beta images
for n = 1:size(vbetas,2),
    if trimsd > 0,
        vbetas(:,n) = trimts(vbetas(:,n), trimsd, []);
    end
    Cout(n)  = corr(mean_roi',vbetas(:,n),'type','Pearson');
end


[~,roiLabel,~] = fileparts(ROI_loc);
% output R correlation results to image

switch correlation_method
    case 'Pearson'
        corr_file = fullfile(sav_dir,['Rcorr_',roiLabel,'_' Events{1},'.nii']);
        writeCorrelationImage(Cout,corr_file, SPM.xVol);
    case 'FischerR'
        % output R correlation results to image
        corr_file = fullfile(sav_dir,['R_atanh_corr_',roiLabel,'_' Events{1},'.nii']);
        writeCorrelationImage(atanh(Cout),corr_file, SPM.xVol);
    case 'FischerZ'
        % output Z correlation results to image
        corr_file = fullfile(sav_dir,['Zcorr_',roiLabel,'_' Events{1},'.nii']);
        writeCorrelationImage((atanh(Cout)*sqrt(length(P)-3)),corr_file, SPM.xVol);
end





