% Function file
% BR Geib 2012
% Designed from beta_series_correlation_nomars by Dennis Thompson
function  idx=GT_Network_Creator...
    (SPM_loc,ROI_dir,ROI_list,Events,sav_dir)

if exist(SPM_loc,'file')
	load(SPM_loc),
	[root,~,~]=fileparts(SPM_loc);
	SPM.swd=root; % Ensure correct path
else
	fprintf('Missing SPM.mat file - exitting\n');
end

if ~exist('trimsd','var'), trimsd = 0; end
if ~iscell(Events), Events = {Events}; end
threshold=0;
covg_limit=0.75;

if ~exist(sav_dir,'dir'), mkdir(sav_dir); end

% locate beta_images related to Event
P = location_of_beta_images_from_event_discription(Events,SPM);
% Get header info for beta data
V = spm_vol(P);

idx=[]; c=1;
for i=1:length(ROI_list)
    % get ROI index and transform matrix
    [XYZ ROImat]= roi_find_index([ROI_dir ROI_list(i).name],threshold);
    % generate XYZ locations for each beta image
    % correcting for alignment issues
    betaXYZ = adjust_XYZ(XYZ, ROImat, V);

    % extract mean of ROI from each beta
    for n = 1:length(betaXYZ),
        foo = spm_get_data(P(n),betaXYZ{n});
        std_roi(n) = std(foo);

        var_roi(n) = var(foo);
        mean_roi(n) = mean(foo(:));
    end
    
    if trimsd > 0,
        mean_roi = trimts(mean_roi, trimsd);
	 end
    
	 perc_covg(i)=length(find((mean_roi~=NaN))) / length(mean_roi);
	 if perc_covg>covg_limit
		beta_series(c,:)=mean_roi; c=c+1; idx=[i idx];
	 end
    clear foo
end
keyboard
% Scan through the ROIs
for i=1:size(beta_series,1)
    fprintf(['  ' n2sp(i,2) ': Processing ' ROI_list(i).name '\n']);
    % Look at the comparions
    for j=1:size(beta_series,1)
        % Now Cout(ROI,ROI)
        Cout(i,j)  = corr(beta_series(i,:)',beta_series(j,:)','type','Pearson');
    end
end

Network=Cout;
save([sav_dir 'Rcorr_' Events{1} '_Network.mat'],'Network'); 
clear Network;
% output R correlation results to image
Network=atanh(Cout);
save([sav_dir 'R_atanh_corr_' Events{1} '_Network.mat'],'Network');
clear Network;
Network=atanh(Cout)*sqrt(length(P)-3); 
save([sav_dir 'Zcorr_' Events{1} '_Network.mat'],'Network');
clear Network;

%=========================================================================%

function [y,ntrimmed] = trimts(y,sd)
% function [y,ntrimmed] = trimts(y,sd)
% y =  one D vector of data 
% sd = number of std - values out of this range are to be replaced
% ntrimmed = number of values replaced
ntrimmed = 0;
idx = find(abs(y) > sd*std(y));
if ~isempty(idx),
        y(idx) = sign(y(idx)) * sd*std(y);
        ntrimmed = length(idx);
end

function P = location_of_beta_images_from_event_discription(Events,SPM)
% P = location_of_beta_images_from_event_discription(Events,SPM)
% Events = Cell Array of Strings that defines a unique set of beta images
% SPM data structure

discription = {SPM.Vbeta.descrip}; % extract discription

% this block finds sets of index where discriptions match Events
Event_name = '';
for n = 1:length(Events),
    idx{n} = strfind(discription,Events{n});
    idx{n} = find(~cellfun('isempty',idx{n})); %strip out non matching results
    Event_name = [Event_name,'_',Events{n}];
end

% this block find the intersection of all sets of index
ref_idx = idx{1};
for n = 2:length(idx),
    reduced_idx{n-1} = intersect(idx{n},ref_idx);
end
final_idx = [];
for n = 1:length(reduced_idx),
    final_idx = union(reduced_idx{n},final_idx);
end

beta_images = {SPM.Vbeta(final_idx).fname}; % name of beta_images for selected events

% make array of location of beta images
for n = 1:length(beta_images)
    P{n} = fullfile(SPM.swd,beta_images{n});
end







