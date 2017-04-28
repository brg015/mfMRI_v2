% BRG 2014 Spring
%
function SL=RSA_model_design(SL)
% Load in stimuli information
data=excel_reader(SL.model.stimuli.file,...
    SL.model.stimuli.col_indx,...
    SL.model.stimuli.col_name);
% By default
stim_files=data{strcmp(SL.model.stimuli.col_name,'filename')}.col;
stim_ID_cell=data{strcmp(SL.model.stimuli.col_name,'ID')}.col;
stim_ID_num=zeros(1,length(stim_ID_cell));
for ii=1:length(stim_ID_cell), stim_ID_num(ii)=str2double(stim_ID_cell{ii}); end
stim_ID_num=1:length(stim_ID_cell);

for m=1:length(SL.model.design)
    if isfield(SL.design,'model')
        N=length(SL.design.model)+1;
    else
        N=1;
    end
    
    switch SL.model.design(m).type
%=========================================================================%
%% Create simple models (category models)
%=========================================================================%
% Simple models apply the value written to the specified column and create
% an RSA model based upon this. - also supports categorical anovas
case 'simple'
% 1) First ensure that we actually have this value
val_cell=data{strcmp(SL.model.stimuli.col_name,SL.model.design(m).name)}.col;
val_num=zeros(1,length(val_cell));
for ii=1:length(val_cell), val_num(ii)=str2double(val_cell{ii}); end
% 2) Set up the DRM
R=nan(length(val_cell),length(val_cell));
A=zeros(length(val_cell),length(val_cell));

% Setting up row by row
for ii=1:length(val_cell)
    % Assign matched values for item ii to mapped value col1 -> col2
    R(ii,(val_num(ii)==val_num))=...
        SL.model.design(m).map(SL.model.design(m).map(:,1)==val_num(ii),2);
    % Assign unmatched values based upon NaN criteria
    R(ii,(val_num(ii)~=val_num))=...
        SL.model.design(m).map(isnan(SL.model.design(m).map(:,1)),2);
    A(ii,(val_num(ii)==val_num))=...
        SL.model.design(m).map(SL.model.design(m).map(:,1)==val_num(ii),1);
end
U_vals=unique(val_num);
for ii=1:length(U_vals),
    anova_map{ii}=zeros(length(val_cell),length(val_cell));
    anova_map{ii}(A==U_vals(ii))=1;
end

% 3) Display drm to check
RSA_pdm(R,ones(1,length(val_num)),SL.model.design(m).name);

% 4) Save the output
if ~exist(SL.model.design(m).save_dir,'dir'),
    mkdir_tree(SL.model.design(m).save_dir);
end

save(fullfile(SL.model.design(m).save_dir,'model.mat'),'R','stim_ID_num','anova_map');

%=========================================================================%
%% Create complex models
%=========================================================================%           
case 'complex'

%=========================================================================%
%% Create custom models
%=========================================================================%
case 'custom' % Requires modification of RSA_custom
    SL=RSA_custom(SL,m);
%=========================================================================%
        otherwise
            error('Would you kindly specify a valid model');
    end
end

end

function SL=RSA_custom(SL,m)

switch SL.model.study
    case 'ER'
        switch SL.model.design(m).name
            case 'unclear'
                % Do stuff
            otherwise
                error('Would you kindly select a valid model');
        end
    otherwise
        error('Would you kindly select an existing study');
end

end


