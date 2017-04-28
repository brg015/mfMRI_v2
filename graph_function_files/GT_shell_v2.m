function GT_shell_v2()
global SL ntwk;
sdisp('WARNING: Still in testing phases...',1);
%=========================================================================%
%% Analysis Setup
%=========================================================================%
switch class(ntwk.analysis)
    case 'cell', loop_net=ntwk.analysis;
    case 'char', loop_net='';
end

ntwk.analysis='';
%=========================================================================%
%% Presets
%=========================================================================%
% By defualt, include all subjects
if ~isfield(SL.dir,'include'), SL.dir.include=ones(1,length(SL.dir.subjects)); end
if ~isfield(ntwk,'rm_hc'),     
    ntwk.rm_hc=0; 
    ntwk.str='';
elseif ntwk.rm_hc==1,
    ntwk.str='_no_hc';
elseif ntwk.rm_hc==2
    ntwk.str='_no_Rhc';
else
    ntwk.str='';
end

if ~exist(fullfile(SL.dir.save,'temp'),'dir'), 
    mkdir(fullfile(SL.dir.save,'temp'));
end
%=========================================================================%
%% Network Atlas
%=========================================================================%
mask_dir='D:\Data\Geib\';
switch ntwk.atlas
    case 'AAL_negX'
        mask_set='WFU_ROIs_resliced_f_negX';
        X=dir(fullfile(mask_dir,mask_set,'*.nii'));
        c=1; for ii=1:length(X),
            if (isempty(findstr('Verm',X(ii).name)) && isempty(findstr('Cerebel',X(ii).name)))
                SL.region.mask{c}=fullfile(mask_dir,mask_set,[SL.ROI.prefix X(ii).name]); c=c+1;
            end
        end  
    case 'AAL_posX'
        error('ROI.set does not exist yet');
    case 'HOA_negX'
        mask_set='Geib_HOA_resliced_negX';
        rm_ROI=[128,376,471,11,72,99,167,238,274,322,324,407]; %N=12
        c=1; for ii=setdiff(1:471,rm_ROI), % 459 ROIs ~30->80 voxels m
            SL.region.mask{c}=fullfile(mask_dir,mask_set,[SL.ROI.prefix 'HOA_' n2sp(ii,3) '.nii']); c=c+1;
        end
    case 'HOA_posX'
        mask_set='Geib_HOA_resliced_poxX';
        rm_ROI=[128,376,471,11,72,99,167,238,274,322,324,407]; %N=12
        c=1; for ii=setdiff(1:471,rm_ROI), % 459 ROIs ~30->80 voxels m
            SL.region.mask{c}=fullfile(mask_dir,mask_set,[SL.ROI.prefix 'HOA_' n2sp(ii,3) '.nii']); c=c+1;
        end
end

%=========================================================================%
%% Analysis
%=========================================================================%
% ntwk.rand (int) -> desribes which random networks run
%   1 -> G_effc
%   2 -> C_coef
%   3 -> Both
for ii=1:length(loop_net)
    ntwk.analysis=loop_net{ii};    
    sdisp(['Network: ' ntwk.analysis],1);
    for jj=1:ntwk.Niterations % is 1 unless random network
        sdisp(['Iteration ' num2str(jj)],1); ntwk.iteration=jj;
        %=================================================================%
        %% Analysis (Iteration & Ntwk)
        %=================================================================%
        switch ntwk.analysis
            case 'reg', 
                ntwk.measures={'C_coef','B_cent','E_cent','P_rank','D_cent','G_effc','C_lvge'};
                ntwk.save_name='network_measures';
            case 'LE'
                ntwk.measures={'L_effc'};
                ntwk.save_name='Subject_Networks_LE';
            case 'rand',    
                ntwk.measures={'C_coef_rand' 'G_effc_rand'};
                ntwk.save_name='Comparison_Networks_rand';
                ntwk.rand=3;
            case 'rand_effc'
                ntwk.measures={'G_effc_rand'};
                ntwk.save_name='Comparison_Networks_rand_G_effc';
                ntwk.rand=1;
            case 'rand_lat'
                ntwk.measures={'C_coef_rand_lat' 'G_effc_rand_lat'};
                ntwk.save_name='Comparison_Networks_lat';
                ntwk.rand=3;
            case 'mod'
                ntwk.measures={'Mod_L' 'Mod_LQ'};
                ntwk.save_name='Mod_Networks';
            case 'path_length'
                ntwk.measures={'PL'};
                ntwk.save_name='Path_Length';
            otherwise
                error('invalid network');
        end
        
        % Setup temporary save files
        if ntwk.Niterations==1
            save_data=fullfile(SL.dir.save,'RegNetworks',[ntwk.save_name ntwk.str '.mat']);
        else
            save_data=fullfile(SL.dir.save,'RandNetworks',...
                [ntwk.save_name '_N' num2str(ntwk.iteration) ntwk.str '.mat']);
        end
        
        if ((~exist(save_data,'file') || SL.dir.overwrite==1) && SL.write.dat==1),
            GT_compute_measures(save_data);
            A=load(save_data);
        else
           if exist(save_data,'file'), A=load(save_data); else error(['DNE: ' save_data]); end
        end
        %=================================================================%
        %% Work on non-rand networks
        %=================================================================%
        %--------------------%
        % Save output
        %--------------------%
        if SL.write.csv==1, GT_write_csv(A); end
        %--------------------%
        % Save output figures
        %--------------------%
        if SL.write.map==1, GT_write_maps(A); end
        %=================================================================%
    end % Iteration loop
    %=====================================================================%
    %% Write random networks
    %=====================================================================%
    if SL.write.rand==1, GT_write_rand; end
    ntwk.measures={}; ntwk.save_name={};
end

%=====================================================================%
%% Compute FSR Measures
%=====================================================================%
if SL.write.FSR==1, GT_FSR; end
if SL.write.gephi==1, GT_subnetwork; end








