%=========================================================================%
% Description
%=========================================================================%
function GT_shell_v3()
global SL ntwk;
sdisp('WARNING: Still in testing phases...',1);
%=========================================================================%
%% Presets
%=========================================================================%
% By defualt, include all subjects
if ~isfield(SL.dir,'include'), SL.dir.include=ones(1,length(SL.dir.subjects)); end

% This is specific to AAL atlas only
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
mask_dir=fullfile(SL.root_dir,'Geib');
switch ntwk.atlas
    case 'AAL_negX'
        mask_set='WFU_ROIs_resliced_f_negX';
        X=dir(fullfile(mask_dir,mask_set,'*.nii'));
        c=1; for ii=1:length(X),
            if (isempty(findstr('Verm',X(ii).name)) && isempty(findstr('Cerebel',X(ii).name)))
                SL.region.mask{c}=fullfile(mask_dir,mask_set,[X(ii).name]); c=c+1;
            end
        end  
    case 'AAL_posX'
        error('ROI.set does not exist yet');
    case 'HOA_negX'
        mask_set='Geib_HOA_resliced_negX';
        rm_ROI=[128,376,471,11,72,99,167,238,274,322,324,407]; %N=12
        c=1; for ii=setdiff(1:471,rm_ROI), % 459 ROIs ~30->80 voxels m
            SL.region.mask{c}=fullfile(mask_dir,mask_set,['HOA_' n2sp(ii,3) '.nii']); c=c+1;
        end
    case 'HOA_posX'
        mask_set='Geib_HOA_resliced_poxX';
        rm_ROI=[128,376,471,11,72,99,167,238,274,322,324,407]; %N=12
        c=1; for ii=setdiff(1:471,rm_ROI), % 459 ROIs ~30->80 voxels m
            SL.region.mask{c}=fullfile(mask_dir,mask_set,['HOA_' n2sp(ii,3) '.nii']); c=c+1;
        end
    case 'Crad200'
        mask_set='\Craddock\tcorr05_2level_all_ERmatch\';
        for ii=1:200
            SL.region.mask{ii}=fullfile(mask_dir,mask_set,['mask' num2str(ii) '.nii']);
        end      
end
%=========================================================================%
%% Analysis
%=========================================================================%
if SL.write.compute==1
    
for ii=1:length(ntwk.analysis)
    sdisp(['Computing: ' ntwk.analysis{ii}],1);
    % GT_write_rand() Defunct atm
    for jj=1:length(ntwk.contrast.nam)
        % Define the save directory
        save_dir_neg=fullfile(SL.dir.save,ntwk.contrast.neg{jj});
        save_dir_pos=fullfile(SL.dir.save,ntwk.contrast.pos{jj});
        save_dir_nam=fullfile(SL.dir.save,ntwk.contrast.nam{jj});
        if ~exist(save_dir_neg,'dir'), mkdir(save_dir_neg); end
        if ~exist(save_dir_pos,'dir'), mkdir(save_dir_pos); end
        if ~exist(save_dir_nam,'dir'), mkdir(save_dir_nam); end
 
        % Compute measures on both sets...
        GT_compute_measures_v2(ntwk.contrast.pos_R{jj},ii,save_dir_pos);
        
        if ~isempty(ntwk.contrast.neg_R)
            GT_compute_measures_v2(ntwk.contrast.neg_R{jj},ii,save_dir_neg);
            % Follow up with a comparison analysis as well - this will be a
            % different meta script to pull in measures from the pos and
            % neg analysis
            
            if ntwk.FSR==1
                for kk=1:length(SL.dir.subjects)
                    R1=load(ntwk.contrast.pos_R{jj}{kk},'R');
                    R2=load(ntwk.contrast.neg_R{jj}{kk},'R');
                    MAT(kk,:)=GT_FSR_v2(R1.R,R2.R);
                end
                save(fullfile(save_dir_nam,'FSR.mat'));
                save_asA=fullfile(save_dir_pos,'FSR_Tval.img');
                [~,p,~,stats]=ttest(MAT);
                [~]=ROI_constructor(SL.region.mask,stats.tstat,'','pure',save_asA,3);  
                data{1}.header='Tval';   data{1}.col=stats.tstat;
                data{2}.header='Zval';   data{2}.col=mean(MAT);
                data{3}.header='SEM';    data{3}.col=std(MAT)./sqrt(size(MAT,1));
                data{4}.header='pvalue'; data{4}.col=p;
                write_struct(data,fullfile(save_dir_nam,['FSR_tstats_full.csv']));
                clear data;
            end
            
            GT_compare_measures(ii,jj);
        end
        
    end
end
       
end
%=====================================================================%
%% Compute FSR Measures
%=====================================================================%
if SL.write.gephi==1, 
    GT_subnetwork; 
end








