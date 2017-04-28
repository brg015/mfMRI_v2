function ROI_shell_v2()
global SL;
%=========================================================================%
% Set defaults
%=========================================================================%
% Save dir is autoset
SL.dir.save=fullfile(SL.dir.outpath,'Group\',SL.ROI.set);
if ~isfield(SL,'corr'), SL.corr='Pearson'; end

%=========================================================================%
% Add in the masks
%=========================================================================%
% AAL -> 90 masks
% HOA -> 471 masks
mask_dir=fullfile(SL.dir.root,'Geib\');
switch SL.ROI.set
    case 'AAL_negX'
        mask_set='WFU_ROIs_resliced_f_negX';
        X=dir(fullfile(mask_dir,mask_set,'*.nii'));
        c=1; for ii=1:length(X),
            if (isempty(findstr('Verm',X(ii).name)) && isempty(findstr('Cerebel',X(ii).name)))
                SL.region.mask{c}=fullfile(mask_dir,mask_set,[X(ii).name]); c=c+1;
            end
        end
    case 'AAL_posX'
        mask_set='WFU_ROIs_resliced_f_posX';
        X=dir(fullfile(mask_dir,mask_set,'*.nii'));
        c=1; for ii=1:length(X),
            if (isempty(findstr('Verm',X(ii).name)) && isempty(findstr('Cerebel',X(ii).name)))
                SL.region.mask{c}=fullfile(mask_dir,mask_set,[X(ii).name]); c=c+1;
            end
        end
    case 'HOA_negX'
        mask_set='Geib_HOA_resliced_negX';
        c=1; for ii=1:471, % 459 ROIs ~30->80 voxels m
            SL.region.mask{c}=fullfile(mask_dir,mask_set,['HOA_' n2sp(ii,3) '.nii']); c=c+1;
        end
    case 'HOA_posX'
        mask_set='Geib_HOA_resliced_posX';
        c=1; for ii=1:471, % 459 ROIs ~30->80 voxels m
            SL.region.mask{c}=fullfile(mask_dir,mask_set,['HOA_' n2sp(ii,3) '.nii']); c=c+1;
        end
    case 'Memex100'     
        mask_dir='D:\Data\Memex\MVPA_analysis_bg\ROI\functional_HOA100\';
        X=dir([mask_dir 'ROI*.nii']);
        for ii=1:length(X)
            SL.region.mask{ii}=fullfile(mask_dir,X(ii).name);
        end
end
%=========================================================================%
%% Network Creation
% Returns:
%   RDM{subject,network}[data X ROI]     => RDM networks
%   Beta{subject,network}[data X ROI]    => Beta networks
%   MINC(subject,mask)                   => Included mask
%=========================================================================%
% RDM_save=fullfile(SL.dir.save,'RDM_data.mat');
beta_save=fullfile(SL.dir.save,'beta_data.mat');
% mean_beta_save=fullfile(SL.dir.save,'mean_beta_data.mat');
if (~exist(beta_save,'file') || SL.analysis.ROI.overwrite==1)

sdisp('Network Analysis',1);
for m=1:length(SL.region.mask)
    [~,mask_name,~]=fileparts(SL.region.mask{m});
    sdisp(SL.region.mask{m},2);
    for ii=1:length(SL.dir.subjects)
        % Loads in R, tmp, design
        LF=fullfile(SL.dir.outpath,SL.dir.subjects{ii},'data',[SL.ROI.prefix mask_name '.mat']);
        % Load in the file...
        if exist(LF,'file'), load(LF); MINC(m,ii)=1; else 
            MINC(m,ii)=0; 
            display([' ' SL.dir.subjects{ii} ': Missing mask']);
            display(['  DNE: ' LF]);
            continue; 
        end
        % Added variables
        % design (from SL.design)
        % R = corr(tmp)
        % tmp [voxels X SL.files]
        cBox=cumsum(design.Box);

        % Here we make the networks
        for jj=1:length(SL.analysis.ROI.save)
            % From Ibox, determine which cells to include, then clearly
            % identify these columns
            I=find(SL.analysis.ROI.Ibox(jj,:));
            cI=[]; % Define column index
            v1=[]; % Define RDM
            for kk=1:length(I)
                if I(kk)==1
                    cI=[cI, 1:cBox(1)];
                else
                    cI=[cI, cBox(I(kk)-1)+1:cBox(I(kk))];
                end
            end        
            % v1=corr(tmp(:,cI));  
            mean_Beta(m,ii,jj)=mean(mean(tmp(:,cI)));
            if m==1, 
%                 RDM{ii,jj}=v1; 
                Beta{ii,jj}=mean(tmp(:,cI),1)';
            else
%                 RDM{ii,jj}=[RDM{ii,jj},v1];
                Beta{ii,jj}=[Beta{ii,jj},mean(tmp(:,cI),1)'];
            end
        end
        clear R tmp design;
    end
end

% Save included masks to output
if ~exist(beta_save,'dir'), mkdir_tree(fullfile(SL.dir.save,'beta')); end
% if ~exist(RDM_save,'dir'), mkdir_tree(fullfile(SL.dir.save,'RDM')); end

save(fullfile(SL.dir.save,'beta','IncludedROIs.mat'),'MINC');
% save(fullfile(SL.dir.save,'RDM','IncludedROIs.mat'),'MINC');
data{1}.header='ID';
data{2}.header='Mask';
for ii=1:length(SL.region.mask),
    data{1}.col(ii)=ii;
    [~,mname,~]=fileparts(SL.region.mask{ii});
    data{2}.col{ii}=mname;
end
write_struct(data,fullfile(SL.dir.save,'ROI_list.csv')); clear data;
save(beta_save,'Beta','-v7.3');
% save(mean_beta_save,'meanBeta','SL.dir.subjects','SL.analysis.ROI.save','-v6');
% save(RDM_save,'RDM','-v7.3');

else
    sdisp('Loading Data...',1);
    load(beta_save); % load(RDM_save);
end
%=========================================================================%
% Micro connectivity script
%=========================================================================%
for jj=1:length(SL.dir.subjects),
    for ii=1:length(SL.analysis.ROI.save)    
        % Beta Network
        Bdir=fullfile(SL.dir.save,'beta',SL.analysis.ROI.save{ii});
        if ~exist(Bdir,'dir'), mkdir(Bdir); end
        switch SL.corr
            case 'Pearson', RBeta=corr(Beta{jj,ii});
            case 'dcorr',   
                X=Beta{jj,ii};
                for aa=1:size(X,2)
                    for bb=1:size(X,2)
                        RBeta(aa,bb)=dcorr(X(:,aa),X(:,bb));
                    end
                end
        end
        R=RBeta; m0=~eye(size(R));
        R1=R.*m0; R1(R1<0)=0; 
        switch SL.corr
            case 'Pearson', save(fullfile(Bdir,[SL.dir.subjects{jj} '.mat']),'R','R1','-v7.3');
            case 'dcorr', save(fullfile(Bdir,[SL.dir.subjects{jj} '_dcorr.mat']),'R','R1','-v7.3');
        end
        clear R R1;
        % RDM Network
%         Rdir=fullfile(SL.dir.save,'RDM',SL.analysis.ROI.save{ii});
%         if ~exist(Rdir,'dir'), mkdir(Rdir); end
%         RDRM=corr(RDM{jj,ii});
%         R=RDRM; m0=~eye(size(R));
%         R1=R.*m0; R1(R1<0)=0; 
%         save(fullfile(Rdir,[SL.dir.subjects{jj} '.mat']),'R','R1','-v7.3');
%         clear R R1;            
    end
end



















