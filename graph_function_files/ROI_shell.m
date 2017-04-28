function ROI_shell()
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
mean_beta_save=fullfile(SL.dir.save,'mean_beta_data.mat');
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
save(mean_beta_save,'meanBeta','SL.dir.subjects','SL.analysis.ROI.save','-v6');
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

% This needs moved/updated at some point

% if SL.analysis.ROI.on==1
% %-------------------------------------------------------------------------%
% %=========================================================================%
% %% ROI Analysis
% % Notes and stuff
% %=========================================================================%
% %-------------------------------------------------------------------------%
% sdisp('ROI Analysis',1);
% % What do we have at this point in time?
% % ANOVA{ROI}[subject X header]
% % AtoB{ROI}[subject X header]
% % Ident{ROI}[subject X header]
% %   Ident_series{ROI,header}[subject X trials]
% % Univariate{ROI}[subject X header]
% % SPEAR{ROI}[subject X header]
% %
% % _h -> has header information
% measures={'SPEAR' 'SPEAR_h' 'Anova_m' 'Anova_std' 'Anova_N' 'Anova_T' ...
%     'Anova_h' 'AtoB' 'AtoB_h' 'Ident_T' 'Ident_MMm' 'Ident_MMstd' ...
%     'Ident_MMN' 'Ident_Mm' 'Ident_Mstd' 'Ident_MN' 'Ident_h' 'Ident_series' 'Ident_delta_m' ...
%     'Univariate_T' 'Univariate_m' 'Univariate_std' 'Univariate_N' 'Univariate_h'};
% for ii=1:length(measures), eval([measures{ii} '={};']); end
% save_file=fullfile(SL.dir.save,'ROI','Compiled_Data.mat');
% if ~exist(fullfile(SL.dir.save,'ROI'),'dir'), mkdir(fullfile(SL.dir.save,'ROI')); end
% if (~exist(save_file,'file') || SL.analysis.ROI.overwrite==1),
% %=========================================================================%
% %% Compile data
% %=========================================================================%    
% [~,mask_name,~]=fileparts(SL.region.mask{1});
% X=load(fullfile(SL.dir.outpath,SL.dir.subjects{1},'data',[SL.ROI.prefix mask_name '.mat']));
% SL.design=X.design;    
% clear X;
% % The design only needs to be loaded such that the calculations are known,
% % it is incredibly critical that this is accurate!
% 
% for m=1:length(SL.region.mask)
%     [~,mask_name,~]=fileparts(SL.region.mask{m});
%     sdisp(SL.region.mask{m},2);
%     
%     % Initialize Counters:
%     SPEAR_c=1;
%     ANOVA_c=1;
%     AtoB_c=1;
%     Identity1_c=1;
%     Univariate_c=1;
% 
%     for ii=1:length(SL.design.calc)
%         for s=1:length(SL.dir.subjects) 
%             % Try to load the specified file, if it don't exist continue
%             Load_File=fullfile(SL.dir.outpath,SL.dir.subjects{s},...
%                 SL.design.calc{ii},[mask_name '_' SL.design.save_str{ii} '.mat']);
%             try
%                 load(Load_File);
%             catch err
%                 continue;
%             end
%             % Compile data
%             switch SL.design.calc{ii}
%                 case 'Spear'                      
%                     SPEAR{m}(s,SPEAR_c)=out_array;
%                     SPEAR_h{SPEAR_c}=[SL.design.save_str{ii}];
%                 case 'Anova1'
%                     Anova_m{m}(s,ANOVA_c)=out_array(2);
%                     Anova_std{m}(s,ANOVA_c)=out_array(3);
%                     Anova_N{m}(s,ANOVA_c)=out_array(4);
%                     Anova_T{m}(s,ANOVA_c)=out_array(5);
%                     Anova_h{ANOVA_c}=[SL.design.save_str{ii}];
%                 case 'AtoB'
%                     AtoB{m}(s,AtoB_c)=out_array;
%                     AtoB_h{AtoB_c}=[SL.design.save_str{ii}];
%                 case 'Identity1'
%                     Ident_T{m}(s,Identity1_c)=out_array(1);
%                     Ident_MMm{m}(s,Identity1_c)=out_array(2);
%                     Ident_MMstd{m}(s,Identity1_c)=out_array(3);
%                     Ident_MMN{m}(s,Identity1_c)=out_array(4);
%                     Ident_Mm{m}(s,Identity1_c)=out_array(5);
%                     Ident_Mstd{m}(s,Identity1_c)=out_array(6);
%                     Ident_MN{m}(s,Identity1_c)=out_array(7);
%                     Ident_delta_m{m}(s,Identity1_c)=out_array(5)-out_array(2);
%                     Ident_h{Identity1_c}=[SL.design.save_str{ii}];
%                     % Let's analyze a tad here
%                     % out_array_2(2,:) => Indices
%                     % out_array_2(1,:) => Values
%                     % Removing for now, not needed, and is difficult to
%                     % have a general solution to this at the moment
%                     % Ident_series{m,Identity1_c}(s,:)=nan(1,size(out_array_2(2,:),2));
%                     % Ident_series{m,Identity1_c}(s,out_array_2(2,:))=out_array_2(1,:);
%                 case 'Univariate'
%                     Univariate_T{m}(s,Univariate_c)=out_array(1);
%                     Univariate_m{m}(s,Univariate_c)=out_array(2);
%                     Univariate_std{m}(s,Univariate_c)=out_array(3);
%                     Univariate_N{m}(s,Univariate_c)=out_array(4);
%                     Univariate_h{Univariate_c}=[SL.design.save_str{ii}];
%                 case 'BetaSeries'
%                     % Handled by networks, no need to do aything   
%                 otherwise
%                     keyboard
%             end % Design Switch
%         end % Subject Loop
%         switch SL.design.calc{ii}
%             case 'Spear',      SPEAR_c=SPEAR_c+1;
%             case 'Anova1',     ANOVA_c=ANOVA_c+1;
%             case 'AtoB',       AtoB_c=AtoB_c+1;
%             case 'Identity1',  Identity1_c=Identity1_c+1;
%             case 'Univariate', Univariate_c=Univariate_c+1;
%         end
%     end % Calc Loop
% end
% 
% save(save_file,measures{:});
% SPEAR_h={};
% 
% % Create some contrast
% for ii=1:length(SL.analysis.ROI.contrast.pos)
%     % Need contrast across any analysis type
%     p=SL.analysis.ROI.contrast.pos{ii};
%     n=SL.analysis.ROI.contrast.neg{ii};
%     
%     % 1) Univariate to start
%     p_idx=strsearch(p,Univariate_h);
%     n_idx=strsearch(n,Univariate_h);
%     if (length(p_idx)==1 && length(n_idx)==1)
%         root_idx=findstr(p,Univariate_h{p_idx});
%         root_str=Univariate_h{p_idx}(1:root_idx-1);
%         idx=length(Univariate_h)+1;
%         Univariate_h{idx}=[root_str SL.analysis.ROI.contrast.nam{ii}];
%         % Calculate the difference
%         for jj=1:length(SL.region.mask)
%             % One sample t-test
%             Univariate_m{jj}(:,idx)=Univariate_m{jj}(:,p_idx)-Univariate_m{jj}(:,n_idx);
%             % Paired t-test? maybe in future
%         end 
%     end
%     clear p_idx n_idx
%     % 2) Anova
%     p_idx=strsearch(p,Anova_h);
%     n_idx=strsearch(n,Anova_h);
%     if (length(p_idx)==1 && length(n_idx)==1)
%         root_idx=findstr(p,Anova_h{p_idx});
%         root_str=Anova_h{p_idx}(1:root_idx-1);
%         idx=length(Anova_h)+1;
%         Anova_h{idx}=[root_str SL.analysis.ROI.contrast.nam{ii}];
%         % Calculate the difference
%         for jj=1:length(SL.region.mask)
%             % One sample t-test
%             Anova_m{jj}(:,idx)=Anova_m{jj}(:,p_idx)-Anova_m{jj}(:,n_idx);
%             % Paired t-test? maybe in future
%         end 
%     end
%     clear p_idx n_idx
%     % 3) AtoB
%     p_idx=strsearch(p,AtoB_h);
%     n_idx=strsearch(n,AtoB_h);
%     if (length(p_idx)==1 && length(n_idx)==1)
%         root_idx=findstr(p,AtoB_h{p_idx});
%         root_str=AtoB_h{p_idx}(1:root_idx-1);
%         idx=length(AtoB_h)+1;
%         AtoB_h{idx}=[root_str SL.analysis.ROI.contrast.nam{ii}];
%         % Calculate the difference
%         for jj=1:length(SL.region.mask)
%             % One sample t-test
%             AtoB{jj}(:,idx)=AtoB{jj}(:,p_idx)-AtoB{jj}(:,n_idx);
%             % Paired t-test? maybe in future
%         end 
%     end
%     clear p_idx n_idx
%     % 4) Spear    
%     p_idx=strsearch(p,SPEAR_h);
%     n_idx=strsearch(n,SPEAR_h);
%     if (length(p_idx)==1 && length(n_idx)==1)
%         root_idx=findstr(p,SPEAR_h{p_idx});
%         root_str=SPEAR_h{p_idx}(1:root_idx-1);
%         idx=length(SPEAR_h)+1;
%         SPEAR_h{idx}=[root_str SL.analysis.ROI.contrast.nam{ii}];
%         % Calculate the difference
%         for jj=1:length(SL.region.mask)
%             % One sample t-test
%             SPEAR{jj}(:,idx)=SPEAR{jj}(:,p_idx)-SPEAR{jj}(:,n_idx);
%             % Paired t-test? maybe in future
%         end 
%     end
%     clear p_idx n_idx
%     % 5) Ident is a tad different
%     p_idx=strsearch(p,Ident_h);
%     n_idx=strsearch(n,Ident_h);
%     if (length(p_idx)==1 && length(n_idx)==1)
%         root_idx=findstr(p,Ident_h{p_idx});
%         root_str=Ident_h{p_idx}(1:root_idx-1);
%         idx=length(Ident_h)+1;
%         Ident_h{idx}=[root_str SL.analysis.ROI.contrast.nam{ii}];
%         % Calculate the difference
%         for jj=1:length(SL.region.mask)
%             % One sample t-test
%             Ident_MMm{jj}(:,idx)=Ident_MMm{jj}(:,p_idx)-Ident_MMm{jj}(:,n_idx);
%             Ident_Mm{jj}(:,idx)=Ident_Mm{jj}(:,p_idx)-Ident_Mm{jj}(:,n_idx);
%             Ident_delta_m{jj}(:,idx)=Ident_Mm{jj}(:,idx)-Ident_MMm{jj}(:,idx);
%             % Paired t-test? maybe in future
%         end 
%     end
%     clear p_idx n_idx
%     clear p n
% end
% 
% save(save_file,measures{:}); % Resave w/ contrast
% 
% else
%     load(save_file);
% end
% 
% %=========================================================================%
% % Analyze that data
% %=========================================================================%
% % Save a barrage of T-maps
% if ~exist(fullfile(SL.dir.save,'ROI','data'),'dir'),
%     mkdir(fullfile(SL.dir.save,'ROI','data'));
% end
% if ~exist(fullfile(SL.dir.save,'ROI','Tmaps'),'dir'),
%     mkdir(fullfile(SL.dir.save,'ROI','Tmaps'));
% end
% calc_maps={'AtoB','Univariate_m','SPEAR','Ident_MMm','Ident_Mm','Ident_delta_m'};
% calc_head={'AtoB_h','Univariate_h','SPEAR_h','Ident_h','Ident_h','Ident_h'};
% for ii=1:length(calc_maps)
%     data=eval(calc_maps{ii});
%     head=eval(calc_head{ii});
%     if ~isempty(head)
%         for jj=1:size(data{1},2)
%             % Write Outmaps
%             save_as=fullfile(SL.dir.save,'ROI','Tmaps',[head{jj} '.nii']);
%             ROI_constructor(SL.region.mask,data,jj,'T',save_as);
%             % Reshape and output values as well
%             save_as_also=fullfile(SL.dir.save,'ROI','data',[head{jj} '.mat']);
%             for kk=1:length(data)
%                out_data{1}(kk,:)=data{kk}(:,jj); 
%             end
%             save(save_as_also,'out_data');
%             clear out_data;
%         end
%     end
% end
% 
% %=========================================================================%
% return;
% % Lots of work needed down here to generalize:
% 
% % Loop through all those masks
% keyboard
% save_dir='D:\Data\Wing\ERMatch_serv1\Analysis\STwa\SL_ER_vols_5vox_sd2_10N_features\ROI\Group\';
% imp_dir=fullfile(save_dir,'Significant');
% 
% % Spearman Settings
% s_model_match={'features_ER_nan','features_unique','features_semantic','gbjet_nan'};
% s_contrast=[0 0 0 -1 1 0 0 0 0;...
%     0 0 0 0 0 -1 0 1 0;...
%     0 0 0 0 0 0 -1 0 1;...
%     0 0 0 0 0 -1 1 0 0;...
%     0 0 0 0 0 0 0 -1 1;...
%     1 0 0 0 0 0 0 0 0;...
%     0 1 0 0 0 0 0 0 0;...
%     0 0 1 0 0 0 0 0 0];
% 
% s_names={'ER_RF' 'R_RF' 'E_RF' 'EvsR_F' 'EvsR_R' 'E' 'R' 'ER'};
% L={'E' 'R' 'ER' 'ER12' 'ER34' 'R12' 'E12' 'R34' 'E34'};
% 
% % Ident Settings
% thresh=[-1 1];
% FModel=load('D:\Data\Wing\ERMatch_serv1\RSA_models\features_ER_nan\model.mat');
% FS_base=nanmean(nanmean(FModel.R));
% data=excel_reader('D:\Data\Wing\ERMatch_serv1\RSA_models\features_semantic_nan\scenes.csv');
% c=1;
% for ii=2:length(data),
%     SModel(c,:)=cell2num(data{ii}.col); c=c+1;
% end
% Sitem=data{1}.col; clear data;
% IdentLegend={'Base' 'High_ER_Sim' 'Low_ER_Sim'};
% 
% for m=1:length(SL.region.mask)
%     [~,name,~]=fileparts(SL.region.mask{m});
%     sdisp(name,1);
%     %=====================================================================%
% %     % Spearman Test
% %     %=====================================================================%
% %     % Let's examine the Spearman feature stuff first
% %     save_spear=fullfile(save_dir,'Spear');
% %     if ~exist(save_spear,'dir'), mkdir(save_spear); end
% %     
% %     %==============================%
% %     % Loop through tested models
% %     %==============================%
% %     for ii=1:length(s_model_match),
% %         sdisp(s_model_match{ii},2);
% %         % Find matching models
% %         I=strsearch(s_model_match{ii},SPEAR_h);
% %         D=SPEAR{m}(:,I);
% %         me=mean(D);
% %         d=std(D);
% %         n=size(D,1);
% %         T=me./(d/sqrt(n));
% %         clear me d I;
% %         
% %         % Compute Basic models
% % %         fig_save=fullfile(save_spear,[name '_Basic_' s_model_match{ii} '.png']);
% % %         figure; bar(T); ylabel('T-value'); ylim([-5 5]); grid;
% % %         set(gca,'XTickLabel',L);
% % %         set(gcf,'position',[0 0 1280 1024]);
% % %         export_fig(fig_save); close all;
% %    
% %         DT(ii,6:8)=T(1:3); % difference
% %         PT(ii,6:8)=T(1:3); % paired T
% %         clear T;
% %         
% %         % Now Compute Complex models
% %         for jj=1:5 % Ignore duplicates from basic
% %            c_vect=s_contrast(jj,:); 
% %            POS=D(:,c_vect==1);
% %            NEG=D(:,c_vect==-1);
% %            % Diff Calc.
% %            Delta=POS-NEG;
% %            Dm=mean(Delta);
% %            Dd=std(Delta);
% %            DT(ii,jj)=Dm/(Dd/sqrt(n));
% %            % Paired Calc.
% %            POSm=mean(POS);
% %            POSd=std(POS);
% %            NEGm=mean(NEG);
% %            NEGd=std(NEG);
% %            PT(jj)=(POSm-NEGm)/(sqrt((POSd/n)+(NEGd/n)));
% %            P(ii,jj)=POSm;
% %            N(ii,jj)=NEGm;
% %            clear c_vect POS NEG Delta Dm Dd;
% %            clear POSm POSd NEGm NEGd;
% %            % Returns DT, PT
% %         end
% %                 
% % %         fig_save=fullfile(save_spear,[name '_Complex_T_' s_model_match{ii} '.png']);
% % %         figure; bar(DT(ii,:)); ylabel('T-value'); ylim([-5 5]); grid;
% % %         set(gca,'XTickLabel',s_names);
% % %         set(gcf,'position',[0 0 1280 1024]);
% % %         export_fig(fig_save); close all;
% %         
% % %         fig_save=fullfile(save_spear,[name '_Complex_PT_' s_model_match{ii} '.png']);
% % %         figure; bar(PT); ylabel('T-value'); ylim([-5 5]); grid;
% % %         set(gca,'XTickLabel',s_names);
% % %         set(gcf,'position',[0 0 1280 1024]);
% % %         export_fig(fig_save); close all;
% %         
% %         clear D PT;
% %         
% %     end % End of model loop
% %     
% %     % Special Pictures:
% %     % P [Model X Contrast]
% %     % N [Model X Contrast]
% %     % P-N -> DT model
% %     % DT [Model X Contrast]
% %     fig_save=fullfile(save_dir,'Spear_Summary',[name '_T.png']);
% %     figure; bar(DT(1:3,1:3)'); ylabel('T-value','fontsize',16); ylim([-5 5]); grid;
% %     set(gca,'XTickLabel',s_names(1:3),'FontSize',16);
% %     legend({'Original' 'Uncommon' 'Common'});
% %     set(gcf,'position',[0 0 1280 1024]);
% %     export_fig(fig_save); close all;
% %    
% %     V=[P(1:3,1:3)' N(1:3,1:3)'];
% %     Z(:,1)=V(:,1);
% %     Z(:,2)=V(:,4);
% %     Z(:,3)=V(:,2);
% %     Z(:,4)=V(:,5);
% %     Z(:,5)=V(:,3);
% %     Z(:,6)=V(:,6);
% %     fig_save=fullfile(save_dir,'Spear_Summary',[name '_Rho.png']);
% %     figure; bar(Z); ylabel('Rho','fontsize',16); ylim([-.1 .1]); grid;
% %     set(gca,'XTickLabel',s_names(1:3),'FontSize',16);
% %     legend({'O-Hit' 'O-Miss' 'C-Hit' 'C-Miss' 'U-Hit' 'U-Miss'});
% %     set(gcf,'position',[0 0 1280 1024]);
% %     export_fig(fig_save); close all;
% 
%     %=====================================================================%
%     % Erik Craziness
%     %=====================================================================%  
%     sdisp('Ident',2);
%     % Let's examine the Spearman feature stuff first
%     save_ident=fullfile(save_dir,'Identity');
%     if ~exist(save_ident,'dir'), mkdir(save_ident); end
%             
%     % Ident_series{m,Identity1_c}(s,:) => {mask X contrast}(subject,:)
%     for jj=3 % For now look at all
%         D=Ident_series{m,jj};
%         ZscoreD=NaN(size(D));
%         n=size(D,1);
%         for ii=1:n
%             temp=D(ii,:);
%             nan_idx=~isnan(temp);
%             ZscoreD(ii,nan_idx)=zscore(D(ii,nan_idx));
%         end
%         % For now ignore NaN 'cause ER diag.
%         Zm=mean(ZscoreD);
%         Zd=std(ZscoreD);
%         ZT=Zm./(Zd/sqrt(n));
%         
%         HI=ZT>thresh(2);
%         LI=ZT<thresh(1);
%         
%         FS_low=nanmean(nanmean(FModel.R(LI,LI)));
%         FS_high=nanmean(nanmean(FModel.R(HI,HI)));
%         
%         FS_vect=[FS_base FS_high FS_low];
%         fig_save=fullfile(save_ident,[name '_' Ident_h{jj} '.png']);
%         figure; bar(FS_vect); ylabel('Spearman Rho'); ylim([0 0.1]); grid;
%         set(gca,'XTickLabel',IdentLegend);
%         set(gcf,'position',[0 0 1280 1024]);
%         export_fig(fig_save); close all;
%     end
%     
%     for jj=1:2 
%         D=Ident_series{m,jj};
%         n=size(D,1);
%         for ii=1:n
%             nan_idx=~isnan(D(ii,:));
%             FS(ii,jj)=nanmean(nanmean(FModel.R(nan_idx,nan_idx)));
%         end
%     end
%     
%     FS_high_m=mean(FS(:,1));
%     FS_high_d=std(FS(:,1));
%     FS_low_m=mean(FS(:,2));
%     FS_low_d=std(FS(:,2));
%     
%     FS_vect=[FS_base FS_high_m FS_low_m];
%     fig_save=fullfile(save_ident,[name '_RemFor_v.png']);
%     figure; bar(FS_vect); ylabel('Spearman Rho'); ylim([0 0.1]); grid;
%     set(gca,'XTickLabel',{'Base' 'Rem_ER_Sim' 'For_ER_Sim'});
%     set(gcf,'position',[0 0 1280 1024]);
%     export_fig(fig_save); close all;
% 
%     T=(FS_high_m-FS_low_m)/(sqrt(FS_high_d/n+FS_low_d/n));
%     
%     fig_save=fullfile(save_ident,[name '_RemFor_T.png']);
%     figure; bar(T); ylabel('T-value'); ylim([-5 5]); grid;
%     set(gcf,'position',[0 0 1280 1024]);
%     export_fig(fig_save); close all;
% end
% 
% 
% 
% % for ii=1:21
% %    I=find(MINC(:,ii));
% % end
% % 
% % 
% 
% 
% 
% end


















