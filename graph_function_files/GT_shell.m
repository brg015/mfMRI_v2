function GT_shell()

global SL;
global ntwk;

if ~isfield(SL.dir,'no_subj_folder'), SL.dir.no_subj_folder=0; end
if ~isfield(SL,'eeg'), SL.eeg=0; end
if ~isfield(SL.dir,'include'), SL.dir.include=ones(1,length(SL.dir.subjects)); end

if ntwk.Niterations==1
    save_data=fullfile(SL.dir.save,[ntwk.save_name '.mat']);
else
    save_data=fullfile(SL.dir.save,'RandNetworks',[ntwk.save_name '_N' num2str(ntwk.iteration) '.mat']);
end

if ~exist(fullfile(SL.dir.save,'temp'),'dir'), 
    mkdir(fullfile(SL.dir.save,'temp'));
end
%=========================================================================%
% Network Analysis
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
    case 'AAL_negX_noHc'
        mask_set='WFU_ROIs_resliced_f_negX';
        X=dir(fullfile(mask_dir,mask_set,'*.nii'));
        keyboard;
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

if ((~exist(save_data,'file') || SL.dir.overwrite==1) && SL.write.dat==1),
%=========================================================================%
% Analysis
%=========================================================================%
for ii=1:length(SL.analyses)
    display(['Running: ' SL.analyses{ii}]);
    display(datestr(now))
    for jj=1:length(SL.dir.subjects)
        % If temp files exist, load them, and then continue
        if ntwk.Niterations==1
            if (exist(fullfile(SL.dir.save,'temp',[ntwk.save_name '_' num2str(ii) '_' num2str(jj) '.mat']),'file') ...
                    && SL.dir.overwrite==0);
                load(fullfile(SL.dir.save,'temp',[ntwk.save_name '_' num2str(ii) '_' num2str(jj) '.mat'])); continue;
            end
        end
        
        display(['  Running: ' SL.dir.subjects{jj}]);
        switch SL.dir.no_subj_folder
            case 0, load_file=fullfile(SL.dir.data,SL.analyses{ii},[SL.dir.subjects{jj} '.mat']);
                    load(load_file);
            % Specific to EEG at the moment, treating it as such
            case 1, load_file=fullfile(SL.dir.data,[SL.dir.subjects{jj} '_' SL.analyses{ii} '.mat']);
                    load(load_file);
                    if SL.eeg==1
                        % R1=dat.R;
                        for kk=1:size(dat.I,1)
                            [r1,c1]=cell_block(kk,kk,dat.Box);
                            R1{kk}=dat.R(r1,c1);
                            clear r1 c1;
                        end
                    end
        end
        
        if iscell(R1), L=length(R1); TR=R1; clear R1; else L=1; end
        
        % If mulitple networks, for example in EEG
        for kk=1:L 
            if L>1, R1=TR{kk}; end
            switch ntwk.analysis
                case 'rand_effc'
                    % Rand1=randmio_und_connected(R1, 10);
                    Rand1 = null_model_und_sign(R1); % Changed 12/10
                    G_effc_rand{ii,jj}=efficiency_wei_node(Rand1).^(-1);  % changed 12/10 => char. path
                    if ntwk.Niterations==1
                        save(fullfile(SL.dir.save,SL.analyses{ii},[SL.dir.subjects{jj} '_rand.mat']),'Rand1');
                    end
                    clear Rand1
                case 'rand'
                    % Rand1=randmio_und_connected(R1, 10);
                    Rand1 = null_model_und_sign(R1); % Changed 12/10
                    C_coef_rand{ii,jj}=clustering_coef_wu(Rand1); 
                    G_effc_rand{ii,jj}=efficiency_wei_node(Rand1).^(-1);  % changed 12/10 => char. path
                    if ntwk.Niterations==1
                        save(fullfile(SL.dir.save,SL.analyses{ii},[SL.dir.subjects{jj} '_rand.mat']),'Rand1');
                    end
                    clear Rand1
                case 'rand_lat'
                    Rand2=latmio_und_connected(R1,10);
                    C_coef_rand_lat{ii,jj}=clustering_coef_wu(Rand2);
                    G_effc_rand_lat{ii,jj}=efficiency_wei(Rand2); 
                    save(fullfile(SL.dir.save,SL.analyses{ii},[SL.dir.subjects{jj} '_latrand.mat']),'Rand2');
                    clear Rand2
                case 'reg'
                    C_coef{ii,jj}=clustering_coef_wu(R1); 
                    B_cent{ii,jj}=betweenness_wei(R1); 
                    E_cent{ii,jj}=eigenvector_centrality_und(R1); 
                    P_rank{ii,jj}=pagerank_centrality(R1, 0.85); 
                    % Degree Centrality (sum row/col divided by 2)
                    D_cent{ii,jj}=sum(R1)./2; 
                    G_effc{ii,jj}=efficiency_wei_node(R1); 
                    [l ~] = wfu_compute_leverage([],R1,[]);
                    C_lvge{ii,jj}=l;
                case 'reg_GE'
                    % C_coef{ii,jj}=clustering_coef_wu(R1); 
                    G_effc{ii,jj}=efficiency_wei_node(R1); 
                    C_coef{ii,jj}=clustering_coef_wu(R1);
                case 'LE'
                    L_effc{ii,jj}=efficiency_wei_node(R1,1); % Takes 1.5 hrs per loop
                case 'mod'
                    % true: https://sites.google.com/site/bctnet/construction
                    R(logical(eye(size(R))))=0; % Fix diag on R
                    [Ci Q]=modularity_louvain_dir(R);
    %                 [Ci Q]=modularity_louvain_und(R1);
                    % Find top nodes
    %                 Rbin=threshold_proportional(R1,0.10); % Threshold
    %                 Rbin(Rbin>0)=1;                       % Binarize
    %                 [Ci Q]=modularity_louvain_und(Rbin);

                    Mod_L{ii,jj}=Ci;
                    Mod_LQ{ii,jj}=Q;
                    % [CIJscore,sn] = score_wu(R1,s)
                    % [Rw] = rich_club_wd(R,20);
                    % R_club{ii,jj}=Rw;
                    clear Ci Q;
            end
            
            if (L>1 && SL.eeg==1)
                % For each ntwk measure
                for aa=1:length(ntwk.measures)
                   ndata=eval([ntwk.measures{aa} '{ii,jj}']);
                   eval([ntwk.measures{aa} '_TimeFreqVal{ii,jj}(dat.I(kk,1),dat.I(kk,2))=mean(ndata);']);
                   clear ndata;
                end
            end
            clear R R1; % Clear original network
        end
        
        % Save progress as we go, cause it seems to like to crash
        if ntwk.Niterations==1
            save(fullfile(SL.dir.save,'temp',[ntwk.save_name '_' num2str(ii) '_' num2str(jj)]),ntwk.measures{:});
        end
    end % Subject Loop
    display(datestr(now))
end % Analy Loop

%=====================%
% Special EEG Save
%=====================%
if (L>1 && SL.eeg==1)
    for aa=1:length(ntwk.measures)
       for bb=1:length(SL.analyses)
           for cc=1:length(SL.dir.subjects)
               D{bb}(cc,:,:)= eval([ntwk.measures{aa} '_TimeFreqVal{bb,cc}']);
           end
           % D{analysis}[Subject X Freq X Time]
       end
       dDm=squeeze(mean(D{1}-D{2}));
       dDs=squeeze(std(D{1}-D{2}));
       dDT=dDm./(dDs/sqrt(size(D{1},1)));
       
       d1m=squeeze(mean(D{1}));
       d1s=squeeze(std(D{1}));
       d1T=d1m./(d1s/sqrt(size(D{1},1)));
       
       d2m=squeeze(mean(D{2}));
       d2s=squeeze(std(D{2}));
       d2T=d2m./(d2s/sqrt(size(D{1},1)));
       
       ddDT=(d1m-d2m)./sqrt(d1s.^2/size(D{1},1)+d2s.^2/size(D{1},1));
    end
end

save(save_data,ntwk.measures{:});
%=========================================================================%
else
   if exist(save_data,'file'),
       load(save_data);
   else 
%        % Resurect these files
%        for kk=1:length(ntwk.measures)
%            for ii=1:length(SL.analyses),
%                save_file=fullfile(SL.dir.save,[SL.analyses{ii} '_' ntwk.measures{kk} '.csv']);
%                data=excel_reader(save_file);
%                for jj=1:length(SL.dir.subjects)
%                    eval([ntwk.measures{kk} '{ii,jj}=cell2num(data{jj}.col);']);
%                end
%            end
%        end
%        save(save_data,ntwk.measures{:});
   end
end
%=========================================================================%
if SL.write.rand==1
%=========================================================================%
for ii=1:ntwk.Niterations
    LV=fullfile(SL.dir.save,'RandNetworks',[ntwk.save_name '_N' num2str(ii) '.mat']);
    RND=load(LV);
    for jj=1:size(RND.C_coef_rand,1)
        for kk=1:size(RND.C_coef_rand,2)
            Clust(ii,jj,kk)=mean(RND.C_coef_rand{jj,kk}); % Average across subject nodes
            PL(ii,jj,kk)=mean(RND.G_effc_rand{jj,kk});    % Average across subject nodes
        end
    end
    clear RND;
end
mClust=squeeze(mean(Clust,1));
mPL=squeeze(mean(PL,1));
for ii=1:size(mClust,1)
    for jj=1:size(mClust,2)
        Cdata{jj}.header=SL.dir.subjects{jj};
        Cdata{jj}.col=mClust(ii,jj);
        Pdata{jj}.header=SL.dir.subjects{jj};
        Pdata{jj}.col=mPL(ii,jj);
    end
    write_struct(Cdata,fullfile(SL.dir.save,SL.analyses{ii},'RandCluster.csv'));
    write_struct(Pdata,fullfile(SL.dir.save,SL.analyses{ii},'RandPathLength.csv'));
    clear Cdata Pdata;
end
        
mmClust=squeeze(mean(mClust,2));
mmPL=squeeze(mean(mPL,2));

keyboard
%=========================================================================%
end
%=========================================================================%
if SL.write.csv==1
%=========================================================================%
% Save output
%=========================================================================%
if exist(fullfile(SL.dir.save,'IncludedROIs'),'file')
    X=load(fullfile(SL.dir.save,'IncludedROIs'));
    X.MINC=logical(X.MINC);
    noROI=0;
else
    noROI=1;
end

% X.MINC=ones(size(X.MINC));
% 1) Check for wanted contrast
N=length(SL.analyses);
for ii=1:length(ntwk.contrast.nam)
    SL.analyses{end+1}=ntwk.contrast.nam{ii};
end

for ii=1:length(SL.analyses)
   for kk=1:length(ntwk.measures)
       save_file=fullfile(SL.dir.save,[SL.analyses{ii} '_' ntwk.measures{kk} '.csv']);
       if (~exist(save_file,'file') || SL.dir.overwrite==1)
           for jj=1:length(SL.dir.subjects)
               if noROI==1
                   if ii<=N % Non contrast analysis
                       data{jj}.header=SL.dir.subjects{jj};
                       data{jj}.col=eval([ntwk.measures{kk} '{ii,jj}']);
                   else
                       PC=strcmp(ntwk.contrast.pos{ii-N},SL.analyses);
                       NC=strcmp(ntwk.contrast.neg{ii-N},SL.analyses);
                       if ~isempty(PC) && ~isempty(NC),
                           data{jj}.header=SL.dir.subjects{jj};
                           data{jj}.col=eval([ntwk.measures{kk} '{PC,jj}'])-eval([ntwk.measures{kk} '{NC,jj}']);
                       end
                   end
               else
                   if ii<=N % Non contrast analysis
                       data{jj}.header=SL.dir.subjects{jj};
                       data{jj}.col(X.MINC(:,jj))=eval([ntwk.measures{kk} '{ii,jj}']);
                       data{jj}.col(~X.MINC(:,jj))=NaN; % Unincluded=NaN
                   else
                       PC=strcmp(ntwk.contrast.pos{ii-N},SL.analyses);
                       NC=strcmp(ntwk.contrast.neg{ii-N},SL.analyses);
                       if ~isempty(PC) && ~isempty(NC),
                           data{jj}.header=SL.dir.subjects{jj};
                           data{jj}.col(X.MINC(:,jj))=...
                               eval([ntwk.measures{kk} '{PC,jj}'])-eval([ntwk.measures{kk} '{NC,jj}']);
                           data{jj}.col(~X.MINC(:,jj))=NaN; % Unincluded=NaN
                       end
                   end
               end
           end % Subject loop for each node and each measure
           write_struct(data,save_file);
           clear data;
       end
   end % Measures Loop
end % Analyis Loop

% Check if Univariate output exists as well, if so, save it
if exist(SL.dir.univ,'dir'),
    F=dir([SL.dir.univ filesep '*.mat']);
    for ii=1:length(F)
        D=load(fullfile(SL.dir.univ,F(ii).name));
        for jj=1:size(D.out_data{1},2),
            data{jj}.header=SL.dir.subjects{jj};
            data{jj}.col=D.out_data{1}(:,jj);
        end
        [~,cor,~]=fileparts(F(ii).name);
        write_struct(data,fullfile(SL.dir.save,['Univ_' cor '.csv']));
        clear data;
    end
end

%=========================================================================%
end
if SL.write.map==1
%=========================================================================%
% Save output
%=========================================================================%
N=length(SL.analyses);
for ii=1:length(ntwk.contrast.nam)
    SL.analyses{end+1}=ntwk.contrast.nam{ii};
end
X=load(fullfile(SL.dir.save,'IncludedROIs'));
X.MINC=logical(X.MINC);

for ii=1:length(ntwk.measures)
    data=eval(ntwk.measures{ii}); 
    sdisp(ntwk.measures{ii},2);
    
    for jj=1:N
        c=1;
        for kk=1:length(SL.dir.subjects)  
            if SL.dir.include(kk)==1
                temp=data{jj,c}; data{jj,c}=[];
                data{jj,c}(X.MINC(:,c))=temp;
                data{jj,c}(~X.MINC(:,c))=NaN;
                sdata{jj,c}=data{jj,c}';
                c=c+1;
            end
        end
    end
    clear data; data=sdata;
    
    for jj=1:length(SL.analyses)
    
        save_asA=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_avg.nii']);
        save_asB=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_T.nii']);
        save_asC=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_shift_avg.nii']);
        save_asD=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_shift_T.nii']);
        save_asE=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_shift_avg.csv']);
        
        save_asF=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_z_avg.csv']);
        save_asG=fullfile(SL.dir.save,SL.analyses{jj},[ntwk.measures{ii} '_z_T.csv']);
        
        if ~exist(fullfile(SL.dir.save,SL.analyses{jj}),'dir'), 
            mkdir(fullfile(SL.dir.save,SL.analyses{jj}));
        end
        
        if jj>N
            PC=strcmp(ntwk.contrast.pos{jj-N},SL.analyses);
            NC=strcmp(ntwk.contrast.neg{jj-N},SL.analyses);
            c=1;
            for kk=1:length(SL.dir.subjects),
                if SL.dir.include(kk)==1
                    data{jj,c}=data{PC,c}-data{NC,c};
                    dataZ{jj,c}=nan_zscore(data{PC,c})-nan_zscore(data{NC,c});
                    c=c+1;
                end
            end
        else
            c=1;
            for kk=1:length(SL.dir.subjects)
                if SL.dir.include(kk)==1
                    dataZ{jj,c}=nan_zscore(data{jj,c});
                    c=c+1;
                end
            end
        end
        
        try
        [~]=ROI_constructor(SL.region.mask,data,jj,'avg',save_asA,1);
        [~]=ROI_constructor(SL.region.mask,data,jj,'T',save_asB,1);
        
        [~]=ROI_constructor(SL.region.mask,dataZ,(jj),'avg',save_asC,1);
        [~]=ROI_constructor(SL.region.mask,dataZ,(jj),'T',save_asD,1);
        catch err
            keyboard;
        end
    end
end

%=========================================================================%
end
if SL.write.sim==1
%=========================================================================%
%% Write Similarity Profile
%=========================================================================%
% Load in Univariate measures
F=dir([SL.dir.univ filesep '*.mat']);
for ii=1:length(F)
    D=load(fullfile(SL.dir.univ,F(ii).name));
    % Silly thing in saving, end is accurate, see ROI_shell for why
    UNI{ii}=D.out_data{1};
    [~,cor,~]=fileparts(F(ii).name);
    ntwk.uni{ii}=cor; 
    clear data cor D;
end

N=length(SL.analyses);
for ii=1:length(ntwk.contrast.nam)
    SL.analyses{end+1}=ntwk.contrast.nam{ii};
end

for ii=1:length(SL.analyses),
    for jj=1:length(ntwk.measures),
        LV=fullfile(SL.dir.save,[SL.analyses{ii} '_' ntwk.measures{jj} '.csv']);
        A=excel_reader(LV);
        for kk=1:length(A), MULTI{ii,jj}(:,kk)=cell2num(A{kk}.col); end
    end
end
TH=10;
% Now we have UNI and MULTI to compare
% UNI{var}(ROI X ID)              13 ntwks
% MULTI{var X ntwk}(ROI X ID)     8  ntwks
load(fullfile(SL.dir.save,'MultiUni','Zdegree.mat')); % Load Z
load(fullfile(SL.dir.save,'MultiUni','Sdegree.mat')); % Subject data
load(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_nod_avg_degree.mat']),'con_rV');
load(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_nod_T_degree.mat']),'con_rVT');
load(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_avg_degree.mat']),'gData');
load(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_T_degree.mat']),'gDataT');
% Write Zmaps as well
% for ii=1:length(SL.analyses)
%     save_asA=fullfile(SL.dir.save,'MultiUni',['Z_' num2str(TH) '_' SL.analyses{ii} '.nii']);
%     save_asB=fullfile(SL.dir.save,'MultiUni',['Z_' SL.analyses{ii} '.nii']);
%     [~,I]=sort(Z{ii},'descend');
%     Zt{ii}=nan(size(Z{ii}));
%     Zt{ii}(I(1:TH))=Z{ii}(I(1:TH));
%     ROI_constructor(SL.region.mask,Zt,ii,'avg',save_asA,2)
%     ROI_constructor(SL.region.mask,Z,ii,'avg',save_asB,2)
% end

%=========================================================================%
% Between Networks
%=========================================================================%
clear I N;

for node_keep=[4 5 6 8]
% 4 5 6 8

aba=1;
for snet={'ant_ret' 'med_ret' 'post_ret'}
% snet='post_ret';
switch snet{aba}
    case 'ant_ret'
        mI=[1 3 7 8]; % c_coef, e_cent, c_lvge, LE
        uI=[1 2 4 5 12 13]; % ret measures
        aI=[4 5 6]; % ret retH retM
        sss=1:4;
    case 'med_ret'
        mI=[1 3 7 8]; % c_coef, e_cent, c_lvge, LE
        uI=[1 2 4 5 12 13]; % ret measures
        aI=[4 5 6]; % ret retH retM
        sss=5:6;
    case 'post_ret'
        mI=[1 3 7 8]; % c_coef, e_cent, c_lvge, LE
        uI=[1 2 4 5 12 13]; % ret measures
        aI=[4 5 6]; % ret retH retM
        sss=7:8;
end
save_net=[snet{aba} '_' num2str(node_keep)];

for TH=[5 10 20 40 80 160 459];   
    for iSubj=1:length(SL.dir.subjects)
        c1=1;
        for ii=aI
            c2=1; 
            for jj=mI
                c3=1;
                for kk=uI
                    v1=UNI{kk}(:,iSubj);  v2=MULTI{ii,jj}(:,iSubj);
                    I1=~isnan(v1);        I2=~isnan(v2);
                    I=and(I1,I2);
                    
                    N=sum(con_rV{node_keep}(sss,:));
                    N(I==0)=-inf;
                    
                    [~,Inodes]=sort(N,'descend');
                    Ikeep=zeros(1,length(N));
                    Ikeep(Inodes(1:TH))=1;
    
                    TV=[v1(and(I,Ikeep')),v2(and(I,Ikeep'))];
    %                 TV=[v1(and(I,joint_d{5}(:,iSubj)>TH)),...
    %                     v2(and(I,joint_d{5}(:,iSubj)>TH))];
                    TVr=corr(TV,'type','Spearman');
                    NET{c1}(iSubj,c3,c2)=TVr(1,2);
                    clear TV v1 v2 I1 I2 I TVr;
                    c3=c3+1;
                end
                c2=c2+1;
            end % Ntwk loop
            c1=c1+1;
        end % Analy loop
    end % Subject loop
    % NET{TimePoint}(ID X UniMeasure X MultMeasure)
    for ii=1:length(NET),
        R=squeeze(mean(NET{ii},1));
        RSA_pdm_GT(R,0.25,SL.analyses{aI(ii)},...
            fullfile(SL.dir.save,'MultiUni',[save_net '_' SL.analyses{aI(ii)} '_' num2str(TH) '.png']));
    end
    clear NET R;
end
end
aba=aba+1;
end
%=========================================================================%
end
if SL.write.con==1
%=========================================================================%
%% Write Connectivity Profile
%=========================================================================%
X=load(fullfile(SL.dir.save,'IncludedROIs'));
X.MINC=logical(X.MINC);

N=length(SL.analyses);
for ii=1:length(ntwk.contrast.nam)
    SL.analyses{end+1}=ntwk.contrast.nam{ii};
end
for ii=1:length(SL.analyses)
    display(['Running: ' SL.analyses{ii}]);
    display(datestr(now))
    for jj=1:length(SL.dir.subjects)
        % This is looking purely at degree.
        display(['  Running: ' SL.dir.subjects{jj}]);
        if ii<=N
            load_file=fullfile(SL.dir.data,SL.analyses{ii},[SL.dir.subjects{jj} '.mat']);
            load(load_file);
        else
            PC=strcmp(ntwk.contrast.pos{ii-N},SL.analyses);
            NC=strcmp(ntwk.contrast.neg{ii-N},SL.analyses);
            load_file=fullfile(SL.dir.data,SL.analyses{PC},[SL.dir.subjects{jj} '.mat']);
            L1=load(load_file);
            load_file=fullfile(SL.dir.data,SL.analyses{NC},[SL.dir.subjects{jj} '.mat']);
            L2=load(load_file);
            R1=L1.R1-L2.R1;
            R=L1.R-L2.R;
        end
        
        V=X.MINC(:,jj);
        [I,~]=find(V==1);
        
        S=sum(R); % Uncorrected Matrix
        joint_d{ii}(V==1,jj)=(S-mean(S))/(std(S));
        joint_d{ii}(V==0,jj)=NaN;

        if ntwk.con.on==1
            for kk=1:length(ntwk.con.roi)      
                [~,~,V1]=find_connections(R,I(ntwk.con.roi(kk)),ntwk.con.per);
                constr_data{kk,jj}(V==1)=zscore(V1);
                constr_data{kk,jj}(V==0)=NaN;
                clear V1
            end
        end

    end
    
    if ntwk.con.on==1
        for kk=1:length(ntwk.con.roi)
            save_asA=fullfile(SL.dir.save,SL.analyses{ii},[ntwk.con.nam{kk} '_avg.nii']);
            save_asB=fullfile(SL.dir.save,SL.analyses{ii},[ntwk.con.nam{kk} '_T.nii']);
            rV(kk,:)=ROI_constructor(SL.region.mask,constr_data,kk,'avg',save_asA,1);
            rVT(kk,:)=ROI_constructor(SL.region.mask,constr_data,kk,'T',save_asB,1);
        end
        save_asC=fullfile(SL.dir.save,SL.analyses{ii},[ntwk.con.sav '_avg.nii']);
        
        gData{ii}=sum(rV);
        gDataT{ii}=mean(rV)./(std(rV)/sqrt(length(SL.dir.subjects)));
        con_rV{ii}=rV;
        con_rVT{ii}=rVT;
        [~]=ROI_constructor(SL.region.mask,gData,ii,'avg',save_asC,2);
        rVs=sum(rV);
        I=isnan(rVs);
        rV(:,I)=[];
        R=corr(rV','type','Spearman');
        RSA_pdm_GT(R,1,SL.analyses{ii},fullfile(SL.dir.save,SL.analyses{ii},[ntwk.con.sav '_Corr.png']));
        clear constr_data rV rVT;
    end
end

% Find the most important nodes
for ii=1:length(joint_d)
    M=mean(joint_d{ii},2);
    S=std(joint_d{ii}');
    Z{ii}=(M-mean(M))./S';
end

if ~exist(fullfile(SL.dir.save,'MultiUni'),'dir'),
    mkdir(fullfile(SL.dir.save,'MultiUni'));
end

save(fullfile(SL.dir.save,'MultiUni','Sdegree.mat'),'joint_d');
save(fullfile(SL.dir.save,'MultiUni','Zdegree.mat'),'Z');
if ntwk.con.on==1
    save(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_nod_avg_degree.mat']),'con_rV');
    save(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_nod_T_degree.mat']),'con_rVT');
    save(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_avg_degree.mat']),'gData');
    save(fullfile(SL.dir.save,'MultiUni',[ntwk.con.sav '_T_degree.mat']),'gDataT');
end


return;
% RetHits
% 56 36 43 55 7
% 56: R Para
% 36: R Fusiform
% 43: L Lignugal
% 55: L Para
% 7:  L Caudate
u=mean(joint_v{8});
s=std(joint_v{8});
T=u./(s./sqrt(21));

[Tval,Tind]=sort(T,'descend');
% Most consistent changes
% 32 => F S Orb L
% 6  => Calc R.
% 63 => Post Central R
% 29 => F S L
% 77 => SM L
%=========================================================================%
end
%=========================================================================%
if SL.write.cns==1
%=========================================================================%
% 1) Load in the ROI list
SL.cns.mask='Hippocampus_R';
ntwk.contrast.pos={'RetHit'};
ntwk.contrast.neg={'RetMiss'};
ntwk.contrast.nam={'RetHvM'};

ROI_list=excel_reader(fullfile(SL.dir.root,'ROI_list.csv'));
I=strcmp(SL.cns.mask,ROI_list{2}.col); % Identify mask of interest

% Load Save Data
D=load(save_data);

% Find the nodes in the contrast
I1=strcmp(SL.analyses,ntwk.contrast.pos{1});
I2=strcmp(SL.analyses,ntwk.contrast.neg{1});
for ii=1:length(SL.dir.subjects)
    % Load in the primary contrast
    R1=load(fullfile(SL.dir.data,ntwk.contrast.pos{1},[SL.dir.subjects{ii} '.mat']));
    R2=load(fullfile(SL.dir.data,ntwk.contrast.neg{1},[SL.dir.subjects{ii} '.mat']));
    % Subtractive contrast (corrected mats)
    r1(:,ii)=R1.R1(I,:);
    r2(:,ii)=R2.R1(I,:);
    clear R1 R2;
    
    % Need degree as well
    d1(:,ii)=D.D_cent{I1,ii};
    d2(:,ii)=D.D_cent{I2,ii};
end
    
% Average Z-shift
Rt1=mean(zscore(r1)-zscore(r2),2);
Rt1s=std(zscore(r1)'-zscore(r2)');
Rt2=Rt1'./(Rt1s./sqrt(length(SL.dir.subjects)));

% Average Z-shift
Dt1=mean(zscore(d1)-zscore(d2),2);
Dt1s=std(zscore(d1)'-zscore(d2)');
Dt2=Dt1'./(Dt1s./sqrt(length(SL.dir.subjects)));

X=[Rt1,Rt2',Dt1,Dt2'];

end







