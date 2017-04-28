function GT_subnetwork()
global SL ntwk;

% Make sure only one contrast is being made
if length(ntwk.contrast.nam)>1, 
    display('Can only make one contrast at a time, try again');
    return; 
end
% Ensure a save directory is hanging out
if ~exist(ntwk.gephi.save,'dir'), mkdir(ntwk.gephi.save); end

%=========================================================================%
% Load in the network
%=========================================================================%
Nnodes=[];
for jj=1:length(ntwk.contrast.pos_R{1})
    La=load(ntwk.contrast.pos_R{1}{jj},'R');
    Lb=load(ntwk.contrast.neg_R{1}{jj},'R');
    % Get the network size
    if jj==1, Nnodes=length(La.R); end
    % Connectivity Calculation
    Nxn(jj,:,:)=La.R-Lb.R;
end
 
% Subject delta values
save1=fullfile(ntwk.gephi.save,['subj_seed_nxn.mat']);
MAT=squeeze(Nxn(:,:,ntwk.gephi.node{1}));
save(save1,'MAT');

% Stats
T=ntwk.gephi.outZ;
for jj=1:Nnodes
    for kk=1:Nnodes
        % Nxn analysis
        % 1) get nxn values
        v=squeeze(Nxn(:,jj,kk));
        % 2) outlier detection
        z=zscore(v); i1=find(abs(z)>T); v(i1)=[];
        % 3) dem stats
        [~,p,~,stats]=ttest(v);
        % 4) save that output
        Smap(jj,kk)=stats.tstat;    % T-value map
        Pmap(jj,kk)=p;              % p-value map
        % 5) save outlier information too
        ONmap(jj,kk)=length(i1);    % Number of outliers
        OVmap(jj,kk)=max(abs(z));  % Maximum value
        clear v z i1 p stats;
    end
    % Degree analysis
    v=squeeze(mean(Nxn(:,jj,:),3));
    z=zscore(v); i1=find(abs(z)>T); v(i1)=[];
    [~,p,~,stats]=ttest(v);
    Dvec(jj)=stats.tstat;     % T-value vector
    Pvec(jj)=p;               % p-value vector
    ONvec(jj)=length(i1);     % Number of outliers
    OVvec(jj)=max(abs(z));   % Maximum value
    clear stats p v z i1;
end  

save1=fullfile(ntwk.gephi.save,['seed_nxn.csv']);
E1{1}.header='ROI';
E1{1}.col=Smap(:,ntwk.gephi.node{1});
write_struct(E1,save1);

%=========================================================================%
% Save some initial output
%=========================================================================%
save1=fullfile(ntwk.gephi.save,['ONmap_z', num2str(T*10)]);
figure(1); pcolor(ONmap); colorbar; 
set(gcf,'position',[0 0 1280 1024]); export_fig(save1); close(1);

save1=fullfile(ntwk.gephi.save,['OVmap_z', num2str(T*10)]);
figure(1); pcolor(OVmap); colorbar; 
set(gcf,'position',[0 0 1280 1024]); export_fig(save1); close(1);

save1=fullfile(ntwk.gephi.save,['ONvec_z', num2str(T*10)]);
figure(1); plot(ONvec); xlabel('ROI'); ylabel('Number Outliers');
set(gcf,'position',[0 0 1280 1024]); export_fig(save1); close(1);

save1=fullfile(ntwk.gephi.save,['OVvec_z', num2str(T*10)]);
figure(1); plot(OVvec); xlabel('ROI'); ylabel('Max Value');
set(gcf,'position',[0 0 1280 1024]); export_fig(save1); close(1);
% clear outlier maps, don't need them anymore
clear ONvec OVvec ONmap OVmap;
% can get rid of nxn as well
clear Nxn;
%=========================================================================%
% Prepare output for GT_gelphi_sub
%=========================================================================%
% Load in ROI location information
switch ntwk.atlas
    case 'AAL_negX'
        N=excel_reader(fullfile(SL.root_dir,'\Geib\WFU_ROIs\ROI_locations.csv'));
    case 'HOA100'
        N=excel_reader(fullfile(SL.root_dir,'\Geib\HOA100_atlas\ROI_locations.csv'));
end

GT.node(1).name='name'; 
for ii=1:Nnodes, GT.node(1).val{ii}=['R' num2str(ii)]; end; 
GT.node(1).class='VARCHAR';

GT.node(2).name='x'; GT.node(2).val=N{8}.col; GT.node(2).class='DOUBLE';
GT.node(3).name='y'; GT.node(3).val=N{9}.col; GT.node(3).class='DOUBLE';
GT.node(4).name='z'; GT.node(4).val=N{10}.col; GT.node(4).class='DOUBLE';
GT.node(5).name='label'; GT.node(5).val=N{11}.col; GT.node(5).class='VARCHAR';

GT.node(6).name='DegreeChange'; 
for ii=1:Nnodes, GT.node(6).val{ii}=num2str(Dvec(ii)); end; 
GT.node(6).class='DOUBLE';

GT.node(7).name='AbsDegreeChange'; 
for ii=1:Nnodes, GT.node(7).val{ii}=num2str(abs(Dvec(ii))); end; 
GT.node(7).class='DOUBLE';

if isfield(ntwk.gephi,'extra')
    c=1;
    a1=excel_reader(ntwk.gephi.extra.file);
    for ii=ntwk.gephi.extra.col
        GT.node(c+7).name=a1{ii}.header{1};
        GT.node(c+7).val=cell2num(a1{ii}.col);
        GT.node(c+7).class='DOUBLE';
        c=c+1;
    end
end

for ii=1:length(ntwk.gephi.node)
    fname=['Network_z' num2str(T*10) '_seed_' num2str(ntwk.gephi.node{ii})];
    %=====================================================================%
    % Integration and Segregation Properites
    %=====================================================================%
    Nlim=sum(Smap(:,ntwk.gephi.node{ii})>ntwk.gephi.thresh);
    
    % 2) for each of these nodes, identify their most connected
    % constinuetns
    for jj=1:Nnodes
        % 3) strength sort
        [y1,i1]=sort(Smap(:,jj),'descend');
        % 4a) calc by number
        i2=i1(2:Nlim+1);
        % 5b) by threshold
        i3=i1(y1>ntwk.gephi.thresh);
        % 6) for each measure computer interconnectedness and degree
        % *degree network
        Sdeg2=Smap([jj i2'],:);
        Sdeg2(:,[jj i2'])=[];
        % *tri network
        Stri2=Smap(:,[jj i2']); Stri2=Stri2([jj i2'],:);
        u2_deg(jj)=mean(mean(Sdeg2));
        u2_tri(jj)=nanmean(nanmean(Stri2));
        
        Sdeg3=Smap([jj i3'],:);
        Sdeg3(:,[jj i3'])=[];
        % *tri network
        Stri3=Smap(:,[jj i3']); Stri3=Stri3([jj i3'],:);
        u3_deg(jj)=mean(mean(Sdeg3));
        u3_tri(jj)=nanmean(nanmean(Stri3)); 
        clear y1 i1 i2 i3 Sdeg2 Stri2 Sdeg3 Stri3;
    end
    
    figure(1);
    temp1=u2_deg(ntwk.gephi.node{ii});
    temp2=sort(u2_deg); pos=find(temp2==temp1);
    subplot(2,2,1); hist(u2_deg); title(['N match deg=' num2str(pos)]); 
    clear temp1 temp2 pos;
    temp1=u2_tri(ntwk.gephi.node{ii});
    temp2=sort(u2_tri); pos=find(temp2==temp1);
    subplot(2,2,2); hist(u2_tri); title(['N match tri=' num2str(pos)]);
    clear temp1 temp2 pos;
    temp1=u3_deg(ntwk.gephi.node{ii});
    temp2=sort(u3_deg); pos=find(temp2==temp1);
    subplot(2,2,3); hist(u3_deg); title(['T match deg=' num2str(pos)]);
    clear temp1 temp2 pos;
    temp1=u3_tri(ntwk.gephi.node{ii});
    temp2=sort(u3_tri); pos=find(temp2==temp1);
    subplot(2,2,4); hist(u3_tri); title(['T match tri=' num2str(pos)]);
    clear temp1 temp2 pos;
    save1=fullfile(ntwk.gephi.save,[fname '_hist_dist.png']);
    set(gcf,'position',[0 0 1280 1024]); export_fig(save1); close(1); 
    
    %=====================================================================%
    % Create gephi networks
    %=====================================================================%
    GT.save_to=fullfile(ntwk.gephi.save,[fname '.gdf']);
    L=GT_gelphi_sub(ntwk.gephi.node{ii},Smap,ntwk.gephi.thresh,GT);
    save1=fullfile(ntwk.gephi.save,[fname '.png']);
    set(gcf,'position',[0 0 1280 1024]); export_fig(save1); close(1); 
    save1=fullfile(ntwk.gephi.save,[fname '.csv']);
    M{1}.header='labels'; M{1}.col=L.labels;
    M{2}.header='ring';   M{2}.col=L.ring;
    write_struct(M,save1);
    clear M;
end

   
