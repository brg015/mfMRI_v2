function GT_FSR()
global SL ntwk;
sdisp('FSR only works for AAL atm...',2);
for ii=1:length(ntwk.contrast.nam)
    for jj=1:length(SL.dir.subjects),
        L1=fullfile(SL.dir.data,ntwk.contrast.pos{ii},[SL.dir.subjects{jj} '.mat']);
        L2=fullfile(SL.dir.data,ntwk.contrast.neg{ii},[SL.dir.subjects{jj} '.mat']);
%         L1d=fullfile(SL.dir.data,ntwk.contrast.pos{ii},[SL.dir.subjects{jj} '_dcorr.mat']);
%         L2d=fullfile(SL.dir.data,ntwk.contrast.neg{ii},[SL.dir.subjects{jj} '_dcorr.mat']);
        La=load(L1);
        Lb=load(L2);      
        for kk=1:90,
            % FSR Calculation
            FSR(jj,kk)=-atanh(corr(La.R(:,kk),Lb.R(:,kk)));
            dFSR(jj,kk)=1-dcorr(La.R(:,kk),Lb.R(:,kk));
            % Connectivity Calculation
            Nxn{kk}(jj,:)=La.R(:,kk)-Lb.R(:,kk);
            MNxn(kk,:,jj)=La.R(:,kk)-Lb.R(:,kk);
            
            % Degree tNcalc
            Ddelt(kk,jj)=mean(La.R(:,kk))-mean(Lb.R(:,kk));
            wNxn{kk}(jj,:)=(La.R(:,kk)-Lb.R(:,kk)).*Ddelt(kk,jj);
            MwNxn(kk,:,jj)=(La.R(:,kk)-Lb.R(:,kk)).*Ddelt(kk,jj);
        end
        N1(jj,:,:)=La.R;
        N2(jj,:,:)=Lb.R;
    end
    keyboard;
    % Save
    Sdir=fullfile(SL.dir.save,'FSR');
    if ~exist(Sdir,'dir'), mkdir(Sdir); end
    
    % Stats
    for kk=1:90
        [~,~,~,stats]=ttest(Nxn{kk});
        tNxn(kk,:)=stats.tstat; clear stats;
        [~,~,~,stats]=ttest(wNxn{kk});
        twNxn(kk,:)=stats.tstat; clear stats;
    end
    
    % gotta get these stats fixed;
    [~,~,~,stats]=ttest(FSR); tFSR=stats.tstat; clear stats;
    save_as=fullfile(Sdir,[ntwk.contrast.nam{ii} '_tFSR.nii']);
    [~]=ROI_constructor(SL.region.mask,tFSR,[],'pure',save_as,3);
    
    [~,~,~,stats]=ttest(dFSR); tdFSR=stats.tstat; clear stats;
    save_as=fullfile(Sdir,[ntwk.contrast.nam{ii} '_tdFSR.nii']);
    [~]=ROI_constructor(SL.region.mask,tdFSR,[],'pure',save_as,3);
    
    [~,~,~,stats]=ttest(Ddelt'); tDdelt=stats.tstat; clear stats;
    save_as=fullfile(Sdir,[ntwk.contrast.nam{ii} '_tDegChange.nii']);
    [~]=ROI_constructor(SL.region.mask,tDdelt,[],'pure',save_as,3);
    
    csvwrite(fullfile(Sdir,[ntwk.contrast.nam{ii} '_CN.csv']),tNxn);
    csvwrite(fullfile(Sdir,[ntwk.contrast.nam{ii} '_CN.csv']),twNxn);
    
    csvwrite(fullfile(Sdir,[ntwk.contrast.nam{ii} '_FSR_raw.csv']),FSR);
    FSR_z=zscore(FSR'); FSR_z=FSR_z';
    csvwrite(fullfile(Sdir,[ntwk.contrast.nam{ii} '_FSR_z_raw.csv']),FSR);
    
    

end


