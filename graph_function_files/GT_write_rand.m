function GT_write_rand()

global SL;
global ntwk;
for ii=1:ntwk.Niterations
    LV=fullfile(SL.dir.save,'RandNetworks',...
        [ntwk.save_name '_N' num2str(ii) ntwk.str '.mat']);
    RND=load(LV);
    fn=fieldnames(RND);
    for jj=1:size(RND.(fn{1}),1)
        for kk=1:size(RND.(fn{1}),2)
            % Load constructed networks
            switch ntwk.rand
                case 1
                    % Compute average over all nodes
                    PL(ii,jj,kk)=mean(RND.G_effc_rand{jj,kk});
                case 2
                    Clust(ii,jj,kk)=mean(RND.C_coef_rand{jj,kk});
                case 3
                    Clust(ii,jj,kk)=mean(RND.C_coef_rand{jj,kk}); % Average across subject nodes
                    PL(ii,jj,kk)=mean(RND.G_effc_rand{jj,kk});    % Average across subject nodes
            end
        end % Subject loop
    end % Network loop
    clear RND;
end

% PL [Iterations X Network X Subject]
switch ntwk.rand
    case 1
        mPL=squeeze(mean(PL,1)); % Average over iterations
        for ii=1:size(mPL,1)
            for jj=1:size(mPL,2)
                Pdata{jj}.header=SL.dir.subjects{jj};
                Pdata{jj}.col=mPL(ii,jj);
            end
            write_struct(Pdata,fullfile(SL.dir.save,SL.analyses{ii},...
                ['RandPathLength' ntwk.str '.csv']));
            clear Pdata;
        end
        mmPL=squeeze(mean(mPL,2))
        msPL=squeeze(std(mPL,0,2))
    case 2
        mClust=squeeze(mean(Clust,1));
        for ii=1:size(mClust,1)
            for jj=1:size(mClust,2)
                Cdata{jj}.header=SL.dir.subjects{jj};
                Cdata{jj}.col=mClust(ii,jj);
            end
        write_struct(Cdata,fullfile(SL.dir.save,SL.analyses{ii},...
            ['RandCluster' ntwk.str '.csv']));
        clear Cdata;
        end       
        mmClust=squeeze(mean(mClust,2))
        mmClust=squeeze(std(mClust,0,2))
    case 3
        mClust=squeeze(mean(Clust,1));
        mPL=squeeze(mean(PL,1));
        for ii=1:size(mClust,1)
            for jj=1:size(mClust,2)
                Cdata{jj}.header=SL.dir.subjects{jj};
                Cdata{jj}.col=mClust(ii,jj);
                Pdata{jj}.header=SL.dir.subjects{jj};
                Pdata{jj}.col=mPL(ii,jj);
            end
        write_struct(Pdata,fullfile(SL.dir.save,SL.analyses{ii},...
            ['RandPathLength' ntwk.str '.csv']));
        write_struct(Cdata,fullfile(SL.dir.save,SL.analyses{ii},...
            ['RandCluster' ntwk.str '.csv']));
        clear Cdata Pdata;
        end
        display('Path (mean/std');
        mmPL=squeeze(mean(mPL,2))
        msPL=squeeze(std(mPL,0,2))
        display('Clust (mean/std)');
        mmClust=squeeze(mean(mClust,2))
        mmClust=squeeze(std(mClust,0,2))
end



