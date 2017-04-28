function GT_compute_measures(save_data)
global ntwk SL;
%=========================================================================%
%% Loop through different network types e.g. encoding & retrieval
%=========================================================================%
keyboard;
for ii=1:length(SL.analyses)
    display(['Running: ' SL.analyses{ii}]);
    display(datestr(now))
    for jj=1:length(SL.dir.subjects)
        % If temp files exist, load them, and then continue, turned off for
        % the time being...
%         if ntwk.Niterations==1
%             if (exist(fullfile(SL.dir.save,'temp',[ntwk.save_name '_' num2str(ii) '_' num2str(jj) '.mat']),'file') ...
%                     && SL.dir.overwrite==0);
%                 load(fullfile(SL.dir.save,'temp',[ntwk.save_name '_' num2str(ii) '_' num2str(jj) '.mat'])); continue;
%             end
%         end
        
        % Load in R and R1
        display(['  Running: ' SL.dir.subjects{jj}]);
        load_file=fullfile(SL.dir.data,SL.analyses{ii},[SL.dir.subjects{jj} '.mat']);
        load(load_file);
        % rm indices 39 40 i.e. AAL hc
        if ntwk.rm_hc==1, 
            R1(39:40,:)=[]; R1(:,39:40)=[]; 
            R(39:40,:)=[]; R(:,39:40)=[]; 
        elseif ntwk.rm_hc==2,
            R1(40,:)=[]; R1(:,40)=[]; 
            R(40,:)=[]; R(:,40)=[]; 
        end
        
        switch ntwk.analysis
            case 'rand_effc'
                % Rand1=randmio_und_connected(R1, 10);
                Rand1 = null_model_und_sign(R1); % Changed 12/10
                G_effc_rand{ii,jj}=efficiency_wei_node(Rand1).^(-1);  % changed 12/10 => char. path
                if ntwk.Niterations==1
                    save(fullfile(SL.dir.save,SL.analyses{ii},[SL.dir.subjects{jj} '_rand_effc' ntwk.str '.mat']),'Rand1');
                end
                clear Rand1
            case 'rand'
                % Rand1=randmio_und_connected(R1, 10);
                Rand1 = null_model_und_sign(R1); % Changed 12/10
                C_coef_rand{ii,jj}=clustering_coef_wu(Rand1); 
                G_effc_rand{ii,jj}=efficiency_wei_node(Rand1).^(-1);  % changed 12/10 => char. path
                if ntwk.Niterations==1
                    save(fullfile(SL.dir.save,SL.analyses{ii},[SL.dir.subjects{jj} '_rand' ntwk.str '.mat']),'Rand1');
                end
                clear Rand1
            case 'rand_lat'
                Rand2=latmio_und_connected(R1,10);
                C_coef_rand_lat{ii,jj}=clustering_coef_wu(Rand2);
                G_effc_rand_lat{ii,jj}=efficiency_wei(Rand2); 
                save(fullfile(SL.dir.save,SL.analyses{ii},[SL.dir.subjects{jj} '_latrand' ntwk.str '.mat']),'Rand2');
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
            case 'path_length'
                % C_coef{ii,jj}=clustering_coef_wu(R1); 
                PL{ii,jj}=efficiency_wei_node(R1).^-1; 
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
        clear R R1; % Clear original network
 
        % Save progress as we go, cause it seems to like to crash, seems a
        % bit excessive to do this with every single loop for random
        % networks
        if ntwk.Niterations==1
            save(fullfile(SL.dir.save,'temp',[ntwk.save_name '_' num2str(ii) '_' num2str(jj) ' ntwk.str ']),ntwk.measures{:});
        end
    end % Subject Loop
    display(datestr(now))
end % Analy Loop

save(save_data,ntwk.measures{:});
%=========================================================================%