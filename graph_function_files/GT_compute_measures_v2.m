function GT_compute_measures_v2(ntwk_files,GT_I,save_dir_pos)
global ntwk SL;
%=========================================================================%
%% Compute measures of interest
%=========================================================================%
% These should've been pre-arranged to directly aline with subjects, but
% this is by no means a gurantee that this was eactually done.
display(datestr(now));
for ii=1:length(ntwk_files)
    load(ntwk_files{ii},'R'); % Load in variable R
    R1=abs(R); R1(logical(eye(size(R))))=0; R1(R1<0)=0;
    % R => original matrix
    % R1=> zero diag + nothing <0
    if SL.transform==1, R1=R1.^4; end
        
    % Remove the Hc if analysis calls for such shenanigans
    if ntwk.rm_hc==1, 
        R1(39:40,:)=[]; R1(:,39:40)=[]; 
        R(39:40,:)=[]; R(:,39:40)=[]; 
    elseif ntwk.rm_hc==2,
        R1(40,:)=[]; R1(:,40)=[]; 
        R(40,:)=[]; R(:,40)=[]; 
    end
    
    % Can be tricky as some of the GT toolbox requires correlation based
    % matrices where other functions require "Length" based matrices i.e.
    % the inverse
    switch ntwk.analysis{GT_I}
        case 'Eig'
            MAT(ii,:)=eigenvector_centrality_und(R1);
        case 'C_Coef', 
            MAT(ii,:)=clustering_coef_wu(R1); 
        case 'Degree', 
            % ignores negative weights
            MAT(ii,:)=sum(R1); 
        case 'wDegree', 
            % ignores negative weights
            MAT(ii,:)=sum(R); 
        case 'Path_L'
            MAT(ii,:)=mean(distance_wei(1./R1));
        case 'Global_Effc', 
            MAT(ii,:)=efficiency_wei_node(R1); 
        case 'wMod1', 
            % Only save out those Q values for now
            for kk=1:10
                [Ci Q]=modularity_louvain_und_sign(R,'sta');
                Qt(kk)=Q;
            end
            MAT(ii)=max(Qt);
         case 'Mod1', 
            for kk=1:10
                [Ci Q]=modularity_louvain_und(R1,1);
                Qt(kk)=Q;
            end
            MAT(ii)=max(Qt);
        case 'Local_Effc',
            MAT(ii,:)=efficiency_wei_node(R1,1);
        case 'BW_cent'
            MAT(ii,:)=betweenness_wei(1./R1);
        case 'PR_cent'
            MAT(ii,:)=pagerank_centrality(R1, 0.85);
        case 'FSR'
            % This needs handled later on...
            continue;   
        otherwise
            display('Analysis DNE atm - sorry');
            return;
    end
end
display(datestr(now));

% Save some sort of .mat file here and output some basic maps as well
save(fullfile(save_dir_pos,[ntwk.analysis{GT_I} '.mat']),'MAT');
save_asA=fullfile(save_dir_pos,[ntwk.analysis{GT_I} '_avg.img']);
if ~(strcmp(ntwk.analysis{GT_I},'Mod1') || strcmp(ntwk.analysis{GT_I},'wMod1'))
    [~]=ROI_constructor(SL.region.mask,mean(MAT),'','pure',save_asA,3);    
end
%=========================================================================%










