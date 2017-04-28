function RSA_noise(SL)
%=========================================================================%
% Description: Though specialized for ER match, I think this could be used
% for other paradigms as well. It's a tad tough to utilize across the board
% in memory studies as the RDM are variable across subjects. Thus while
% encoding can typically be properly managed. Retrieval sets are too
% variable across subjects to be overly useful.
%=========================================================================%
N=96;

for ii=1:length(SL.noise)
    %=====================================================================%
    % Process the incoming data
    %=====================================================================%
    % 1) Make all of the RDMs equal to the size of the minimum required
    % such that all subjects are comparable and have the maximal number of
    % possible trails
    
    % Special calc. cause I forgot to grab ID idxs
    for cursub=1:length(SL.noise(ii).subjects)
        cd(strcat(SL.dir.stpath,'\',SL.dir.subjects{cursub}));
        load('SPM.mat');
        SL.design.ID_idx=[];
        SL=RSA_find_betas(SPM,SL,...
            strcat(SL.dir.stpath,'\',SL.dir.subjects{cursub}),...
            SL.dir.subjects{cursub});
        % Finding missing trials
        ID_idx{cursub}=setdiff(1:N,SL.design.ID_idx(1:length(SL.design.ID_idx)/2));
        % Need to know what the subject has included (and order)
        ID{cursub}=SL.design.ID_idx(1:length(SL.design.ID_idx));
    end
    BAD_trials=unique([ID_idx{:}]);
    
    for cursub=1:length(SL.noise(ii).subjects)
        % Psuedo design
        D=nan(length(ID{cursub}));
        D(1:length(D)/2,1:length(D)/2)=1;
        % Bad trial
        [~,I_bt,~]=intersect(ID{cursub}(1:end/2),BAD_trials);    
        for kk=1:length(I_bt)
            D(I_bt(kk),:)=0;
            D(:,I_bt(kk),:)=0;
        end
        % Off diagnol
        nan_idx=tril(ones(length(D)),-1)==0;
        D(nan_idx)=nan;
        D(length(D)/2+1:end,:)=nan;
        D(:,length(D)/2+1:end)=nan;
        % Format vector
        rep_vector{cursub}=D(~isnan(D(:)));
    end
    keyboard
    
    %===========================%
    % Replotting Vector
    %===========================%
    % This should have a size equal to that of rep_vector as is this where
    % values are going to be saved into.
    L=(N-length(I_bt));
    Dtemp=ones(L);
    nan_idx=tril(ones(L),-1)==0;
    Dtemp(nan_idx)=nan;
       
    for jj=1:length(SL.noise(ii).subjects)
        if SL.noise(ii).ROI==0
            R{jj}=load(fullfile(SL.noise(ii).dir,SL.noise(ii).subjects{jj},SL.noise(ii).file));
            if jj==1
                % Need to make a mask based upon the 1st subject in order to
                % make memory managable (even on serv1) -> ~2% per subj
                X=mean(R{jj}.NOISE,2);
                I=~isnan(X);
            end
            RR(jj,:,:)=R{jj}.NOISE(I,logical(rep_vector{jj}));
            R{jj}.NOISE=[];
        else
            R{jj}=load(fullfile(SL.noise(ii).dir,SL.noise(ii).subjects{jj},SL.noise(ii).file));
            % Uses for design matrix to operate on
            R{jj}.R(isnan(R{jj}.design.matrix{1}))=nan;
            V=R{jj}.R(~isnan(R{jj}.R));
            RR(jj,:)=V(logical(rep_vector{jj}));
            % Subject specific maps
            Dsubj{jj}=Dtemp;
            Dsubj{jj}(Dsubj{jj}==1)=RR(jj,:);
            if jj==1
                Template=Dtemp;
                DM=R{jj}.design.matrix{1}(~isnan(R{jj}.design.matrix{1}));
                DM=DM(logical(rep_vector{jj}));
                Template(Dtemp==1)=DM;
            end    
        end
    end
    
    % Hate to use a for loop, but don't think matrix operations have enough
    % memory to work well here
    m=0; 
    for jj=1:size(RR,2), 
       if SL.noise(ii).ROI==0, 
           V=squeeze(RR(:,jj,:)); 
       else
           V=RR;
           Result=Dtemp;
           Result(Dtemp==1)=mean(V,1);
       end
       
       mV=mean(V,1);
       if sum(isnan(mV))==0
           % I don't think the mask is quite as good as I thought it was,
           % think this is causing NaN erros, trying to avoid this, but
           % will have to see how many voxels we loose.
           for kk=1:size(V,1)
               % Mean not inlcuding the trial tested
               tmV=mean(V(setdiff(1:size(V,1),kk),:));
               corr_high(kk)=corr(mV(:),V(kk,:)','type','Spearman');
               corr_low(kk)=corr(tmV(:),V(kk,:)','type','Spearman');
               if SL.noise(ii).ROI==1
                   corr_mod(kk)=corr(DM,V(kk,:)','type','Spearman');
               end
           end
           LV=mean(corr_low);
           HV=mean(corr_high);
       else
           m=m+1;
           LV=nan;
           HV=nan;
       end
       out_lin_low(jj)=LV;
       out_lin_hgh(jj)=HV;
       
       if SL.noise(ii).ROI==1, break; end
       
    end

    

    if SL.noise(ii).write==1
        out_low=nan(size(I));
        out_hgh=nan(size(I));
        out_low(I)=out_lin_low;
        out_hgh(I)=out_lin_hgh;
        
        SL.V=spm_vol(SL.design.ID_file{1});
        % Now just need to reshape and save
        [~,name_save,~]=fileparts(SL.noise(ii).file);
        for jj=1:2
            switch jj
                case 1, data=out_low; save_as=[name_save '_low'];
                case 2, data=out_hgh; save_as=[name_save '_high'];
            end
            data = reshape(data,R{1}.DIM);
            SL.V.dim=R{1}.DIM;
            SL.V.fname = fullfile(SL.dir.outpath,[save_as '.img']);
            SL.V.descrip = name_save;
            spm_write_vol(SL.V,data);
        end
    end
    keyboard
end
            
        
            