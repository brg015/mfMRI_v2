
% BRG Spring 2014
% [Box,POS_idx]=RSA_find_betas(SPM,SL)
% RSA specific
%
% Description
%   Inputs:
%       SPM : Loaded SPM stucture
%       SL  : SL structure (see master file for details)
%       data_dir : location of subject data
%       subject : subject ID
%   Outputs:
%       Box     (vector)->Indicates size of conditioned cells
%
function [SPM]=RSA_find_betas(SPM,data_dir,subject)
global SL;

debug=1; ignore_diff=0; C=[];
%========================%
% RETURNS SL (only!)
%========================%
% SL.design.ID_idx (matrix)         => ID Code
% SL.design.ID_descrip (cell array)
% SL.design.ID_file (cell array)
% SL.design.Box (Dimensions)
%=========================================================================%
% Mumford Hack
%=========================================================================%
% SPM is [], we need it fix it to look like it actually exists. Pretty sure
% all we will need is to set up the Vbeta descrip field
if SL.design.mumford==1,
    % Collect file list
    SPM_files=dir(fullfile(SL.dir.stpath,subject,SL.design.mumford_dir,['*' SL.design.mumford_str]));
    if isempty(SPM_files),
        display('Error in RSA_find_betas');
        error('Can''t find mumford files');
    end
    for ii=1:length(SPM_files),
        [~,d,~]=fileparts(SPM_files(ii).name);
        SPM.Vbeta(ii).descrip=d;
        SPM.Vbeta(ii).fname=fullfile(SL.design.mumford_dir,SPM_files(ii).name);
        clear d;
    end
end
%=========================================================================%
% Sort data files
%=========================================================================%
% Delineate which files need to be loaded & initialize structures
ID_code=nan(length(SPM.Vbeta),1);
descrip_match=nan(length(SPM.Vbeta),length(SL.design.cond_str));

%========================%
% Collect Codes
%========================%
for ii=1:length(SPM.Vbeta)
    % Pull unique Vbeta IDs if possible
    if SL.ID.check>0
        % Pull ID code for each file
        indx1=findstr(SL.ID.match{1},SPM.Vbeta(ii).descrip);
        indx2=findstr(SL.ID.match{2},SPM.Vbeta(ii).descrip);
        if ~isempty([indx1,indx2])
            try
                ID_code(ii)=str2double(SPM.Vbeta(ii).descrip(indx1+length(SL.ID.match{1}):indx2-1));
            catch err
                ID_code(ii)=NaN;
            end
        end
    end
    % Also find description matches
    for jj=1:length(SL.design.cond_str)
        descrip_match(ii,jj)=mstrcmp(SPM.Vbeta(ii).descrip,SL.design.cond_str{jj});         
    end
end

%========================%
% DRM Hack - too hard to integrate
%========================%
if SL.ID.check==2
    RSA_find_betas_hemi(descrip_match,ID_code,SPM,subject,data_dir);
    % Some quick visualization stuff here as well...
    % -- find FF counts -- % add to SL.design
    display('Warning: Supposing on inputs...');
    SL.design.hemi_trial_counts_cfmh=zeros(length(SL.design.ID_idx),4);
    % 1 5: CR
    % 2 6: FA 
    % 3 7: M
    % 4 8: H
    jb=cumsum(SL.design.Box);
    for ii=1:72
        I=(SL.design.ID_idx==ii); % Index list
        % Stats please 
        It=zeros(1,length(SL.design.ID_idx));
        Ic=It; Ic(1:jb(1))=1;
        If=It; If(jb(1)+1:jb(2))=1;
        Im=It; Im(jb(2)+1:jb(3))=1;
        Ih=It; Ih(jb(3)+1:jb(4))=1;
        % Compute them
        C=sum(and(Ic,I));
        F=sum(and(If,I));
        M=sum(and(Im,I));
        H=sum(and(Ih,I));
        % Save it away
        SL.design.hemi_trial_counts_cfmh(I,1)=C;
        SL.design.hemi_trial_counts_cfmh(I,2)=F;
        SL.design.hemi_trial_counts_cfmh(I,3)=M;
        SL.design.hemi_trial_counts_cfmh(I,4)=H;
        SL.design.hemi_trial_cfmh(ii,:)=[C,F,M,H];
    end % for each list
    return;
end

% Determine position ind. of ID codes
c=1; POS_idx=nan(length(SPM.Vbeta),1);
Box=nan(1,4);
for ii=1:size(descrip_match,2),
    Box(ii)=length(find(descrip_match(:,ii)==1));
    location=c:c+Box(ii)-1;
    POS_idx(descrip_match(:,ii)==1)=location;
    c=c+Box(ii);
end
SL.design.Box=Box;
% SL.design.Box is updated later on to account for overlap
%========================%
% Classify
%========================%
if SL.ID.check>0 && isfield(SL.ID,'overlap')
    % Check to make sure that overlap occurs
    for ii=1:length(SL.ID.overlap)
        % Location of description matches
        ID_temp_loc=descrip_match(:,SL.ID.overlap{ii}==1)==1;   
        ID_align_loc=descrip_match(:,SL.ID.overlap{ii}==2)==1;

        % Pull relevent code values
        ID_temp=ID_code(ID_temp_loc);   % First set (enc oft)
        ID_align=ID_code(ID_align_loc); % Second set (ret oft)
        
        % IDs should overlap, if don't error is present
        switch SL.ID.setcheck(ii)
            case 1 % Max Overlap
                ignore_diff=1; 
                % Save C which indicates overlapping trials. File output
                % will be trimmed to only include items in C
                [tmp,~,~]=intersect(ID_align,ID_temp);
                C=[C; tmp];
        end
        
        % Sort IDs numerically
        if isempty(setdiff(ID_temp,ID_align)) || ignore_diff==1
            % Now we can just sort the values
            [~,sort_temp]=sort(ID_temp);
            [~,sort_align]=sort(ID_align);
            % And modify the offsets based upon the ordering - horribly
            % ugly, but gets the job done...
            POS_values_temp=POS_idx(ID_temp_loc);   val_temp=nan(1,length(sort_temp));
            POS_values_align=POS_idx(ID_align_loc); val_align=nan(1,length(sort_align));
            for jj=1:length(sort_temp) 
                val_temp(sort_temp(jj))=POS_values_temp(jj);  
            end
            for jj=1:length(sort_align)
                val_align(sort_align(jj))=POS_values_align(jj);
            end
            POS_idx(ID_temp_loc)=val_temp;
            POS_idx(ID_align_loc)=val_align;   
        else
            if debug==1, 
                display('Confusion in RSA_find_betas.m - bad overlap');
                keyboard; 
                SPM.Vbeta(ID_align_loc).descrip
            end
            error('Badly defined overlap');
        end
    end
end
% KeyVar
%   POS_idx
%   descrip_match
% Display sorted files
% keyboard;
c=1; 
for ii=1:size(descrip_match,2)
    idx=find(descrip_match(:,ii)==1);
    [~,pos_sort]=sort(POS_idx(idx));
    bc=1; 
    for jj=pos_sort'
        POS=idx(jj);
        if ~isempty(find(C==ID_code(POS), 1)) || ignore_diff==0
            SL.design.ID_idx(c)=ID_code(POS);
            SL.design.ID_descrip{c}=SPM.Vbeta(POS).descrip;
            SL.design.ID_file{c,1}=fullfile(data_dir,SPM.Vbeta(POS).fname);    
            
            if SL.run.include>=1
                indx1=findstr(SL.run.match{1},SL.design.ID_descrip{c});
                indx2=findstr(SL.run.match{2},SL.design.ID_descrip{c});
                SL.design.run(c)=str2num(SL.design.ID_descrip{c}(indx1(end)+length(SL.run.match{1}):indx2(end)-1));
            end

            if SL.regress.on==1
               for kk=1:length(SL.regress.str)
                   indx1=findstr(SL.regress.str{kk}{1},SL.design.ID_descrip{c});
                   indx2=findstr(SL.regress.str{kk}{2},SL.design.ID_descrip{c});
                   SL.regress.val(c,kk)=str2num(SL.design.ID_descrip{c}(indx1(end)+length(SL.regress.str{kk}{1}):indx2(end)-1));
               end
            end  
            B(ii)=bc; bc=bc+1; c=c+1; 
        end
    end
end

if SL.run.include>=1
   for ii=1:length(SL.design.run)
       for jj=1:length(SL.design.run)
           SL.run.matrix(ii,jj)=(SL.design.run(ii)==SL.design.run(jj));
       end
   end
end

% Needs to fix up regression parameters here
cBox=cumsum(B);
if SL.regress.on==1
    if SL.regress.match>0 
        for ii=1:length(SL.ID.overlap)
            r1I=find(SL.ID.overlap{ii}==1);
            r2I=find(SL.ID.overlap{ii}==2);
            if r1I==1, r1=1:cBox(1); else r1=cBox(r1I-1)+1:cBox(r1I); end
            if r2I==1, r2=1:cBox(1); else r2=cBox(r2I-1)+1:cBox(r2I); end
            if SL.regress.match==1
                SL.regress.val(r2,:)=SL.regress.val(r1,:);
            elseif SL.regress.match==2
                SL.regress.val(r1,:)=SL.regress.val(r2,:);
            end
        end
    end
end

if ~exist('B','var'),
    display('....FAILED: no files exist'); 
    SL.err=1; return;
end
SL.design.Box=B; % Update Box to include only added values
%=========================================================================%
%% Look at preprocessed values, not beta values
%=========================================================================%
% This section simply examines the SPM description and grabs preprocessed
% images corresponding to the beta files
if SL.preprocess.on==1
    c1=1;
    for ii=1:length(SPM.Sess)
        for jj=1:length(SPM.Sess(ii).U)
            event{c1}=SPM.Sess(ii).U(jj).name{1};
            event_time(c1)=SPM.Sess(ii).U(jj).ons;
            c1=c1+1;
        end
    end

    for ii=1:length(event)
        % Find correct index
        index=strsearch(event{ii},SL.design.ID_descrip);
        if ~isempty(index)
            % Find which session
            sess_idx=strfind(SL.design.ID_descrip{index},'Sn(');
            sess=str2double(SL.design.ID_descrip{index}(sess_idx+3));
            % Find index of event 
            for jj=1:length(SL.preprocess.advance.N)
                FileN=round(event_time(ii)+SL.preprocess.advance.N(jj));
                % Define new file to read
                SL.design.ID_file{index,jj}=fullfile(SL.preprocess.dir,subject,...
                    SL.preprocess.sub_dir{sess},...
                    [SL.preprocess.prefix.name n2sp(FileN,SL.preprocess.prefix.num) ...
                    SL.preprocess.prefix.img]);
            end
        end
    end
end


    
    
    
    
    
    
    
    
    