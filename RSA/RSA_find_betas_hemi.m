function RSA_find_betas_hemi(descrip_match,ID_code,SPM,subject,data_dir)  
global SL;
% Step 1) Define Box, it tells us what is included
for ii=1:length(SL.ID.overlap)
    % Location of description matches
    ID_temp_loc=descrip_match(:,SL.ID.overlap{ii}==1)==1;   
    ID_align_loc=descrip_match(:,SL.ID.overlap{ii}==2)==1;

    % Pull relevent code values
    ID_temp=ID_code(ID_temp_loc);   % First set (enc oft)
    ID_align=ID_code(ID_align_loc); % Second set (ret oft)

    % Ensure overlap
    uID_temp=unique(ID_temp);
    uID_align=unique(ID_align);
    C=setdiff(uID_align,uID_temp);
    for jj=1:length(C)
        ID_temp(ID_temp==C(jj))=[];
        ID_temp_loc(ID_temp==C(jj))=[];
        ID_align(ID_align==C(jj))=[];
        ID_align_loc(ID_align==C(jj))=[];
    end

    % Design Box based upon longest length
    % Box(SL.ID.overlap{ii}==1)=max([length(ID_temp),length(ID_align)]);
    % Box(SL.ID.overlap{ii}==2)=max([length(ID_temp),length(ID_align)]);
    % ID_align is ret, we only need to line up encoding with retrieval
    % trials, if anything, above should've been min...
    Box(SL.ID.overlap{ii}==1)=length(ID_align);
    Box(SL.ID.overlap{ii}==2)=length(ID_align);
end
BC=cumsum(Box);

% Step 2) Define matches based upon longest
for ii=1:length(SL.ID.overlap)
    % Location of description matches
    ID_temp_loc=descrip_match(:,SL.ID.overlap{ii}==1)==1;   
    ID_align_loc=descrip_match(:,SL.ID.overlap{ii}==2)==1;

    % Pull relevent code values
    ID_temp=ID_code(ID_temp_loc);   % First set (enc oft)
    ID_align=ID_code(ID_align_loc); % Second set (ret oft)

    uID_temp=unique(ID_temp);
    uID_align=unique(ID_align);
    C=setdiff(uID_align,uID_temp);
    for jj=1:length(C)
        ID_temp(ID_temp==C(jj))=[];
        ID_temp_loc(ID_temp==C(jj))=[];
        ID_align(ID_align==C(jj))=[];
        ID_align_loc(ID_align==C(jj))=[];
    end

    Tcode=ID_temp;
    Tp=find(ID_temp_loc);

    Mcode=ID_align;
    Mp=find(ID_align_loc);

    L=length(ID_align);

    i1=BC(SL.ID.overlap{ii}==1)-Box(SL.ID.overlap{ii}==1);
    i2=BC(SL.ID.overlap{ii}==2)-Box(SL.ID.overlap{ii}==1);

    % Using the Mcode (match code from retrieval) we need to identify
    % items to include
    for jj=1:L, 
        % First let's identify included codes...
        idx=find(Tcode==Mcode(jj));
        if ~isempty(idx)
            POS=Tp(idx);
            SL.design.ID_idx(jj+i1)=Tcode(idx);
            SL.design.ID_descrip{jj+i1}=SPM.Vbeta(POS).descrip;
            SL.design.ID_file{jj+i1,1}=fullfile(data_dir,SPM.Vbeta(POS).fname);

            % Lay down template code (ret)
            SL.design.ID_idx(jj+i2)=Mcode(jj);
            SL.design.ID_descrip{jj+i2}=SPM.Vbeta(Mp(jj)).descrip;
            SL.design.ID_file{jj+i2,1}=fullfile(data_dir,SPM.Vbeta(Mp(jj)).fname);
        end 
    end

end
SL.design.Box=Box;
write_vbeta(subject);

% Set regress values to be equal to the condition
for ii=1:length(BC)
    if ii==1, SL.regress.val(1:BC(1))=ii; else SL.regress.val(BC(ii-1)+1:BC(ii))=ii; end
end
SL.regress.val=SL.regress.val';
