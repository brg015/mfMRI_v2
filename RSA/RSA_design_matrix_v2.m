% BRG Spring 2014
% [Box,POS_idx]=RSA_design_matrix(SPM,SL)
% RSA specific
%
% Description
%   Inputs:
%       Box (vector) : from RSA_find_betas.m output
%       SL           : SL structure (see master file for details)
%       indx         : index of design matrix
%   Outputs:
%       SL.design.matrix
%
function RSA_design_matrix_v2(indx)
global SL;
%=========================================================================%
% Generic Model Setup
%=========================================================================%
Box=SL.design.Box;
% If no interactions or weights are setup, include all of them, set 
if isempty(SL.design.interactions{indx}),
    SL.design.interactions{indx}=combn(1:length(SL.design.cond_str),2);
end
% Setup the design matrix (items to include*)
SL.design.matrix{indx}=nan(sum(Box));

if SL.design.custom(indx)==1        
    load(SL.design.model{indx}); % Models should be saved as R vars        
        % Based upon the items saved in the model, we adjust the
        % model to align with the subjects data and remove data
        % that is not included for the given subject
    if ~exist('stim_ID_num','var') || ~exist('R','var'),
        keyboard;
    end  
end
cBox=cumsum(Box);
%=========================================================================%
% Set AtoB inputs if needed
%=========================================================================%
if strcmp(SL.design.calc{indx},'AtoB')
    A=SL.design.E{indx}; Av=[];
    for ii=1:length(A)
        if A(ii)==1, Av=[Av,1:cBox(1)]; else Av=[Av,cBox(A(ii)-1)+1:cBox(A(ii))]; end
    end     
    B=SL.design.R{indx}; Bv=[];
    for ii=1:length(B)
        if B(ii)==1, Bv=[Bv,1:cBox(1)]; else Bv=[Bv,cBox(B(ii)-1)+1:cBox(B(ii))]; end
    end  
    SL.design.Ev{indx}=Av;
    SL.design.Rv{indx}=Bv;
end
clear Av Bv A B;
%=========================================================================%
% Set regress inputs if needed (for identity1)
%=========================================================================%
if strcmp(SL.design.calc{indx},'Identity1')
    if SL.design.anova{indx}.regress.on==1
        if SL.design.anova{indx}.regress.E.uni==1
            A=SL.design.anova{indx}.regress.E.pos; Av=[];
            for ii=1:length(A)
                if A(ii)==1, Av=[Av,1:cBox(1)]; else Av=[Av,cBox(A(ii)-1)+1:cBox(A(ii))]; end
            end  
            SL.design.anova{indx}.regress.E.v=Av; clear A Av;
        end
        if SL.design.anova{indx}.regress.R.uni==1
            A=SL.design.anova{indx}.regress.R.pos; Av=[];
            for ii=1:length(A)
                if A(ii)==1, Av=[Av,1:cBox(1)]; else Av=[Av,cBox(A(ii)-1)+1:cBox(A(ii))]; end
            end 
            SL.design.anova{indx}.regress.R.v=Av; clear A Av;
        end
    end
end
%=========================================================================%
if strcmp(SL.design.calc{indx},'Identity1')
    % Initiliaze f matrizes
    for jj=1:length(SL.design.anova{indx}.names)
        SL.design.anova{indx}.f{jj}=zeros(sum(Box));
    end
end
%=========================================================================%        
% For each specified interaction
for ii=1:size(SL.design.interactions{indx},1)   
    % Setup box interactions
    if SL.design.interactions{indx}(ii,1)==1, r=1:sum(Box(1:SL.design.interactions{indx}(ii,1)));
    else r=sum(Box(1:SL.design.interactions{indx}(ii,1)-1))+1:sum(Box(1:SL.design.interactions{indx}(ii,1)));
    end
    if SL.design.interactions{indx}(ii,2)==1, c=1:sum(Box(1:SL.design.interactions{indx}(ii,2)));
    else c=sum(Box(1:SL.design.interactions{indx}(ii,2)-1))+1:sum(Box(1:SL.design.interactions{indx}(ii,2)));
    end  
    
    % Box should be same same size if Identity model is being run
    if strcmp(SL.design.calc{indx},'Identity1')
        if length(r)~=length(c)
            % SL.err=1;
            sdisp('Warning Identity1 is misset - continuing...',1);
            % Hemi breaks this, no matter though... the code handles this
            % exception just fine.
        end
    end
    %=====================================================================%
    % Arrange anova matrices
    %=====================================================================%
    if SL.design.custom(indx)==0
        SL.design.matrix{indx}(r,c)=1;
        % Setup the ANOVA mask for the design matrix - this is fine to do
        % internally as we'll only filling r/c rectangles
        if strcmp(SL.design.calc{indx},'Identity1') 
            if ismember(SL.design.interactions{indx}(ii,:),SL.design.anova{indx}.diag,'rows');
                SL.design.anova{indx}.f{1}(r,c)=eye(length(r),length(c));
                SL.design.anova{indx}.f{2}(r,c)=double(~eye(length(r),length(c)));
            else
                SL.design.anova{indx}.f{2}(r,c)=double(ones(length(r),length(c)));
            end
        end
    %=====================================================================%
    % Arrange custom design matrices
    %=====================================================================%
    elseif SL.design.custom(indx)==1
        % Loop through all positions
%         keyboard
        if isequal(stim_ID_num,SL.design.ID_idx)
            SL.design.matrix{indx} = R;
        else
            for x=r
                for y=c
                    SL.design.matrix{indx}(x,y)=...
                        R(stim_ID_num==SL.design.ID_idx(x),...
                        stim_ID_num==SL.design.ID_idx(y));
                end
            end
        end
    end
end
% Parameter is really only pertinent to continous models...
SL.design.power{indx}=sum(SL.design.matrix{indx}(~isnan(SL.design.matrix{indx}))/2);

%=========================================================================%
% Remove the diagnol
%=========================================================================%
SL.design.matrix{indx}=tril(SL.design.matrix{indx},-1);
nan_idx=tril(ones(length(SL.design.matrix{indx})),-1)==0;
SL.design.matrix{indx}(nan_idx)=nan;
if SL.run.include==1
    % Within run comparisons set to NaN
    SL.design.matrix{indx}(SL.run.matrix==1)=NaN;
elseif SL.run.include==2
    % Between run comparisons set to NaN
    SL.design.matrix{indx}(SL.run.matrix==0)=NaN;
end

if strcmp(SL.design.calc{indx},'Identity1');
    if isfield(SL.design.anova{indx},'f')
        for jj=1:length(SL.design.anova{indx}.f)
            SL.design.anova{indx}.f{jj}(nan_idx)=0; % Anova map, so ==0 
            if SL.run.include==1, SL.design.anova{indx}.f{jj}(SL.run.matrix==1)=0; end
            
            if std(std(SL.design.anova{indx}.f{jj}))==0
                SL.design.anova{indx}.f{jj}=[];
                display('Warning: Anova Matrices have no variability, removing');
            end
        end
        SL.design.anova{indx}.f=SL.design.anova{indx}.f(~cellfun('isempty',SL.design.anova{indx}.f));
    end
end
%=========================================================================%
% Linearize Identity1 output
%=========================================================================%
ii=indx;
if strcmp(SL.design.calc{indx},'Identity1');
    SL.design.anova{ii}.f{1}=SL.design.anova{ii}.f{1}.*SL.design.matrix{ii};
    SL.design.anova{ii}.f{2}=SL.design.anova{ii}.f{2}.*SL.design.matrix{ii};
    SL.design.anova{ii}.f{1}(isnan(SL.design.anova{ii}.f{1}))=0;
    SL.design.anova{ii}.f{2}(isnan(SL.design.anova{ii}.f{2}))=0;
    SL.design.anova{ii}.fIon=[];
    SL.design.anova{ii}.fIof=[];
    if strcmp(SL.design.anova{ii}.type,'Identity1')
        MposC=find(sum(SL.design.anova{ii}.f{1},1)==1); % Positions with ones (columns)
        MposR=find(sum(SL.design.anova{ii}.f{1},2)==1); % Positions with ones (rows)
        for jj=1:length(MposR)
            A1=zeros(size(SL.design.anova{ii}.f{1})); % I matrix
            A1(MposR(jj),MposC(jj))=1;                % Hit Index
            A2=zeros(size(SL.design.anova{ii}.f{2})); % O matrix
            switch SL.design.anova{ii}.row
                case 0, A2(MposR,MposC(jj))=1;
                case 1, A2(MposR(jj),MposC)=1; 
                case 2, A2(MposR,MposC(jj))=1;
                        A2(MposR(jj),MposC)=1;
            end
            A2=A2-A1; % Remove identity
            SL.design.anova{ii}.fIon{jj}=reshape(A1,1,[]);
            SL.design.anova{ii}.fIof{jj}=reshape(A2,1,[]);
            clear A1 A2;
        end
    end
end

if strcmp(SL.design.calc{ii},'Anova1')
    SL.design.f{ii}=reshape(SL.design.matrix{ii},1,[]);
end
