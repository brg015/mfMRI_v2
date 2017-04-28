function [out_index,out_array]=RSA_model_Identity1_v2(R,ii)
global SL;
% ii      => model index
% tmp     => correlation input matrix
% R       => correlation matrix
% outname =>
% 
%-------------------------------------------------------------------------%
% Match vs. Mismatch
%-------------------------------------------------------------------------%
% 1) Grab wanted values
y_data=[R(SL.design.anova{ii}.f{1}==1)' R(SL.design.anova{ii}.f{2}==1)'];
m_data=[ones(1,length(R(SL.design.anova{ii}.f{1}==1))) zeros(1,length(R(SL.design.anova{ii}.f{2}==1)))];
% Initialize
out_array={NaN NaN NaN NaN(1,sum(m_data)) NaN(1,sum(m_data))};
% Even: MM  (mismatch)  => 1
% Odd : M   (match)     => 0
MMmean=mean(y_data(m_data==0));
Mmean=mean(y_data(m_data==1));
Iseries=y_data(m_data==1); %added, same as Iseries(jj)=v1(1) below

out_array{1}=MMmean; out_array{2}=Mmean;
out_index={[SL.design.save_str{ii} '_MMmean'],[SL.design.save_str{ii} '_Mmean']};
%-------------------------------------------------------------------------%
% Z series
%-------------------------------------------------------------------------%
% Note that the row and column distinction is handled over in
% RSA_output_maps.m - the indicies are set according to this distinction
% and we no longer care by this point in time.

% computation of z-score (below) commented out

Rl=reshape(R,1,[]);
for jj=1:length(SL.design.anova{ii}.fIon)
    % Pull values from Ion and Iof
    v1=[Rl(SL.design.anova{ii}.fIon{jj}==1) Rl(SL.design.anova{ii}.fIof{jj}==1)]; 
    % z-score
    v2=(v1-mean(v1))./std(v1);
    % save values
    Zseries(jj)=v2(1);
    Iseries(jj)=v1(1); 
end
out_array{3}=mean(Zseries);

if SL.design.Identity1{ii}.FourD==1
    out_array{4}=Zseries;
    out_array{5}=Iseries;
end
   

out_index=[out_index ...
    [SL.design.save_str{ii} '_Zmean'] ...
    [SL.design.save_str{ii} '_Zseries'] ...
    [SL.design.save_str{ii} '_Iseries']];

%=========================================================================%
% Regression calculations (added 12/15/14)
%=========================================================================%
% if SL.design.anova{ii}.regress.on==1
% n={};
% % Collect univariate measures
% A=[]; c=1;
% if SL.design.anova{ii}.regress.E.uni==1
%     A(:,c)=mean(tmp(:,SL.design.anova{ii}.regress.E.v)); c=c+1;
%     n=[n 'E'];
% end
% if SL.design.anova{ii}.regress.R.uni==1
%     A(:,c)=mean(tmp(:,SL.design.anova{ii}.regress.R.v)); 
%     n=[n 'R']; c=c+1;
% end
% if SL.design.anova{ii}.regress.ER.uni==1
%     A(:,c)=A(:,1).*A(:,2);
%     n=[n 'ER'];
% end
% if ~isempty(A), A=zscore(A); end
% 
% % Collect Zseries measures
% B=[]; c=1;
% if SL.design.anova{ii}.regress.Cseries==1
%     B(:,c)=zscore(Cseries); c=c+1;
%     n=[n 'Cseries'];
% end
% if SL.design.anova{ii}.regress.Zseries==1
%     B(:,c)=zscore(Zseries);
%     n=[n 'Zseries'];
% end
% 
% C=[]; c=1;
% if sum(~logical(SL.design.anova{ii}.regress.measures_prd))>0
%     C=SL.regress.val(Mpos,~logical(SL.design.anova{ii}.regress.measures_prd))==1;
%     n=[n SL.regress.name(~logical(SL.design.anova{ii}.regress.measures_prd))];
% end
%   
% % Behav         => Y
% % This needs a better fix later on, but being lazy for the moment.
% try
%     Y=SL.regress.val(length(v1),logical(SL.design.anova{ii}.regress.measures_prd));
% catch err
%     Y=SL.regress.val(logical(SL.design.anova{ii}.regress.measures_prd),length(v1));
% end
% 
% Y2=Y;
% if ~isempty(SL.design.anova{ii}.regress.map)
%     for jj=1:size(SL.design.anova{ii}.regress.map,1)
%         % Prevent indexing from ruining
%         Y2(Y==SL.design.anova{ii}.regress.map(jj,1))=SL.design.anova{ii}.regress.map(jj,2);
%     end
% end
% X=[A,B,C]; 
% mdl=fitglm(X,Y2,'Distribution',SL.regress.type);
% BE=mdl.Coefficients.Estimate;
% if ~isempty(outname)
%     s=find(strcmp([SL.design.save_str{ii} '_Iseries'],outname));
%     for ii=2:length(BE), 
%         out_array=[out_array BE(ii)]; 
%         out_index=[out_index [s+ii-1]];
%     end
% else
%     for ii=2:length(BE), 
%         out_array=[out_array BE(ii)]; 
%         out_index=[out_index n{ii-1}];
%     end 
% end

%=========================================================================%    
end
