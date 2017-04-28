function [out_index,out_array]=RSA_model_Identity1_hemi(R,tmp,ii,outname)
global SL;


% 1) Grap wanted values
y_data=[]; m_data=[];
for jj=1:length(SL.design.anova{ii}.f)
    y_data=[y_data R(SL.design.anova{ii}.f{jj}==1)'];
    L=sum(sum(SL.design.anova{ii}.f{jj}));
    m_data=[m_data ones(1,L)*rem(jj,2)];
end
% Even: MM  (mismatch)  => 1
% Odd : M   (match)     => 0
MMmean=mean(y_data(m_data==0));
% MMstd=std(y_data(m_data==0));
% MMn=length(y_data(m_data==0));
Mmean=mean(y_data(m_data==1));
% Mstd=std(y_data(m_data==1));
% Mn=length(y_data(m_data==1));
% PairT=(Mmean-MMmean) / (sqrt(Mstd/Mn+MMstd/MMn));
out_array={NaN NaN};
out_array={MMmean Mmean};
out_index={[SL.design.save_str{ii} '_MMmean'] ...
    [SL.design.save_str{ii} '_Mmean']};
% 4D Timeseries outputs (added 12/10/14)
Mpos=find(sum(SL.design.anova{ii}.f{1})==1);  
LMpos=length(Mpos); 

Zseries=NaN;
Iseries=NaN;

% New edits (6/5/15) - this can surely be streamlined...
% switch SL.design.anova{ii}.row
%     case 0
%         for jj=1:LMpos
%             % Collect index values
%            Ii=logical(SL.design.anova{ii}.f{1}(:,Mpos(jj))); % Ident trial jj
%            Io=logical(SL.design.anova{ii}.f{2}(:,Mpos(jj))); % Off   trial jj
%            I=or(Ii,Io);                       % Either
%            v4d=R(I,Mpos(jj));                       % 4d vector
%            v4o=R(Io,Mpos(jj));                      % SET Vector
%            zv4d=zscore(v4d);                  % Zscore 4d
%            Ihit=find(Ii)==find(I);
%            Cseries(jj)=mean(v4o);    % column series
%            Zseries(jj)=zv4d(Ihit);   % Z series
%            Iseries(jj)=v4d(Ihit);    % Iseries
%         end
%     case 1

        for jj=1:LMpos
            % Collect index values along column finds...
           Ii=logical(SL.design.anova{ii}.f{1}(:,Mpos(jj)));     % Ident trial jj
           
           % Grab the relvent row
           Io=logical(SL.design.anova{ii}.f{2}(Ii,:));           % Off trial jj *rows*
           
           I=Io; I(Mpos(jj))=1;
           v4d=R(Ii,I);                       % 4d vector
%            v4o=R(Ii,Io);                      % SET Vector
           zv4d=zscore(v4d);                  % Zscore 4d
           % Find ident match equal to ident match or any?
           % Ihit=find(Ii)==find(or(Ii,Iotrue)); before 10/12
           Ihit=Mpos(jj)==find(I);
%            Cseries(jj)=mean(v4o);    % row series
           Zseries(jj)=zv4d(Ihit);   % Z series
           Iseries(jj)=v4d(Ihit);    % Iseries
        end  
% end
% Zzseries=zscore(Zseries)';

out_array=[out_array mean(Zseries) Iseries];
   
out_index=[out_index ...
    [SL.design.save_str{ii} '_Zmean'] ...
    [SL.design.save_str{ii} '_Iseries']];
%=========================================================================%
% Regression calculations (added 12/15/14)
%=========================================================================%
% if SL.design.anova{ii}.regress.on==1
% 
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
%     Y=SL.regress.val(Mpos,logical(SL.design.anova{ii}.regress.measures_prd));
% catch err
%     Y=SL.regress.val(logical(SL.design.anova{ii}.regress.measures_prd),Mpos);
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
% 
% %=========================================================================%    
% end
