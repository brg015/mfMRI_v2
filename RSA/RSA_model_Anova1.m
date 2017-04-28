function [out_index,out_array]=RSA_model_Anova1(R,ii)
global SL;
y_data=[]; % input data
m_data=[]; % Condition effect
ms=[];

for jj=1:length(SL.design.anova{ii}.f)
    y_data=[y_data R(SL.design.anova{ii}.f{jj}==1)'];
    L=sum(sum(SL.design.anova{ii}.f{jj}));
    m_data=[m_data ones(1,L)*SL.design.anova{ii}.cond(jj)];
    ms((jj-1)*3+1)=mean(R(SL.design.anova{ii}.f{jj}==1)); % Value
    ms((jj-1)*3+2)=std(R(SL.design.anova{ii}.f{jj}==1));  % Std
    ms((jj-1)*3+3)=sum(sum(SL.design.anova{ii}.f{jj}==1));  % Nsamples
end

% Running anova and saving the output
try
    [~,T,STATS,~]=anovan(y_data',{m_data'},...
        'varnames',{'Cond'},'model','interaction','display','off');
catch err
    display('I broke again');
    return;
    % keyboard;
end
out_array{1}=T{2,6};
out_index={'Key'};

for kk=1:length(SL.design.anova{ii}.names)
    out_array=[out_array ms((kk-1)*3+1)];
    out_index=[out_index [SL.design.save_str{ii} '_' SL.design.anova{ii}.names{kk} '_v']];

    out_array=[out_array ms((kk-1)*3+2)];
    out_index=[out_index [SL.design.save_str{ii} '_' SL.design.anova{ii}.names{kk} '_sd']];

    out_array=[out_array ms((kk-1)*3+3)];
    out_index=[out_index [SL.design.save_str{ii} '_' SL.design.anova{ii}.names{kk} '_N']];

    M=ms((kk-1)*3+1); SD=ms((kk-1)*3+2); N=ms((kk-1)*3+3);
    out_array=[out_array M/(SD/sqrt(N))];   
    out_index=[out_index [SL.design.save_str{ii} '_' SL.design.anova{ii}.names{kk} '_T']];
end

if isfield(SL.design.anova{ii},'coeffs'),
    for jj=1:length(SL.design.anova{ii}.coeffs.idx)
        out_index=[out_index SL.design.anova{ii}.coeffs.names{jj}];
        out_array=[out_array STATS.coeffs(SL.design.anova{ii}.coeffs.idx(jj))];
    end
end  
