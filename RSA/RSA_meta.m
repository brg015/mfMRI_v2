function model=RSA_meta(cond_code,save_str,model)
global SL; 
%=========================================================================%
% Setup models for all types of interactions wanted...
% Inputs:
%   cond_code: 1 = encoding low memory
%   cond_code: 2 = encoding high memory
%   cond_code: 3 = retrieval low memory
%   cond_code: 4 = retrieval high memory
%   save_str = prefix of model name;
%=========================================================================%
% Determine model cells from cond_code
%=========================================================================%
LL=setdiff(setdiff(combn(1:length(cond_code),2),...
    combn(find(cond_code>2),2),'rows'),combn(find(cond_code<3),2),'rows');

LL_ME=[];
R=find(cond_code>2);
E=find(cond_code<3);

for jj=1:length(R)
    LL_ME=[LL_ME; E(jj) R(jj); R(jj) E(jj)];
end

%=========================================================================%  
S{1}=combn(find(cond_code==2),2);
S{2}=combn(find(cond_code==1),2);
S{3}=intersect(LL,combn(find(or(cond_code==2,cond_code==4)),2),'rows');
S{4}=intersect(LL,combn(find(or(cond_code==1,cond_code==3)),2),'rows');
S{5}=combn(find(cond_code==4),2);
S{6}=combn(find(cond_code==3),2);
S{7}=combn(find(or(cond_code==1,cond_code==2)),2);
S{8}=LL;
S{9}=combn(find(or(cond_code==3,cond_code==4)),2);

head={'EncHit' 'EncMiss' 'ERHit' 'ERMiss' 'RetHit' 'RetMiss' 'Enc' 'ER' 'Ret'};
%=========================================================================%
% AtoB (x4)
%=========================================================================%
% SL.design.interactions{model}=[S{1}; S{5}];
% SL.design.save_str{model}=[save_str '_AtoB_ERHit'];
% SL.design.custom(model)=0;
% SL.design.calc{model}='AtoB';
% SL.design.E{model}=E;
% SL.design.R{model}=R;
% model=model+1;
% 
% SL.design.interactions{model}=[S{2}; S{6}];
% SL.design.save_str{model}=[save_str '_AtoB_ERMiss'];
% SL.design.custom(model)=0;
% SL.design.calc{model}='AtoB';
% SL.design.E{model}=E;
% SL.design.R{model}=R;
% model=model+1;
%=========================================================================%
% Setup Univaraite Models (x4)
%=========================================================================%
% for ii=[1 2 5 6]
%     SL.design.interactions{model}=[S{ii}];
%     SL.design.save_str{model}=[save_str '_Univariate_' head{ii}];
%     SL.design.custom(model)=0;
%     SL.design.calc{model}='Univariate';
%     model=model+1;
% end
%=========================================================================%
% Mean Series (x4)
%=========================================================================%
for ii=[1 2 3 4 5 6 7 8 9]
    SL.design.interactions{model}=[S{ii}];
    SL.design.save_str{model}=[save_str '_Anova1_' head{ii}];
    SL.design.custom(model)=0;
    SL.design.calc{model}='Anova1';
    
    SL.design.anova{model}.cond=ones(size(SL.design.interactions{model},1),1);
    SL.design.anova{model}.names{1}='ME';
    SL.design.anova{model}.type='Anova1';
    model=model+1;
end
%=========================================================================%
% Setup Identity Models (x6)
%=========================================================================%
% E,R,ER X Hit/Miss
for ii=[8 4 3]
    % Use LL_ME because we want to compare across memory conditions
    SL.design.interactions{model}=S{ii};
    SL.design.save_str{model}=[save_str '_Identity1_' head{ii}];
    SL.design.custom(model)=0;
    SL.design.calc{model}='Identity1';
    
    SL.design.anova{model}.type='Identity1';
    SL.design.anova{model}.names={'On' 'Off'};
    SL.design.anova{model}.diag=LL_ME;
    if ii==8
        SL.design.interactions{model}=LL_ME;
        SL.design.anova{model}.regress.on=1;
        SL.design.anova{model}.regress.E.uni=1;
        SL.design.anova{model}.regress.E.pos=E;
        SL.design.anova{model}.regress.R.uni=1;
        SL.design.anova{model}.regress.R.pos=R;
        SL.design.anova{model}.regress.ER.uni=1;
        SL.design.anova{model}.regress.Cseries=1;
        SL.design.anova{model}.regress.Zseries=1;
        SL.design.anova{model}.regress.measures_prd=[1];
        SL.design.anova{model}.regress.map=[1 -1; 2 1; 3 -1; 4 1];
    else
        SL.design.anova{model}.regress.on=0;
    end
    model=model+1;
end




