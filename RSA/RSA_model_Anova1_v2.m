function [out_index,out_array]=RSA_model_Anova1_v2(R,ii)
global SL;

RL=reshape(R,1,[]);
v=mean(RL(SL.design.f{ii}==1));
% vsd=std(RL(SL.design.f{ii}==1));

out_array={v};
out_index={[SL.design.save_str{ii} '_key']};


