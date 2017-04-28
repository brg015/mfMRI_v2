function [out_index,out_array]=RSA_model_univariate(tmp,ii)
global SL;

Events=nansum(SL.design.matrix{ii})>0;
m=mean(mean(tmp(:,Events)));
s=std(std(tmp(:,Events)));
N=sum(Events);
T=m/(s/sqrt(N));
out_array={T m s N};
out_index={'Key' [SL.design.save_str{ii} '_mean'] ...
[SL.design.save_str{ii} '_std'] ...
[SL.design.save_str{ii} '_N']};