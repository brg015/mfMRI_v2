function [out_index,out_array]=RSA_model_AtoB(R,ii)
global SL;
EE=R(SL.design.Ev{ii},SL.design.Ev{ii}); 
RR=R(SL.design.Rv{ii},SL.design.Rv{ii});
v1=EE(~isnan(SL.design.matrix{ii}(SL.design.Ev{ii},SL.design.Ev{ii})));
v2=RR(~isnan(SL.design.matrix{ii}(SL.design.Rv{ii},SL.design.Rv{ii})));
RHO=corr([v1,v2]);
out_array={RHO(1,2)};
out_index={'Key'};   