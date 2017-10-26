function [out_index,out_array]=RSA_model_MVPA_SVM(tmp,ii)
global SL;
%-------------------------------------------------------------------------%
% Bare bones MVPA (slow ~.82 seconds a loop on random test)
%-------------------------------------------------------------------------%
SVMModel=fitcsvm(tmp',SL.design.classID{ii},'KernelFunction','linear',...
    'Standardize',true,'Leaveout','on');
out_array={1-SVMModel.kfoldLoss-.5}; % Loss is misclass
out_index={[SL.design.save_str{ii} '_key']};




    
    
    