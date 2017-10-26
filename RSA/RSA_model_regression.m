function [out_index,out_array]=RSA_model_regression(R,ii)
global SL;

% Pull out the predictor from brain data
nan_idx=isnan(reshape(SL.design.matrix{ii}{1},1,[]));
y_data=reshape(R,[],1); y_data(nan_idx)=[];
m_data=[SL.design.MR{ii} ones(size(SL.design.MR{ii},1),1)]; % include int

% Initialize Output
% out_array={NaN(1,length(SL.design.matrix{ii}))}; % Is NaNs
out_index=SL.design.MRname{ii};                % Is Preset

% Run multiple regression
b=regress(y_data,m_data);
out_array=num2cell(b(1:end-1));   % don't include intercept