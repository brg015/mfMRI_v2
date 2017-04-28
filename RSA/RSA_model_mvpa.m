function RSA_model_mvpa()

%=============================================================%
% MVPA Calculation (binary only)
%=============================================================%  
% 1) Collect data & Label (also setup splits while labeling)
Data = SL.files(SL.LOC(curLOC).box,:)';
Data=Data(:,~isnan(sum(Data)));              % Rm NaN
Data=Data-mean(Data,2)*ones(1,size(Data,2)); % Subtract mean
%         if size(Data,2)<40, continue; end % arbit limit on voxel number

Label = zeros(1,sum(SL.design.Box));
Set = zeros(1,sum(SL.design.Box));
for jj=1:length(SL.design.Box)
    [I,~]=cell_block(jj,jj,SL.design.Box);
    Label(I)=SL.design.mvpa{ii}.cond(jj);
    Set(I)=Randi(SL.design.mvpa{ii}.iter,[1,length(I)]);
end
% 2) SVM testing

% Random data
Data=Randi(100,size(Data));
Label(1:2:end)=1;
Label(2:2:end)=2;
% 'boxconstraint' lower
% 'kktviolationlevel'positive
% 'tolkkt'  
for jj=1:SL.design.mvpa{ii}.iter
    Itrain=Set~=jj;
    Itest=Set==jj;
    % Train
    % Training: [observations X features]
    % Group:    [label of observations]
    SVMStruct=svmtrain(Data(Itrain,:),Label(Itrain),'kktviolationlevel',0.5);
    % Test
    % SVMStruct: from above
    % Sample:    [observations X features]
    % Group:     [predicted value]
    Guess=svmclassify(SVMStruct,Data(Itest,:));
    % Evaluate feedback
    acc(jj)=mean(Guess==Label(Itest)');
    clear SVMStruct Guess
end
v=mean(acc); T=(v-0.5)/(std(acc)/sqrt(length(acc)));
out{strcmp([SL.design.save_str{ii} '_v'],out_name)}.mat(SL.LOC(curLOC).voi,1)=v;
out{strcmp([SL.design.save_str{ii} '_T'],out_name)}.mat(SL.LOC(curLOC).voi,1)=T;