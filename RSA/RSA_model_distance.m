function [out_index,out_array]=RSA_model_distance(tmp,ii)
global SL;
%-------------------------------------------------------------------------%
% Model distance...
%-------------------------------------------------------------------------%
if strcmp(SL.design.calc{ii},'mahalanobis')
    for jj=1:length(SL.design.classID{ii})
       Xtest=tmp'; Xtest(jj,:)=[]; % Clear test sample
       Ctest=SL.design.classID{ii};  Ctest(jj)=[];   % Clear test case

       X1=Xtest(Ctest==1,:);  % Class1
       X2=Xtest(Ctest==0,:);  % Class0
       clear Xtest; Xsamp=tmp(:,jj)';       

       D1(jj)=pdist2(Xsamp,mean(X1),SL.design.distance{ii},nancov(X1));
       D2(jj)=pdist2(Xsamp,mean(X2),SL.design.distance{ii},nancov(X2)); 
       % Number of stdevs from a high-dimensional space
       clear X1 X2 Xsamp;
    end   
    % Use chi-squre distribution to determine liklihood
    % p1=chi2pdf(D1.^2,size(tmp,1)); p1m=mean(p1(SL.design.classID{ii}==1));
    % p2=chi2pdf(D2.^2,size(tmp,1)); p2m=mean(p2(SL.design.classID{ii}==0));
else
    for jj=1:length(SL.design.classID{ii})
       Xtest=tmp'; Xtest(jj,:)=[]; % Clear test sample
       Ctest=SL.design.classID{ii};  Ctest(jj)=[];   % Clear test case

       X1=Xtest(Ctest==1,:);  % Class1
       X2=Xtest(Ctest==0,:);  % Class0
       clear Xtest; Xsamp=tmp(:,jj)';       

       D1(jj)=mean(pdist2(X1,Xsamp,SL.design.distance{ii}));
       D2(jj)=mean(pdist2(X2,Xsamp,SL.design.distance{ii}));
       clear X1 X2 Xsamp;
    end
end

Fclass=(D1<D2)'; % if D1(class1) is closer then Flcass is true i.e. 1


out_array={mean(Fclass==SL.design.classID{ii})-.5}; 
out_index={[SL.design.save_str{ii} '_key']};