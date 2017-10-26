function [ACC,D1,D2,Dclass1,Dclass2]=RSA_model_distance_v2(tmp,ii)
% Give me...
% distance measures
% templates...
global SL;
%-------------------------------------------------------------------------%
% Model distance...
%-------------------------------------------------------------------------%
Dclass1=tmp(:,SL.design.classID{ii}==1);
Dclass2=tmp(:,SL.design.classID{ii}==0);
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

Fclass=(D1<D2)'; % if D1(class1) is closer then Flcass is true i.e. 1
ACC=mean(Fclass'==SL.design.classID{ii});
