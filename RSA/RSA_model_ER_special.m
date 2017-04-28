function out=RSA_model_ER_special(R,vox,tmp,out)
global SL;
%=========================================================================%
% Optimized for ROIs atm
%=========================================================================%
% Find similarity structure via tmp
ztmp=zscore(tmp'); ztmp=ztmp';
dtmp=detrend(tmp,'constant');

dR=corr(dtmp); dRl=reshape(dR,1,[]);
zR=corr(ztmp); zRl=reshape(zR,1,[]);
for ii=1:96
    I=ii:96:2016;
    a1=mean(tmp(:,I)');
    a2=std(tmp(:,I)');
end

Rl=reshape(R,1,[]); 
% For each phase
c=1;
% +96  -> Ret
% +192 -> ER
for ii=1 % Just encoding for now, due to time...
    for jj=1:96
        vOn=Rl(SL.design.anova{1}.fl_on{192+c});  
        vOnP{ii,jj}=vOn;
        vOf=Rl(SL.design.anova{1}.fl_of{192+c});
        c=c+1;
        Zt=zscore([vOn,vOf]);
        u=mean(Zt(1:420)); s=std(Zt(1:420)); T=u/(s/sqrt(420));
        B1(jj)=T;
        B2(jj)=mean(vOn);
        B3(jj)=mean(vOf);
        clear Ion Iof vOn vOf Zt u s T;  
    end
    out{1+(ii-1)*7}.mat(vox,:)=mean(B2);
    out{2+(ii-1)*7}.mat(vox,:)=mean(B3);
    out{3+(ii-1)*7}.mat(vox,:)=mean(B1)/(std(B1)/sqrt(96));
    out{4+(ii-1)*7}.mat(vox,:)=B1;
    % clear B1 B2 B3;
%     vOn=Rl(SL.design.anova{1}.fat_on{ii});   
%     vOf=Rl(SL.design.anova{1}.fat_off{ii});   
%     out{5+(ii-1)*7}.mat(vox,:)=mean(vOn);
%     out{6+(ii-1)*7}.mat(vox,:)=mean(vOf);
%     [h,p,ci,stats]=ttest2(vOn,vOf,'Vartype','unequal');
%     out{7+(ii-1)*7}.mat(vox,:)=stats.tstat;
    % clear vOn vOf;
end

if isnan(mean(B2)), keyboard; end
v={mean(B2), mean(B3), mean(B1)/(std(B1)/sqrt(96)), B1};
