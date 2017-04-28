function GT_compare_measures(I1,I2)
global ntwk SL;
% I1 => ntwk analysis
% I2 => contrast analysis
save_dir_neg=fullfile(SL.dir.save,ntwk.contrast.neg{I2});
save_dir_pos=fullfile(SL.dir.save,ntwk.contrast.pos{I2});
save_dir_nam=fullfile(SL.dir.save,ntwk.contrast.nam{I2});

F1=load(fullfile(save_dir_pos,[ntwk.analysis{I1} '.mat']));
F2=load(fullfile(save_dir_neg,[ntwk.analysis{I1} '.mat']));

if strcmp(ntwk.analysis{I1},'Mod1') || strcmp(ntwk.analysis{I1},'wMod1')
   save1=fullfile(save_dir_nam,[ntwk.analysis{I1} '.png']);
   [~,p,~,stats]=ttest(F1.MAT-F2.MAT);
   display('Mod Results');
   display(['  p-value: ' num2str(p)]);
   display(['  T-value: ' num2str(stats.tstat)]);
   h=notBoxPlot([F1.MAT;F2.MAT]');
   title(['p-value: ' num2str(p)]); ylabel('Q-value');
   set(gca,'xticklabel',{'Positive' 'Negative'});
   set(gcf,'position',[0 0 1280 1024]); export_fig(save1);
   close all;
   return;
end

% Load in competing representations
for ii=1:size(F1.MAT,2)
    [p,~,stats]=ranksum(F1.MAT(:,ii),F2.MAT(:,ii));
    A1(ii)=p;
    if isfield(stats,'zval'), A2(ii)=stats.zval; end
end

% Ranksum Test
MAT=A1;
save(fullfile(save_dir_nam,[ntwk.analysis{I1} '.mat']),'MAT');
save_asA=fullfile(save_dir_nam,[ntwk.analysis{I1} '_ranksum_p.img']);
[~]=ROI_constructor(SL.region.mask,MAT,'','pure',save_asA,3);  
if exist('A2','var');
    MAT=A2;
    save(fullfile(save_dir_nam,[ntwk.analysis{I1} '.mat']),'MAT');
    save_asA=fullfile(save_dir_nam,[ntwk.analysis{I1} '_zval.img']);
    [~]=ROI_constructor(SL.region.mask,MAT,'','pure',save_asA,3); 
end
% Ttest area
MAT=F1.MAT-F2.MAT;
[~,p,~,stats]=ttest(MAT);
save(fullfile(save_dir_nam,[ntwk.analysis{I1} '.mat']),'MAT');
save_asA=fullfile(save_dir_nam,[ntwk.analysis{I1} '_Ttest_p.img']);
[~]=ROI_constructor(SL.region.mask,p,'','pure',save_asA,3);  
save(fullfile(save_dir_nam,[ntwk.analysis{I1} '.mat']),'MAT');
save_asA=fullfile(save_dir_nam,[ntwk.analysis{I1} '_tval.img']);
[~]=ROI_constructor(SL.region.mask,stats.tstat,'','pure',save_asA,3);  

% More fancier T-Test stuff
G1=mean(F1.MAT); G2=std(F1.MAT)./sqrt(size(F1.MAT,1));
H1=mean(F2.MAT); H2=std(F2.MAT)./sqrt(size(F2.MAT,1));
data{1}.header='POSvalue';   data{1}.col=G1;
data{2}.header='POSsem';     data{2}.col=G2;
data{3}.header='NEGvalue';   data{3}.col=H1;
data{4}.header='NEGsem';     data{4}.col=H2;
data{5}.header='Tstat';      data{5}.col=stats.tstat;
data{6}.header='p-value';    data{6}.col=p;
write_struct(data,fullfile(save_dir_nam,[ntwk.analysis{I1} '_tstats_full.csv']));
clear data;
% 
clear p;

% Permutation test.
for ii=1:size(F1.MAT,2),
    p(ii)=run_permutation(F1.MAT(:,ii),F2.MAT(:,ii),1,save_dir_nam,0,[ntwk.analysis{I2} '_perm'],300000);
end

MAT=p;
% csvwrite(fullfile(save_dir_nam,[ntwk.analysis{I1} '_permutation.csv']),MAT)
save(fullfile(save_dir_nam,[ntwk.analysis{I1} '_permutation.mat']),'MAT');
save_asA=fullfile(save_dir_nam,[ntwk.analysis{I1} '_permutation.img']);
[~]=ROI_constructor(SL.region.mask,MAT,'','pure',save_asA,3);  




