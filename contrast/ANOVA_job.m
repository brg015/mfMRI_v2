
function matlabbatch=ANOVA_job(save_dir,level,design)
if ~exist(save_dir,'dir'), mkdir(save_dir); end

matlabbatch{1}.spm.stats.factorial_design.dir = {save_dir};
for ii=1:length(level)
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(ii).name = level(ii).name;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(ii).levels = level(ii).factor;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(ii).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(ii).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(ii).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(ii).ancova = 0;
end
%-------------------------------------------------------------------------%
% Design levels
%-------------------------------------------------------------------------%
for ii=1:prod([level(:).factor])
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(ii).scans = design.inimgs(:,ii);
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(ii).levels=[design.V(ii,:)];
end
   
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(save_dir,'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

