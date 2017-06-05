%-----------------------------------------------------------------------
% Job saved on 24-May-2016 18:48:01 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

numpairs = 19;
basepath = 'F:\Data2\SchemMemYA\Analysis\STnwa_art_mum_anatmask_tvals\mvpa_spm\STnwa_art_anatmaskbetas\';
folder1 = '\crossdec_congruence_SL\';
folder2 = '\crossdec_congruence_SL\';

subjects={'20278' '20287' '20364' '20414' '20502' '20545' '20284' '20346' ... 
                 '20399' '20491' '20512' '20288' '20465' ...
                 '20577' '20602' '20612' '20626' '20641' '20643'};  

imgname1 = 's8mm_res_sensitivity_minus_chance.nii';
imgname2 = 's8mm_res_specificity_minus_chance.nii';


matlabbatch{1}.spm.stats.factorial_design.dir = {'F:\Data2\SchemMemYA\Analysis\STnwa_art_mum_anatmask_tvals\mvpa_spm\STnwa_art_anatmaskbetas\secondlevel\congruence\pairedt19_crossdec_SensSpec_s8mm'};


for q = 1:numpairs
    vol1 = strcat(basepath,subjects{q},folder1,'\',imgname1);
    vol2 = strcat(basepath,subjects{q},folder2,'\',imgname2);
    matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(q).scans = {vol1; vol2};
end


matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'F:\Data2\SchemMemYA\Analysis\STnwa_art_mum_anatmask_tvals\mvpa_spm\STnwa_art_anatmaskbetas\mask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);
