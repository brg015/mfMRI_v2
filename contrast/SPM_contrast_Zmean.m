clc;
root_dir='F:\Data2\CODE\';
code_dir{1}=fullfile(root_dir,'mfMRI_v2-master');
code_dir{2}=fullfile(root_dir,'xjview');
code_dir{3}=fullfile(root_dir,'function_files-master');
addpath('D:\Data\fmri\spm8');
for ii=1:length(code_dir)
    addpath(genpath(code_dir{ii}));
end
% spm fmri
spm_jobman('initcfg'); 
%-------------------------------------------------------------------------%
save_dir=fullfile('F:\Data2\Geib\ER_Match\Hyper_Results\Mem_Zintaxn\');

wrk_dir=fullfile('F:\Data2\Geib\ER_Match\Hyper_Results\');
subjects={'13100' '13205' '13220' '13448' ...
    '13510' '13523' '13617' '13743' '13779' '13793' '13807' '13845' ...
    '13491' '13549' '13552' '13562' '13655' '13658' '13683' ...
    '13693' '13720'};

imgnames = {'ERS_v1_Zmean' 'ERS_v2_Zmean' 'ERS_v3_Zmean' 'ERS_v4_Zmean'};

numsubjs = length(subjects);
numcells = length(imgnames);

inimgs = cell(length(subjects),length(imgnames));
for i = 1:length(imgnames)
    for s = 1:length(subjects)
        inimgs(s,i) = {strcat(wrk_dir,subjects{s},filesep,[imgnames{i} '.img'])};
    end
end
%-------------------------------------------------------------------------%
% level(1).name='ID/OffID';
% level(1).factor=2;
level(1).name='Memory';
level(1).factor=4;

design.inimgs=inimgs;
design.numsubjs=numsubjs;
design.numcells=numcells;
design.V=[1; 2; 3; 4];

matlabbatch=ANOVA_job(save_dir,level,design);
spm_jobman('run',matlabbatch);
