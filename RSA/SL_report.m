function SL_report(out_name,subjN)
global SL;
% Setup the QA directory

QA_dir=fullfile(SL.dir.outpath,SL.dir.QA,filesep);
if ~exist(QA_dir,'dir'), mkdir(QA_dir); end

QA_subj=fullfile(QA_dir,SL.dir.subjects{subjN});
if ~exist(QA_subj,'dir'), mkdir(QA_subj); end

% File specific information 
sdisp('I/O',1);
display(['Reading from: ' SL.dir.stpath]);
display(['Saving to   : ' SL.dir.outpath]);
switch SL.dir.overwrite
    case 0, display('Overwrite   : OFF');
    case 1, display('Overwrite   : ON');
end
switch SL.region.use_mask
    case 0, display('Masking     : OFF');
    case 1, display('Masking     : ON');
end
%=========================================================================%
% Save some general output
%=========================================================================%
write_vbeta(subjN,QA_subj);

if SL.run.include==1
    RSA_pdm(SL.run.matrix,SL.design.Box,'Run Matrix',...
        fullfile(QA_subj,['run_matrix']));
end
    
for ii=1:length(SL.design.save_str)
    % Each file gets own model folder
    if ~exist(fullfile(QA_subj,SL.design.save_str{ii}),'dir')
       mkdir(fullfile(QA_subj,SL.design.save_str{ii}));
    end
    sdisp(['Model ' n2sp(ii,2) ': Type ' SL.design.calc{ii} ' : ' SL.design.save_str{ii}],1);
    model_type=SL.design.calc{ii};
    if strcmp(model_type, 'Anova1'),    model_type='Anova'; end
    if strcmp(model_type, 'Identity1'), model_type='Anova'; end
    if strcmp(model_type, 'Spear'),     model_type='cont'; end
    if strcmp(model_type, 'Euclid'),    model_type='cont'; end
    if strcmp(model_type, 'Kendall'),   model_type='cont'; end
    RSA_pdm(SL.design.matrix{ii},SL.design.Box,'Design',...
        fullfile(QA_subj,SL.design.save_str{ii},'Design'));
    
    switch SL.design.calc{ii}
        case 'Anova1'
            RSA_pdm(SL.design.matrix{ii},SL.design.Box,'AnovaModel',...
                fullfile(QA_subj,SL.design.save_str{ii},'AnovaModel'));
        case 'Identity1'
            RSA_pdm(SL.design.anova{ii}.f{1},SL.design.Box,SL.design.Identity1{ii}.names{1},...
                fullfile(QA_subj,SL.design.save_str{ii},SL.design.Identity1{ii}.names{1}));
            RSA_pdm(SL.design.anova{ii}.f{2},SL.design.Box,SL.design.Identity1{ii}.names{2},...
                fullfile(QA_subj,SL.design.save_str{ii},SL.design.Identity1{ii}.names{2}));
        otherwise
            RSA_pdm(SL.design.matrix{ii},SL.design.Box,'ContModel',...
                fullfile(QA_subj,SL.design.save_str{ii},'ContModel'));
    end
end
     
if SL.region.noSL==0
    display('--Save Files S(1)--');
    display('DNE: Does not exist');
    display('FND: Found');
    display('--Save Files S(1)--');
    for ii=1:length(out_name)
        f=fullfile(SL.dir.outpath,SL.dir.subjects{1},[out_name{ii} '.img']);
        switch exist(f,'file')
            case 0, display(['DNE: ' f]);
            case 2, display(['FND: ' f]);
        end
    end
end               
                        



