% Incredibly simple contrast manager via GT_contrast
function RSA_contrast(SL)

spm_jobman('initcfg');
for c=1:length(SL.con)
    for ii=1:length(SL.con(c).subjects)
        display([' Subject: ' SL.con(c).subjects{ii} ' - ' SL.con(c).save_name]);
        if exist(fullfile(SL.con(c).dir,SL.con(c).subjects{ii},SL.con(c).save_name),'file')
%             continue;
        end
        if isempty(SL.con(c).func), SL.con(c).func=1; end
        try
        switch SL.con(c).func
            case 'PT'     
                I1m=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},[SL.con(c).files{1} '_v.img']); % Image 1 mean
                I2m=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},[SL.con(c).files{2} '_v.img']); % Image 2 mean
                I1s=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},[SL.con(c).files{1} '_sd.img']); % Image 1 std
                I2s=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},[SL.con(c).files{2} '_sd.img']); % Image 2 std
                load(fullfile(SL.con(c).dir,SL.con(c).subjects{ii},[SL.con(c).files{1} '_v.mat']));
                N1=N; clear N; % Samples in 1 mean
                load(fullfile(SL.con(c).dir,SL.con(c).subjects{ii},[SL.con(c).files{2} '_v.mat']));
                N2=N;  % Samples in 2 mean

                ms=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},'temp1.img');   % Name of mean image
                s1=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},'temp2.img');    % Name of std1
                s2=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},'temp3.img');    % Name of std2
                s=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},'temp4.img');     % Name of summed std
                ps=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},'temp5.img');    % Name of pooled std
                T=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},SL.con(c).save_name);     % Name of T image

                GT_contrast({I1m,I2m},[1 -1],ms)
                GT_contrast({I1s},1,s1,'divide',N1)
                GT_contrast({I2s},1,s2,'divide',N2)
                GT_contrast({s1,s2},[1 1],s)
                GT_contrast({s},1,ps,'function','sqrt')
                GT_contrast({ms,ps},1,T,'divide_image')
     
            otherwise
                for jj=1:length(SL.con(c).files)
                    f{jj}=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},SL.con(c).files{jj});
                end
                save_name=fullfile(SL.con(c).dir,SL.con(c).subjects{ii},SL.con(c).save_name);
                GT_contrast(f,SL.con(c).w,save_name);     
%                        GT_contrast(f,SL.con(c).w,save_name,'function',SL.con(c).func);
                
        end
        catch err
            display('  Couldn''t Process - Skipped');
        end
    end
end
