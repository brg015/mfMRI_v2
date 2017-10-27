function [out] = searchlight(subject,cursub)
global SL;
%=========================================================================%
%% RSA Style Searchligth
%=========================================================================%
% Consider
%   How to handle NaN?
%   => Currently are just removed (believe they are voxel specific)
%   => May want to limit how many can be removed and validness retained
%   Deviant values => [y,ntrimmed] = trimts(y,sd)?
%   => Currently ignored
%   How to handle statistics?
%   => Bootstrap like?
%   We have about 80,000 voxels to examine...
%
% Important Variables
%   SL.files  => [Voxels X files] Matrix
%   SL.LOC    => Voxel index and surrounding voxel indices
%=========================================================================%
%% Setup output matrices
%=========================================================================%
% Kriegskorte 2008
% Initialize out matrix [voxel X design]
% NOTES: These arrays are highly interwined with variables in
% RSA_bootstrap. The variable outname serves as a reference to the out
% matrices. RSA_bootstrap save values with the expectation of matching them
% to out_name. Thus if out_name is changed OR RSA_bootstrap out_index is
% changed, the two scripts can no longer talk.
%=========================================================================%
if (SL.region.noSL==0 && SL.skip_maps==0) %commented out for ROI
    [out,out_name]=RSA_output_maps_v2(subject);
    if SL.err==1, SL.err=0; return; end
else % End of making output maps (needed for searchlight)
    out={};
    out_name={};
end
%=========================================================================%
%% Check Design
%=========================================================================%
if SL.dir.QAcheck==1
    SL_report(out_name,cursub); % Generate folders, report plans
end

if isfield(SL.dir,'QAonly') 
    return
end

%=========================================================================%
%% Modify volumes if searchlight is turned off
%=========================================================================%
% MVPA based ROIs will not even make it here...
if SL.region.noSL==1
    % Set the box of the item to equal all vois included
    SL.LOC(1).box=[SL.LOC(:).voi]; 
    SL.analysis.loop=1;
    SL.LOC(2:end)=[];
end
%=========================================================================%
% ROI RSA Calc
%=========================================================================%
if (SL.region.noSL==1)
    % Assume only one model being run
    if ~exist(fullfile(SL.dir.outpath,'ROI',subject),'dir'),
        mkdir(fullfile(SL.dir.outpath,'ROI',subject));
    end    
end
%=========================================================================%
%% Volvume Calculations
%=========================================================================%
c=1;
% Approx 36.5K loops
INC_events=size(SL.files,2);
INC_events_min=INC_events/16;

L=SL.analysis.loop(1:int32(length(SL.LOC)*SL.analysis.sparse_sampling));    
for curLOC = L
    R=[];
    %==========
    % 1st level correlation
    %==========
    % Pull values from maps in given region;
    tmp = SL.files(SL.LOC(curLOC).box,:);    
    tmp(isnan(mean(tmp,2)),:)=[]; % Remove NaN voxels
    
    N2=size(tmp,1);  % Voxel w/ removed
   
    % Continue loop if too few values found
    % out{5}.mat(SL.LOC(curLOC).voi,:)=N2;
    if N2<SL.analysis.voxel_per, continue; end
    
    % Else calcule the primary correlation
    R=corr(tmp,'type','Pearson'); % Comparison design
    
    if sum(isnan(R(:,1)))>INC_events_min, continue; end
    % out{6}.mat(SL.LOC(curLOC).voi,:)=sum(isnan(R(:,1)));
    
    if SL.fisher==1, R=0.5*log((1+R)./(1-R)); end
    
    %==========
    % 2nd level correlation & bootstrap
    %==========
    if SL.region.noModel==0
        for ii=1:length(SL.design.matrix)
            out=RSA_bootstrap(subject,R,SL.LOC(curLOC).voi,out,ii,out_name,tmp);
        end
    end

    % If searchlight analysis - save key parameters
    if (SL.region.noSL==1)
        design=SL.design;
        [~,mname,~]=fileparts(SL.region.mask);
        if isempty(SL.dir.QA)
            ROI_save_file=fullfile(SL.dir.outpath,'ROI',subject,'data',[SL.dir.QA mname,'.mat']);
        else
            ROI_save_file=fullfile(SL.dir.outpath,'ROI',subject,'data',[SL.dir.QA '_' mname,'.mat']);  
        end
        
        if (exist(ROI_save_file,'file') && SL.dir.overwrite==0)
            fprintf('.....ROI data saved, likely run, skipping\n');
            return; 
        end
        
        if ~exist(fullfile(SL.dir.outpath,'ROI',subject,'data'),'dir'),
            mkdir(fullfile(SL.dir.outpath,'ROI',subject,'data'));
        end
        save(ROI_save_file,'R','design','tmp','out');
        if SL.region.noModel==1, return; end;
    end
    
    c=c+1;
    if rem(c,1000)==0,
        fprintf(['|---Voxels: ' num2str(c) ' @ ' datestr(now) '---|\n']);
    end
end   
%=========================================================================%
%% Searchlight Data
%=========================================================================%
% Save DataMap
if SL.region.use_mask==0, 
    SL.files=[]; 
    QA_dir=fullfile(SL.dir.outpath,SL.dir.QA,filesep);
    QA_subj=fullfile(QA_dir,SL.dir.subjects{cursub});
    if ~exist(QA_subj,'dir'), mkdir(QA_subj); end
    save(fullfile(QA_subj,'SL.mat'),'SL');
end

if (SL.region.noSL==0)
    fprintf(strcat('....Saving for\t',subject,'\n'));
    if ~exist(fullfile(SL.dir.outpath,subject),'dir'), 
        mkdir(fullfile(SL.dir.outpath,subject));
    end
    
    for ii=1:length(out)
        data=out{ii}.mat;
        SL.V.fname = fullfile(SL.dir.outpath,subject,[out{ii}.name '.img']);
        SL.V.descrip = out{ii}.name;
        if (~exist(SL.V.fname,'file') || SL.dir.overwrite==1)       
            if size(out{ii}.mat,2)==1
                % 3D file
                data = reshape(data,[SL.V.dim(1) SL.V.dim(2) SL.V.dim(3)]);
                SL.V.dim=[SL.V.dim(1) SL.V.dim(2) SL.V.dim(3)];
                spm_write_vol(SL.V,data);
            else
                % 4D file
                data = reshape(data,[SL.V.dim(1) SL.V.dim(2) SL.V.dim(3) size(out{ii}.mat,2)]);
                SL.V.dim=[SL.V.dim(1) SL.V.dim(2) SL.V.dim(3) size(out{ii}.mat,2)];
                % spm_write_vol(SL.V,data); Only for 3D volumes
                get_nii=load_nii(SL.design.ID_file{1}); % Get template
                get_nii.img=data;
                get_nii.fileprefix=fullfile(SL.dir.outpath,subject,[out{ii}.name]);
                get_nii.hdr.dime.dim(1:5)=[4 SL.V.dim];
                
                design.Box=SL.design.Box;
                design.ID_file=SL.design.ID_file;
                design.ID_descrip=SL.design.ID_descrip;
                design.ID_idx=SL.design.ID_idx;
                
                save_nii(get_nii,SL.V.fname);
                save(fullfile(SL.dir.outpath,subject,[out{ii}.name '.mat']),'design');
            end
        else
            display([SL.V.fname ': Already Exists - skipping']);
        end
  
  
    end
    
end

end % Function end
%=========================================================================%
%% Debug notes
%=========================================================================%
% Exist in an attempt to reconstruct the data matrix R from prior data sets
% that were run in a different fashion
% bug.prior_dir='X:\Erik\ERMatch_X\Analysis\STwa\SL_vols_5vox';
% bug.analysis={'Viv1_ID' 'Viv1_RAND' ...
%     'Viv2_ID' 'Viv2_RAND' ...
%     'Viv3_ID' 'Viv3_RAND' ...
%     'Viv4_ID' 'Viv4_RAND'};  
% bug.file_string={'ID*' 'RAND*'...
%     'ID*' 'RAND*'...
%     'ID*' 'RAND*'...
%     'ID*' 'RAND*'};
% bug.on=1;
function bug_test(SL,R,subject,curLOC,bug)

% Load old image and reshape - being a taad risky here in that we
% assume that all images are in the same normalized reference frame
% - witch theoritcally they should be...
c=1;
for kk=1:length(bug.analysis)
    file_list=dir(fullfile(bug.prior_dir,bug.analysis{kk},...
        [subject '_' bug.analysis{kk}(1:4) '_' bug.file_string{kk} '*.img']));
    for f=1:length(file_list)
        [~,fn,~]=fileparts(file_list(f).name);
        ID_idx(c)=str2double(fn(end-1:end));
        if isnan(ID_idx(c));
            ID_idx(c)=str2double(fn(end));
        end
        file_load{c}=fullfile(bug.prior_dir,bug.analysis{kk},file_list(f).name);
        file_type{c}=bug.analysis{kk}(findstr('_',bug.analysis{kk})+1:end);
        c=c+1;

    end
end

% Load all the prior values
for kk=1:length(file_load)
    prior_val(kk)=nan; 
    if ~isnan(ID_idx(kk))
        img_file=fullfile(file_load{kk});  
        img_data=spm_read_vols(spm_vol(img_file));
        img_line=reshape(img_data,[prod(SL.V.dim) 1]); % in 1D
        prior_val(kk)=img_line(SL.LOC(curLOC).voi);
    end
end 

% Replicate calculations for MEs
% 1) Multiply the data DRM by 
RSA_pdm(SL.design.matrix{2},SL.design.Box,'NonMatch');
Match=SL.design.matrix{2};
Match(Match==-1)=NaN;
RSA_pdm(Match,SL.design.Box,'Match');
NonMatch=SL.design.matrix{2};
NonMatch(NonMatch==1)=NaN;
NonMatch(NonMatch==-1)=1;
RSA_pdm(NonMatch,SL.design.Box,'NonMatch');
% Dimension should not matter here
Match_value=nanmean(Match.*R);
NonMatch_value=nanmean(NonMatch.*R);

% Output results
for kk=find(strcmp('ID',file_type))
    if ~isnan(ID_idx(kk))
        ID_val=ID_idx(kk);
        fprintf(['IDX: ' num2str(ID_val) '\n']);
        v_idx=find(SL.design.ID_idx==ID_val);
        zz=1;
            fprintf([' OLD: ' num2str(prior_val(kk)) '\n']);
            display(['  File: ' file_load{kk}]);
            fprintf([' NEW: ' num2str(Match_value(v_idx(zz))) '\n']);
            display(['  File: ' SL.design.ID_file{v_idx(zz)}]);
            display(['  Dscp: ' SL.design.ID_descrip{v_idx(zz)}]);
            display(['  Post: ' num2str(SL.design.ID_idx(v_idx(zz)))]);
            fprintf([' DIF: ' num2str(Match_value(v_idx(zz))-prior_val(kk)) '\n'])
    end
end
for kk=find(strcmp('RAND',file_type))
   if ~isnan(ID_idx(kk))
        ID_val=ID_idx(kk);
        fprintf(['IDX: ' num2str(ID_val) '\n']);
        v_idx=find(SL.design.ID_idx==ID_val);
        zz=1;
            fprintf([' OLD: ' num2str(prior_val(kk)) '\n']);
            display(['  File: ' file_load{kk}]);
            fprintf([' NEW: ' num2str(NonMatch_value(v_idx(zz))) '\n']);
            display(['  File: ' SL.design.ID_file{v_idx(zz)}]);
            display(['  Dscp: ' SL.design.ID_descrip{v_idx(zz)}]);
            display(['  Post: ' num2str(SL.design.ID_idx(v_idx(zz)))]);
            fprintf([' DIF: ' num2str(NonMatch_value(v_idx(zz))-prior_val(kk)) '\n'])
    end
end
        
end
