function mask_value=ROI_constructor(mask_list,data,contrast,fxn,save_as,form)
% Inputs
%  mask_list => list of included mask
%  data      => data{ROI}[subject X contrast]
%  contrast  => contrast to plot
%  fxn       => how to summarize across subjects
%  save_as   => file output
%  form      => format of data

% Configure data if needed - really just reshapes it.
if exist('form','var'),
    switch form
        case 1 % data{contrast X subject}(ROI)
            A=length(data{1,1}); % # ROIs
            B=size(data,1);      % # of contrast
            C=size(data,2);      % # of subjects
            for ii=1:A
                for jj=1:B
                    for kk=1:C
                        new_data{ii}(kk,jj)=data{jj,kk}(ii);
                    end
                end
            end
            clear data; data=new_data; clear new_data;
        case 2 % data{contrast}(ROI)
            A=length(data{1});   % # ROIs
            B=length(data);      % # of contrast
            C=1;                 % # of subjects
            for ii=1:A
                for jj=1:B
                    new_data{ii}(1,jj)=data{jj}(ii);
                end
            end
            clear data; data=new_data; clear new_data;
        case 3
            PI=data; 
    end
end

% 1) Setup a template
T=load_nii(mask_list{1}); s=size(T.img);
v=zeros(1,numel(T.img));

if ~strcmp(fxn,'pure')
    N=size(data{1},1); % Number of subjects
end

% 2) Loop through masks
for ii=1:length(mask_list) 
    switch fxn
        case 'avg'
            mask_data=data{ii}(:,contrast);
            mask_value(ii)=mean(mask_data);
        case 'T'
            mask_data=data{ii}(:,contrast);
            v_mne=mean(mask_data);
            v_std=std(mask_data);
            mask_value(ii)=v_mne/(v_std/sqrt(N)); 
            clear v_mne v_std;
        case 'pure'
            mask_value(ii)=PI(ii);
        otherwise
            error('Invalid fxn');
    end
    m=load_nii(mask_list{ii});
    v1=logical(reshape(m.img,1,[]));
    v(v1)=mask_value(ii);   
    clear m mask_data v1;
end

% Save the raw mask values as well
 [rut,fil,~]=fileparts(save_as);
if ~strcmp(fxn,'pure')
    for ii=1:length(mask_list)
        for kk=1:size(data{1}(:,contrast),1)
            out_data{kk}.col(ii)=data{ii}(kk,contrast);
            out_data{kk}.header=['Subject_' num2str(kk)];
        end
    end
    write_struct(out_data,fullfile(rut,[fil '_ind.csv']));
end

csv_save{1}.header='Values';
csv_save{1}.col=mask_value;
write_struct(csv_save,fullfile(rut,[fil '.csv']));

T.img=reshape(v,s);
[rut,nam,~]=fileparts(save_as);

% Update some hdr components
T.fileprefix=fullfile(rut,nam);
T.hdr.dime.glmax=max(mask_value);
T.hdr.dime.glmin=min(mask_value);
T.hdr.dime.datatype=16;
T.hdr.dime.bitpix=32;
T.filetype=1;

save_nii(T,save_as);
