function GT_write_maps(A)
global SL;
global ntwk;

N=length(SL.analyses);
for ii=1:length(ntwk.contrast.nam)
    SL.analyses{end+1}=ntwk.contrast.nam{ii};
end
X=load(fullfile(SL.dir.save,'IncludedROIs'));
X.MINC=logical(X.MINC);

for ii=1:length(ntwk.measures)
    data=eval(['A.' ntwk.measures{ii}]); 
    sdisp(ntwk.measures{ii},2);
    
    for jj=1:N
        c=1;
        for kk=1:length(SL.dir.subjects)  
            if SL.dir.include(kk)==1
                temp=data{jj,c}; data{jj,c}=[];
                data{jj,c}(X.MINC(:,c))=temp;
                data{jj,c}(~X.MINC(:,c))=NaN;
                sdata{jj,c}=data{jj,c}';
                c=c+1;
            end
        end
    end
    clear data; data=sdata;
    
    for jj=1:length(SL.analyses)
        save_asA=fullfile(SL.dir.save,SL.analyses{jj},ntwk.measures{ii},'avg.nii');
        save_asB=fullfile(SL.dir.save,SL.analyses{jj},ntwk.measures{ii},'T.nii');
        save_asC=fullfile(SL.dir.save,SL.analyses{jj},ntwk.measures{ii},'shift_avg.nii');
        save_asD=fullfile(SL.dir.save,SL.analyses{jj},ntwk.measures{ii},'shift_T.nii');
           
        if ~exist(fullfile(SL.dir.save,SL.analyses{jj}),'dir'), 
            mkdir(fullfile(SL.dir.save,SL.analyses{jj}));
        end
        
        if jj>N
            PC=strcmp(ntwk.contrast.pos{jj-N},SL.analyses);
            NC=strcmp(ntwk.contrast.neg{jj-N},SL.analyses);
            c=1;
            for kk=1:length(SL.dir.subjects),
                if SL.dir.include(kk)==1
                    data{jj,c}=data{PC,c}-data{NC,c};
                    dataZ{jj,c}=nan_zscore(data{PC,c})-nan_zscore(data{NC,c});
                    c=c+1;
                end
            end
        else
            c=1;
            for kk=1:length(SL.dir.subjects)
                if SL.dir.include(kk)==1
                    dataZ{jj,c}=nan_zscore(data{jj,c});
                    c=c+1;
                end
            end
        end
        
        try
        [~]=ROI_constructor(SL.region.mask,data,jj,'avg',save_asA,1);
        [~]=ROI_constructor(SL.region.mask,data,jj,'T',save_asB,1);
        
        [~]=ROI_constructor(SL.region.mask,dataZ,(jj),'avg',save_asC,1);
        [~]=ROI_constructor(SL.region.mask,dataZ,(jj),'T',save_asD,1);
        catch err
            keyboard;
        end
    end
end
