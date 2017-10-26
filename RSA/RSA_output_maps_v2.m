function [out,out_name]=RSA_output_maps_v2(subject)
global SL;

out={};
out_name={};
N_design=length(SL.design.matrix);
% Setup "Key" contrast. These represent the main value that is output by a
% given operation
for ii=1:N_design
    M=length(out);
    switch SL.design.calc{ii}
        case 'Anova1'
            out{M+1}.mat=nan(prod(SL.V.dim),1);  
            out{M+1}.name=[SL.design.save_str{ii} '_key'];
        case 'Identity1'
            out{M+1}.mat=nan(prod(SL.V.dim),1);
             out{M+1}.name=[SL.design.save_str{ii} '_MMmean'];
             out{M+2}.mat=nan(prod(SL.V.dim),1);
             out{M+2}.name=[SL.design.save_str{ii} '_Mmean'];
             out{M+3}.mat=nan(prod(SL.V.dim),1);
             out{M+3}.name=[SL.design.save_str{ii} '_Zmean'];
             if SL.design.Identity1{ii}.FourD==1
                 out{M+4}.mat=nan(prod(SL.V.dim),sum(nansum(SL.design.matrix{ii})>0));
                 out{M+4}.name=[SL.design.save_str{ii} '_Zseries'];
                 out{M+5}.mat=nan(prod(SL.V.dim),sum(nansum(SL.design.matrix{ii})>0));
                 out{M+5}.name=[SL.design.save_str{ii} '_Iseries'];
             end
        case 'MVPA'
            out{M+1}.mat=nan(prod(SL.V.dim),1);  
            out{M+1}.name=[SL.design.save_str{ii} '_key'];
        case 'distance'
            out{M+1}.mat=nan(prod(SL.V.dim),1);  
            out{M+1}.name=[SL.design.save_str{ii} '_key'];
        case 'MRegression'
            for jj=1:length(SL.design.matrix{ii})
                % Setup output maps for each beta image
                [~,nam,~]=fileparts(SL.design.model{ii}{jj});
                SL.design.MRname{ii}{jj}=[SL.design.save_str{ii} '_' nam '_key'];
                out{M+jj}.mat=nan(prod(SL.V.dim),1);  
                out{M+jj}.name=[SL.design.save_str{ii} '_' nam '_key'];
                V(:,jj)=reshape(SL.design.matrix{ii}{jj},1,[]);
            end
            % Now, let's make linear vectors for each feature while
            % ignoring NaNs
            V=V(~isnan(sum(V')),:);
            if SL.design.ortho(ii)==1
               % Ortho V to the first column 
               % Not yet implemented
            end
            SL.design.MR{ii}=V;
            clear V;
        otherwise % Assumes custom model
            out{M+1}.mat=nan(prod(SL.V.dim),1);  
            out{M+1}.name=[SL.design.save_str{ii} '_key'];
    end
end

for ii=1:length(out), out_name{ii}=out{ii}.name; end
%=========================================================================%
if (exist(fullfile(SL.dir.outpath,subject,[out{ii}.name '.img']),'file') && SL.dir.overwrite==0)
    display(['Found: ' fullfile(SL.dir.outpath,subject,[out{ii}.name '.img'])]);
    display(' Advancing to next subject');
    SL.err=1;
end