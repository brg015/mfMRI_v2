function [out,noise_out,out_name]=RSA_output_maps_hemi(subject)
global SL;

out={};
out_name={};
noise_out={};
N_design=length(SL.design.matrix);

% Setup "Key" contrast. These represent the main value that is output by a
% given operation
if ~strcmp(SL.Hemi_set.model,'Memory')
    for ii=1:N_design
        out{ii}.mat=nan(prod(SL.V.dim),1);  
        out{ii}.name=[SL.design.save_str{ii} '_key'];
    end
    out{3}.mat=nan(prod(SL.V.dim),1);
    out{3}.name=[SL.design.save_str{2} '_offkey'];
end

% Setup ANOVA-like contrast
if isfield(SL.design,'anova')
    for ii=1:length(SL.design.anova)
        if ~isempty(SL.design.anova{ii})
           M=length(out); % Set current open model
           if strcmp(SL.design.anova{ii}.type,'Identity1')
                 out{M+1}.mat=nan(prod(SL.V.dim),1);
                 out{M+1}.name=[SL.design.save_str{ii} '_MMmean'];
%                  out{M+2}.mat=nan(prod(SL.V.dim),1);
%                  out{M+2}.name=[SL.design.save_str{ii} '_MMstd'];
%                  out{M+3}.mat=nan(prod(SL.V.dim),1);
%                  out{M+3}.name=[SL.design.save_str{ii} '_MMn'];
                 out{M+2}.mat=nan(prod(SL.V.dim),1);
                 out{M+2}.name=[SL.design.save_str{ii} '_Mmean'];
%                  out{M+5}.mat=nan(prod(SL.V.dim),1);
%                  out{M+5}.name=[SL.design.save_str{ii} '_Mstd'];
%                  out{M+6}.mat=nan(prod(SL.V.dim),1);
%                  out{M+6}.name=[SL.design.save_str{ii} '_Mn'];
                 %========================================================%
                 % Added outputs (12/9/2014)
                 %========================================================%
                 % Z id seres
                 out{M+3}.mat=nan(prod(SL.V.dim),1);
                 out{M+3}.name=[SL.design.save_str{ii} '_Zmean'];
%                  out{M+8}.mat=nan(prod(SL.V.dim),1);
%                  out{M+8}.name=[SL.design.save_str{ii} '_Zstd'];
%                  out{M+9}.mat=nan(prod(SL.V.dim),1);
%                  out{M+9}.name=[SL.design.save_str{ii} '_Zn'];
                 % Column means
%                  out{M+10}.mat=nan(prod(SL.V.dim),1);
%                  out{M+10}.name=[SL.design.save_str{ii} '_Cmean'];
%                  out{M+11}.mat=nan(prod(SL.V.dim),1);
%                  out{M+11}.name=[SL.design.save_str{ii} '_Cstd'];
%                  out{M+12}.mat=nan(prod(SL.V.dim),1);
%                  out{M+12}.name=[SL.design.save_str{ii} '_Cn'];
                 % 4D maps as well
                 
%                  out{M+13}.mat=nan(prod(SL.V.dim),sum(nansum(SL.design.matrix{ii})>0));
%                  out{M+13}.name=[SL.design.save_str{ii} '_Cseries'];
%                  out{M+14}.mat=nan(prod(SL.V.dim),sum(nansum(SL.design.matrix{ii})>0));
%                  out{M+14}.name=[SL.design.save_str{ii} '_Zseries'];
%                  out{M+15}.mat=nan(prod(SL.V.dim),sum(nansum(SL.design.matrix{ii})>0));
%                  out{M+15}.name=[SL.design.save_str{ii} '_zZseries'];
                 if sum(sum(SL.design.anova{ii}.f{1}))==0
                     out{M+4}.mat=nan(prod(SL.V.dim),1);
                 else
                     out{M+4}.mat=nan(prod(SL.V.dim),sum(sum(SL.design.anova{ii}.f{1})));
                 end
                 out{M+4}.name=[SL.design.save_str{ii} '_Iseries'];
                 
                 % Beta estimates - right now only one regression model
                 % exists
                 if SL.design.anova{ii}.regress.on==1
                     cc=1;
                     if SL.design.anova{ii}.regress.E.uni==1
                         out{M+16+cc}.mat=nan(prod(SL.V.dim),1);
                         out{M+16+cc}.name=[SL.design.save_str{ii} '_zBeta_Enc_Uni'];
                         cc=cc+1;
                     end
                     if SL.design.anova{ii}.regress.R.uni==1
                         out{M+16+cc}.mat=nan(prod(SL.V.dim),1);
                         out{M+16+cc}.name=[SL.design.save_str{ii} '_zBeta_Ret_Uni'];
                         cc=cc+1;
                     end
                     if SL.design.anova{ii}.regress.ER.uni==1
                         out{M+16+cc}.mat=nan(prod(SL.V.dim),1);
                         out{M+16+cc}.name=[SL.design.save_str{ii} '_zBeta_ER_Uni'];
                         cc=cc+1;
                     end
                     if SL.design.anova{ii}.regress.Cseries==1
                         out{M+16+cc}.mat=nan(prod(SL.V.dim),1);
                         out{M+16+cc}.name=[SL.design.save_str{ii} '_zBeta_zSET'];
                         cc=cc+1;
                     end
                     if SL.design.anova{ii}.regress.Zseries==1
                         out{M+16+cc}.mat=nan(prod(SL.V.dim),1);
                         out{M+16+cc}.name=[SL.design.save_str{ii} '_zBeta_zID'];
                         cc=cc+1;
                     end
                     A=find(~SL.design.anova{ii}.regress.measures_prd);
                     for kk=1:length(A)
                         out{M+16+cc}.mat=nan(prod(SL.V.dim),1);
                         out{M+16+cc}.name=[SL.design.save_str{ii} SL.regress.name{A(kk)}];
                         cc=cc+1;
                     end
                     clear A;
                 end
                 %========================================================%
            elseif (strcmp(SL.design.anova{ii}.type,'Anova1') || strcmp(SL.design.anova{ii}.type,'Anova2'))
                for jj=1:length(SL.design.anova{ii}.names)
                    % Safer to perallocate, but more difficult to code
                    out{M+(jj-1)*4+1}.mat=nan(prod(SL.V.dim),1);
                    out{M+(jj-1)*4+1}.name=[SL.design.save_str{ii} '_' SL.design.anova{ii}.names{jj} '_v'];
                    out{M+(jj-1)*4+2}.mat=nan(prod(SL.V.dim),1);
                    out{M+(jj-1)*4+2}.name=[SL.design.save_str{ii} '_' SL.design.anova{ii}.names{jj} '_sd'];
                    out{M+(jj-1)*4+3}.mat=nan(prod(SL.V.dim),1);
                    out{M+(jj-1)*4+3}.name=[SL.design.save_str{ii} '_' SL.design.anova{ii}.names{jj} '_N'];
                    out{M+(jj-1)*4+4}.mat=nan(prod(SL.V.dim),1);
                    out{M+(jj-1)*4+4}.name=[SL.design.save_str{ii} '_' SL.design.anova{ii}.names{jj} '_T'];
                end
           end
        end
        if isfield(SL.design.anova{ii},'coeffs')
            M=length(out);
            for jj=1:length(SL.design.anova{ii}.coeffs.idx)
                out{M+jj}.mat=nan(prod(SL.V.dim),1);
                out{M+jj}.name=SL.design.anova{ii}.coeffs.names{jj};
            end
        end
    end
end

if isfield(SL.design,'mvpa')
    for ii=1:length(SL.design.mvpa)
        if ~isempty(SL.design.mvpa{ii})
            M=length(out);
            out{M+1}.mat=nan(prod(SL.V.dim),1);
            out{M+1}.name=[SL.design.save_str{ii} '_T'];
        end
    end
end

% Setup univariate contrast as well
for ii=1:length(SL.design.calc)
   if strcmp(SL.design.calc{ii},'Univariate'),
       M=length(out);
       out{M+1}.mat=nan(prod(SL.V.dim),1);
       out{M+1}.name=[SL.design.save_str{ii} '_mean'];
       out{M+2}.mat=nan(prod(SL.V.dim),1);
       out{M+2}.name=[SL.design.save_str{ii} '_std'];
       out{M+3}.mat=nan(prod(SL.V.dim),1);
       out{M+3}.name=[SL.design.save_str{ii} '_N'];
   end    
end

for ii=1:length(out), out_name{ii}=out{ii}.name; end
%=========================================================================%
if (exist(fullfile(SL.dir.outpath,subject,[out{ii}.name '.img']),'file') && SL.dir.overwrite==0)
    display(['Found: ' fullfile(SL.dir.outpath,subject,[out{ii}.name '.img'])]);
    display(' Advancing to next subject');
    SL.err=1;
end