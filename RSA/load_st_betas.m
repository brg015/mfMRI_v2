% BRG modified 2014 (Spring)
% Updated to limit loading of data based upon which files are needed.
function [betadata] = load_st_betas(SL)
    bugger=0;
    V = spm_vol(SL.design.ID_file{1});
    nvox = V.dim(1)*V.dim(2)*V.dim(3);
    betadata = zeros(nvox,sum(SL.design.Box));

%     fprintf('%-40s: %40s','Reading images','...working')

    for i = 1:size(SL.design.ID_file,1)
        dprint(['Trial ' n2sp(i,3) ': ' SL.design.ID_descrip{i}],bugger);
        for j= 1:size(SL.design.ID_file,2)
            dprint([' Image ' num2str(j) ': ' SL.design.ID_file{i,j}],bugger);
            % Read from Vbeta index (aligns w/ POS_idx)
            tmp = spm_read_vols(spm_vol(SL.design.ID_file{i,j}));
            % Assign to POS_idx
            betadata(:,i,j) = reshape(tmp,[nvox 1]);    
        end
    end
%     fprintf('%s%40s\n',repmat(sprintf('\b'),1,40),'...done')

% Now lets consider weights if preprocessing data is being utilized
if SL.preprocess.on==1
   for ii=1:length(SL.preprocess.advance.w)
        betadata(:,:,ii)=betadata(:,:,ii).*SL.preprocess.advance.w(ii);
   end
   betadata=mean(betadata,3);
  
end