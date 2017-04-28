% Only works w/ SPM8
%
% Inputs:
%	img_cor   => cell array of images
%	contrast  => array of manipulations e.g. [1 1 -1] means to add the first
%					 two images and subtract the third
%	save_name => where to save the file adn what to name
%  varargin  => typically this can be set to 1, which is the defualt unless
%					 a value is passed in. All this allows for is the division
%					 of the images values which can be helpful for normalization
%					 and such
function GT_contrast(img_cor,contrast,save_name,varargin)

[out_dir,out_fil,suf]=fileparts(save_name);

for ii=1:length(img_cor)
	matlabbatch{1}.spm.util.imcalc.input{ii} = [img_cor{ii} ',1'];
end

matlabbatch{1}.spm.util.imcalc.output =[out_fil suf];
matlabbatch{1}.spm.util.imcalc.outdir = {[out_dir '/']};


if ~isempty(varargin), 
	switch varargin{1}
		case 'divide' 
			sdz=varargin{2};
			ivector=['(' (contrast2ivector(contrast)) ') / ' num2str(sdz)];
		case 'function'
			ivector=[varargin{2} '(' contrast2ivector(contrast) ')'];
        case 'divide_image'
            ivector='i1 ./ i2';
	end
else
	ivector=contrast2ivector(contrast);
end
			
matlabbatch{1}.spm.util.imcalc.expression = ivector;

matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
save([out_dir '/' out_fil '_jobfile.mat'],'matlabbatch');
spm_jobman('run',matlabbatch);

end

