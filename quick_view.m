% Function file
% quick_view
% B.R. Geib
%
% Purpose:
%   Save screenshots from xjview with premade settings
%    
% Inputs: quick_view(type,dir_file,dir_save)
%   type      => file name e.g. con_00001.img
%   dir_file  => root directory
%   dir_save  => where to save (doesn't need suffix e.g. .jpg)
%
% Uses functions
%   xjview_mod.m => Changed to allow preset p and cluster values
%
function quick_view(map_file,dir_save,Ylim)

% Save image
if exist(map_file,'file')==2
   display(['\tFile Found: ' map_file ': \n']);
else
   display(['\tFile: ' map_file ' DNE\n']);
   return
end

xjview(map_file); 

% Modify intensity outside of xjview
if exist('Ylim','var')
    hObject=gcf;
    handles = guidata(hObject);
    Ylow=find(handles.currentDisplayIntensity{1}<Ylim(1));
    handles.currentDisplayIntensity{1}(Ylow)=Ylim(1);
    Yhigh=find(handles.currentDisplayIntensity{1}>Ylim(2));
    handles.currentDisplayIntensity{1}(Yhigh)=Ylim(2);

    [handles.hReg, handles.hSection, handles.hcolorbar] = ...
        Draw(handles.currentDisplayMNI, handles.currentDisplayIntensity, hObject, handles);
end

% Save Figure
if exist('dir_save','var')
    if ~isempty(dir_save)
        if ~exist(dir_save,'file')
            fg = spm_figure('FindWin','Graphics');
            % Where to save screenshots
            win_print(['Saved ' dir_save '.jpg\n']);
            print(fg, '-noui', '-djpeg', dir_save);
        end  
    end
end


