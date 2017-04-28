function [filePath, fileName, fileExt, fileVer] = wfu_fileparts(inFileName,extFlag,nameFlag)

%   Returns the path, name, and extension of a file 
%   FORMAT [path,name,ext,ver] = wfu_fileparts(file,extFlag,nameFlag) 
%   WFU_FILEPARTS is platform dependent
%
%   path            - directory path
%   name            - filename 
%   ext             - file extension
%
%   extFlag         - if == 1, name includes extension
%   nameFlag        - if == 1, only name is returned
%
%   See also FILEPARTS, FULLFILE, PATHSEP, FILESEP.
%   ***********************************************
%   Edited 3/8/12 to make compatible with MATLAB2011b and later
%_____________________________________________________________

if nargin<1, return; end
if nargin<2 || isempty(extFlag)
    extFlag = false; 
end
if nargin<3 || isempty(nameFlag) 
    nameFlag = false; 
end

if iscell(inFileName)
    inFileName = inFileName{1};
end

if ~isempty(inFileName)
    [filePath, fileName, fileExt] = fileparts(inFileName);
else
    filePath = [];
    fileName = [];
    fileExt = [];
    
end
fileVer = [];

if extFlag
    fileName = [fileName,fileExt]; 
end
if nameFlag
    filePath = fileName;
end
