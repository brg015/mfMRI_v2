function bool = wild_strfind(str,patt)
% This function is a working progress and is likely to change.
% Currently, it can accept normal string finds, AND (*) finds, or OR (|)
% finds.  Functionality will be added as requested.
%
% I.E.:
%       'Cnd1*Sac1' - returns true for a string with both attributes
%       'Cnd1|Cnd2' - returns true for either condition
%       'Cnd1'      - returns true for strings containing condition


if nargin ~= 2
    error('mind your input arguments - you need two'); 
end

if ~isempty(strfind(patt,'|'))
    type = 'or';
else 
    type = 'and';
end



switch type
    case 'and'
        
        loc = strfind(patt,'*');
        nflags = length(loc)+1;
        flag = [];

        if ~isempty(loc)

            loc = [0 loc length(patt)+1];

            for i = 2:length(loc)
                if isempty(strfind(str,patt(loc(i-1)+1:loc(i)-1)))
                    break
                else
                    flag(i-1) = strfind(str,patt(loc(i-1)+1:loc(i)-1));
                end
            end

        else

            flag = strfind(str,patt);

        end

        if length(flag) == nflags
            bool = 1;
        else
            bool = 0;
        end
        
        
    case 'or'
        
        loc = strfind(patt,'|');
        nflags = length(loc) + 1;
        flag = [];
        
        if ~isempty(loc)

            loc = [0 loc length(patt)+1];

            for i = 2:length(loc)
                if ~isempty(strfind(str,patt(loc(i-1)+1:loc(i)-1)))
                    flag(i-1) = strfind(str,patt(loc(i-1)+1:loc(i)-1));
                end
            end

        else
            flag = strfind(str,patt);
        end
        
        if ~isempty(flag)
            bool = 1;
        else
            bool = 0;
        end
        
        
end