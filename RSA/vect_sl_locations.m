
function [LOC] = vect_sl_locations(SPM,SL)

if ~isfield(SL,'region'), SL.region.use_mask=0; end

switch SL.region.use_mask
    case 0
        file_template=SPM.Vbeta(1).fname;
    case 1
        file_template=SL.region.mask;
end

brain_data = spm_read_vols(spm_vol(file_template));
[x,y,z] = size(brain_data);
LOC = struct('voi',cell(1),'box',zeros(SL.design.SSL^3,1));
nn = 1;

switch SL.region.use_mask
    case 0 % Brain
        for zz = 1:z   
            for yy = 1:y    
                for xx = 1:x
                    if ~isnan(brain_data(xx,yy,zz))                
                        [combos] = get_combos(xx,yy,zz,SL.design.SSL,x,y,z);
                        LOC(nn).voi = ((x*y*(zz-1)) + (x*(yy-1)) + xx);                
                        box = zeros(length(combos),1);
                        for ccombo = 1:length(combos)
                            box(ccombo,1) = ((x*y*(combos(ccombo,3)-1)) + (x*(combos(ccombo,2)-1)) + combos(ccombo,1));
                        end
                        LOC(nn).box = box;
                        nn = nn + 1;
                    end
                end
            end
        end
    case 1 % Mask
        for zz = 1:z   
            for yy = 1:y    
                for xx = 1:x
                    if brain_data(xx,yy,zz)>=1
                        [combos] = get_combos(xx,yy,zz,SL.design.SSL,x,y,z);
                        LOC(nn).voi = ((x*y*(zz-1)) + (x*(yy-1)) + xx);                
                        box = zeros(length(combos),1);
                        for ccombo = 1:length(combos)
                            box(ccombo,1) = ((x*y*(combos(ccombo,3)-1)) + (x*(combos(ccombo,2)-1)) + combos(ccombo,1));
                        end
                        LOC(nn).box = box;
                        nn = nn + 1;
                    end
                end
            end
        end
end

function [locs] = get_combos(xx,yy,zz,SSL,x,y,z)

plusminus = floor(SSL/2);
locs = zeros(SSL^3,3);
n = 1;
for i = (xx-plusminus):(xx+plusminus)
    for j = (yy-plusminus):(yy+plusminus)
        for k = (zz-plusminus):(zz+plusminus)
            locs(n,:) = [i j k];
            n = n + 1;   
        end
    end
end


ind = zeros(length(locs),4);
ind(:,1) = sum((locs<1),2);
ind(:,2) = locs(:,1)>x;
ind(:,3) = locs(:,2)>y;
ind(:,4) = locs(:,3)>z;

dropind = logical(sum(ind,2));
locs(dropind,:) = [];

