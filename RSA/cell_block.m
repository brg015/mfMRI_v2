function [r,c]=cell_block(x,y,Box)

if x==1
    r=1:Box(1);
else
    r=sum(Box(1:x-1))+1:sum(Box(1:x));
end

if y==1
   c=1:Box(1);
else
    c=sum(Box(1:y-1))+1:sum(Box(1:y));
end

