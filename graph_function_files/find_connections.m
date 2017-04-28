function [con,v,V]=find_connections(R,I,N)

V=R(I,:);
keep=round(length(V)*N);
[y,ind]=sort(V,2,'descend');
con=ind(1:keep);
v=y(1:keep);
