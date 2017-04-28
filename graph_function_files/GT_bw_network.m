function R2=GT_bw_network(N,R,adj)

switch adj
    case 'abs', Radj=abs(R);
    otherwise,  Radj=R;
end

for ii=1:size(Radj,1)
    for jj=1:size(Radj,2)
        if (ii~=jj && ii~=N && jj~=N)
            n1=Radj(ii,N);
            n2=Radj(jj,N);
            R2(ii,jj)=(n1+n2)/2;
        else
            R2(ii,jj)=0;
        end
    end
end