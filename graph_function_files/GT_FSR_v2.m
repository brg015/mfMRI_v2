function FSR=GT_FSR_v2(R1,R2)

N=length(R1);
for ii=1:N
    A1(ii)=corr(R1(ii,:)',R2(ii,:)');
end
A2=-atanh(A1);
A3=(A2-mean(A2))./std(A2);

FSR=A3;

end


