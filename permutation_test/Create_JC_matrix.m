function matrix = Create_JC_matrix(epimem, semmem,folderOut,thresh)
% create similarity matrix from two different matrices
% code is based on Create_JC_Network based from WFU
%   Detailed explanation goes here
if nargin< 4
    thresh = 0;
end

if nargin < 3
     error('Need at least 3 inputs!');
end


mkdir(folderOut);

allmem = [epimem;semmem];

numberSubjects = size(allmem,1);



Origname = [folderOut, '/MI_all_JC.mat'];
k=1;

for i = 1:numberSubjects
%     fprintf('\t%d of %d\n', i, numberSubjects)
    data1 = allmem(i,:);
    
   
    data1 = data1 .* (data1 > thresh);
    for j = 1:i
%         fprintf('\t\t%d of %d\n', j, i);
        if j==i
            data2 = data1;
        else
            data2 = allmem(j,:);
          data2 = data2 .* (data2 > thresh);
        end
        % Calculate Czekanowski Index

        yy = (double(data1(:)) - double(data2(:)));
        x = sum(abs((double(data1(:)) - double(data2(:)))));
        y = sum(data1(:) + data2(:));
        temp = 1 - (sum(abs((double(data1(:)) - double(data2(:)))))/sum(data1(:) + data2(:)));
       % Jaccardize it
        temp = temp/(2- temp);
        MI_all(i,j) = temp;
        MI_all(j,i) = MI_all(i,j);
        
    end
end
JI_all = MI_all;
JI_orig = JI_all;
save(Origname,'JI_all','JI_orig');
matrix = Origname;