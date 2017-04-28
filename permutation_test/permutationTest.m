function p = permutationTest(fJacc,numGroup1, numGroup2, sameSubjectsFlag,name, numPerm, typeOption)
%    
%   This function takes in a statistical matrix and permutates it
%   numPerm times, calculating a desired test statistic for each
%   permutation. 
%   
%   A figure is created showing a histogram of the test statistic values. 
% 
%   The p-value returned is calculated by the number of test statistic
%   values greater than the first test statistic value divided by numPerm.
% 
%   Parameters:
%            fJacc - name of .mat file containing JI_all from scaled incl.
%            numGroup1 - the number of subjects in group 1
%            numGroup2 - the number of subjects in group 2
%            sameSubjectsFlag - flag for if you are comparing the same
%                               subjects across different tasks or runs
%                               1 - YES
%                               0 - NO (DEFAULT)
%            numPerm - number of permutations 
%            typeOption - 1 = jaccard
%                         2 = k stat
% 
%    Robert Lyday - 4/30/12
%       Updated for distribution 6/20/12

% checks for bad inputs and sets them to default
fprintf(fJacc);
if nargin < 5
    numPerm = 3000;
end
if nargin < 4
    sameSubjectsFlag = 0;
end
if nargin < 1
    error('Need an input matrix');
end
load(fJacc);
maxNum = length(JI_all);
JI_orig = JI_all;
if nargin<3
    fprintf('Setting size of group 1 and group 2 to %d',maxNum/2);
    numGroup1 = maxNum/2;
    numGroup2 = numGroup1;
end

load(fJacc);
counter = 0;
fprintf('\n');
for i = 1:numPerm
%     clc;
%     fprintf('%.2f',i/numPerm);
    if mod(i, 50000)== 0
        fprintf('%.2f percent complete\n', (i/numPerm*100));
    end
    temp1 = JI_all(1:numGroup1, 1:numGroup1);
    temp2 = JI_all(numGroup1+1:maxNum, numGroup1+1:maxNum);
    a = [temp1(temp1>0);temp2(temp2>0)];
    within = mean([temp1(temp1>0);temp2(temp2>0)]);

    
    temp1 = JI_all(numGroup1+1:maxNum, 1:numGroup1);
    temp2 = JI_all(1:numGroup1,numGroup1+1:maxNum);
    between = mean([temp1(temp1>0);temp2(temp2>0)]);%mean([temp1(:);temp2(:)]);
    if typeOption == 1
        stats(i) = within/between;
    else
        stats(i) = between/within;
    end
    if stats(i)>=stats(1)
        counter=counter+1;
    end
    
    % Option 2 - Switch between groups only
    if sameSubjectsFlag
        temp1 = randi(2,numGroup1,1);
        newOrder = zeros(maxNum,1);
        for j = 1:length(temp1)
            if temp1(j) == 2
                newOrder(j) = numGroup1 + j;
                newOrder(numGroup1 + j) = j;
            else
                newOrder(j) = j;
                newOrder(numGroup1 + j) = numGroup1 + j;
            end
        end
    % Option 1 -Completely Random
    else
        newOrder = randsample(maxNum, maxNum);
    end
    
    
    JI_all=JI_all(newOrder,newOrder);

end

stats(1);
counter;
p =(counter-1)/numPerm; % remove one to not count the intitial state

% *********************************
% Uncomment for plot
% *********************************
% [n xout] = hist(stats,100);
% maxN = max(n);
% figure; hold on;
% h =bar(xout,n); 
% i = 1;
% while xout(i) < stats(1) && i < numel(xout)
%     i = i + 1;
% end
% Put a line where the initial value is
% line([xout(i) xout(i)],[0 maxN], 'Color','r');

hold off;
% **********************************
% END PLOT
% **********************************


[pth fn] = wfu_fileparts(fJacc);
imgName = strcat(pth, '\permDistribution_',name,'.jpg');
saveas(gcf, imgName);
imgName = strcat(pth, '\permDistribution_',name,'.fig');
saveas(gcf, imgName);
save(fJacc, 'JI_all', 'p', 'numPerm', 'counter', 'JI_orig');
