function p = run_permutation(group1Flist, group2Flist, sameSubjectFlag, folderOut, hasNaN,name, numPerm)
% 
%   This is a wrapper for two functions, the first being Create_JC_Network,
%   which creates a similarity matrix for the two groups with the
%   Jaccardized Chekanowski index.
% 
%   The second code computes the permutation test using the 
% 
%   Parameters:
%          group1Flist - matrix 1 of subjects and their nodes
%          group2Flist - matrix 2 of subjects and their nodes
%          sameSubjectFlag - do the two groups represent the same people
%               1 - (YES) enusre that Group 1 Subject 1 is the same person
%                   as Group 2 Subject 1
%               0 - NO
%          numPerm - how many permutations (Default 300000)
%          folderOut - Where to put the output files
%          hasNaN = 1 clears out of NaN data and keep track of nodes with NaN data in a removed_node file
%                   0 skips it
%          name - name of file of its plots/values to be saved
%          modified from xnet_permutation, WFU code         


% sets to default 30000
if nargin < 9
    numPerm = 300000;
end


% runs function that clears out NaN in both matrix
if hasNaN ==1
    group1Flist = clear_NaN(group1Flist,folderOut);
    group2Flist = clear_NaN(group2Flist,folderOut);
end
% removes ID of subject, uncomment
%  group1Flist(:,1) =[];
%  group2Flist(:,1) =[];



%makes matrix
fJacc = Create_JC_matrix(group1Flist, group2Flist, folderOut, 0);

numGroup1 = size(group1Flist,1);
numGroup2 = size(group2Flist,1);
%runs permutation function
p = permutationTest(fJacc, numGroup1, numGroup2, sameSubjectFlag, name,numPerm, 1)