% BR Geib Summer 2015
%
% Description: 
%
% Inputs
%   R1 (N x N - numeric) => Network 1
%   R2 (N x N - numeric) => Network 2
%   npermpute (numeric)  => number of iterations
% Outputs
%   FSR_emp (1 x N - numeric) => emperical FSR values
%   FSR_stt (1 x N - numeric) => FSR statistical values
%   FSR_z   (1 x N - numeric) => z-transformed emperical values
%
% Notes:
%   1) code could clearly be optimized in many sections
function [FSR_emp,FSR_z]=GT_FSR_calc(R1,R2,npermute)

% Define basic parameters
N=size(R1,1); % Number of nodes

%=========================================================================%
% Compute theoritical distribution
%=========================================================================%
% c=1; % intialize counter
% for ii=1:npermute
%     I1=randperm(N); % Random assignment R1
%     I2=randperm(N); % Random assignment R2
%     R1_tmp=R1(I1,I1); % Scramble R1
%     R2_tmp=R2(I2,I2); % Scramble R2
%     for jj=1:N
%         FSR_thr(c)=1-corr(R1_tmp(:,jj),R2_tmp(:,jj)); % FSR calc
%         c=c+1;                                        % iterate
%     end % node pair loop
%     clear I1 I2 R1_tmp R2_tmp; % clear loop vars
% end % permutation loop
% L=length(FSR_thr);
% FSR_thr (1 x X) -> theory distribution
% L (X)           -> number of values computed
%
% Observations: theory distribution looks fairly uniform, as could be
% expected based upon the randomness involved
%=========================================================================%
% Compute emperical distribution
%=========================================================================%
c=1; % reset counter
for jj=1:N
    FSR_emp(c)=1-corr(R1(:,jj),R2(:,jj)); % FSR calc
    c=c+1;                                % iterate
end % node pair loop
FSR_z=zscore(FSR_emp);
% FSR_emp (1 x N) -> emperical distrubiton
%
% Observations: emperical distribution isn't quite gaussian, z-transform is
% a bit suspicious here
%=========================================================================%
% Compute significance
%=========================================================================%
% [FSR_thr_sort,~]=sort(FSR_thr); % Sort thr values
% for jj=1:N
%     [~,I]=min(abs(FSR_thr_sort-FSR_emp(jj)));
% end


















