function [er_term,NSRp,FSR,ev,ii]=GT_NSR(X1,X2,iter)
% Inputs
%   X1   => Adjacency matrix 1
%   X2   => Adjacency matrix 2
%   iter => Maximum number of iterations
% Outputs
%   er_term => accuracy over time
%   NSRp    => predicted NSR
%   FSR     => first step reorganization
%   ev      => evolution
%   ii      => iterations required
tol=1e-4;        % Preset tolerance
N=length(X1); % Number of nodes
%=========================================================================%
% Original
%=========================================================================%
% 1) Intialize NSRg (guess)
NSRg=ones(N,1);
% 2) Calculate NSRp (predicted)
for ii=1:iter
    for jj=1:N
       temp=weightedcorrs([X1(setdiff(1:N,jj),jj),X2(setdiff(1:N,jj),jj)],NSRg(setdiff(1:N,jj)));
       NSRp(jj,1)=1-temp(1,2); clear temp;
    end
    if ii==1, FSR=NSRp; end
    % 3) Determine the error term RSS
    er_term(ii)=sum((NSRg-NSRp).^2);
    % 4) If error term below set tolerance, exit
    % if er_term(ii)<tol, return; end;
    % 5) Else, determine a new guess NSRg based on the mean
    % NSRg=mean([NSRg,NSRp],2);
    % 6) Update is NSRg
    NSRg=NSRp;
    temp=corr([NSRp,FSR],'type','Spearman');
    ev(ii)=temp(1,2); clear temp;
end