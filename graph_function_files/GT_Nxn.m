function [er_term,NXNp,ii]=GT_Nxn(X1,iter)
% Inputs
%   X1   => Adjacency matrix 1
%   X2   => Adjacency matrix 2
%   iter => Maximum number of iterations
% Outputs
%   er_term => accuracy over time
%   NXNp    => predicted NSR
%   ii      => iterations required
tol=1e-4;     % Preset tolerance
N=length(X1); % Number of nodes
I=logical(eye(N,N));

% 1) Intialize NSRg (guess)
NXNg=X1;
% 2) Calculate NSRp (predicted)
for ii=1:N
    DC(ii)=mean(NXNg(setdiff(1:N,ii))); % -> make this fixed
end

for ii=1:iter
    for jj=1:N
       for kk=1:N
           if jj==kk
               NXNp(jj,kk)=1;
           else
               NXNp(jj,kk)=NXNg(jj,kk)*DC(jj);
           end
       end
    end
    % Fix the prediction based upon the weight...
     for jj=1:N
        NXNp(jj,setdiff(1:N,jj))=(NXNp(jj,setdiff(1:N,jj))./mean(NXNp(jj,setdiff(1:N,jj)))).*DC(jj);
     end
    
    % 3) Determine the error term RSS
    er_term(ii)=sum((reshape(NXNg,1,[])-reshape(NXNp,1,[])).^2);
    % 4) If error term below set tolerance, exit
    if er_term(ii)<tol, return; end;
    % 5) Else, determine a new guess NSRg based on the mean
    NXNg=mean([reshape(NXNg,1,[]);reshape(NXNp,1,[])],1);
    NXNg=reshape(NXNg,N,N);
    % 6) Now we need to renormalize based upon DC
    keyboard;
    for jj=1:N
        NXNg(jj,setdiff(1:N,jj))=(NXNg(jj,setdiff(1:N,jj))./mean(NXNg(jj,setdiff(1:N,jj)))).*DC(jj);
    end
    NXNg(I)=1;
end