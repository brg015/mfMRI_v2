function L=GT_gelphi_sub(H,Cmat,T,GT)
% BRG 2015 (function file)
%
% Description: accepts an adjacency matrix (Cmat) and tresholds it based
% upon GT.pmap based upon either 'p' or 'T' value. If GT.pmap does not
% exist, then Cmat is used and is thresholded based upon variable 'T'.
%
% Input variables
%   H (numeric)          -> identifies the seed node
%   Cmat (square matrix) -> adjacency matrix to reduce, connection
%                           strengths in this matrix are written out
%   T (numeric)          -> threshold to adjust GT.pmap
%   GT (structure - optional)
%       GT.pmap (stat. square matrix) -> if exist then is used as the
%                                        threshold map instead of Cmat
%       GT.sign ('T' | 'p' 'negT')    -> threshold type, if 'T' looks for
%                                        greater than, if 'p' looks for
%                                        less than (wrt to var T)
%       GT.ring (search depth)        -> Max is 4, min is 3 (defualt)
%   GT (structure - nonoptional) * this is essential to writing out .gdf
%   files. It writes out nodal properties based upon .gdf format and the
%   nodes that are retained in the network.
%       GT.node(x).name   -> variable name
%       GT.node(x).val    -> variable values
%       GT.node(x).class  -> variable class
%   * the first variable name must be 'name'. This is the name of the node,
%   nodes should be named Rx where 'x' is the nodal position.
%   * val values must be strings! this is written to a text file
%       GT.save_to        -> fullfile path of .gdf file to save
%
% Output variables 
%   L (cell str)    -> list of ROIs included
% Planned updates.
% 1) allow for multiple nodes to be input at the same time
% 2) allow for writing of edge properties as well

if ~isfield(GT,'pmap'), GT.pmap=Cmat; end
minmap=GT.pmap; 
if ~isfield(GT,'sign'), GT.sign='T'; end
if ~isfield(GT,'ring'), GT.ring=3; end

% Setup a bunch of null matrices
Cshr=zeros(size(Cmat));
tCshr=zeros(size(Cmat));
Cshr1=zeros(size(Cmat));
Cshr2=zeros(size(Cmat));
if GT.ring>2, Cshr3=zeros(size(Cmat)); end
if GT.ring>3, Cshr4=zeros(size(Cmat)); end

% start the seed
Cshr1(H,H)=1;

switch GT.sign
    case 'p', LT=1;
    case 'T', LT=0;
    case 'negT', LT=1;
end

% Contrast the map to find connections
if LT==1, I1=minmap(:,H)<T; else I1=minmap(:,H)>T; end
for aa=find(I1)'
    Cshr2(H,aa)=1; Cshr2(aa,H)=1;
    
    if LT==1, I2=minmap(:,aa)<T; else I2=minmap(:,aa)>T; end
    
    if GT.ring>2
        for bb=find(I2)'
            Cshr3(aa,bb)=1; Cshr3(bb,aa)=1;
                
            if GT.ring>3
                if LT==1, I3=minmap(:,bb)<T; else I3=minmap(bb,:)>T; end
                for cc=find(I3)'
                    Cshr4(bb,cc)=1;
                    Cshr4(cc,bb)=1;
                end
            end
            
        end
    end
end    
%=========================================================================%
% Derive networks from contrast
%=========================================================================%
% Chsr1 => Central Ring
% Cshr2 => Secondary Ring...
if GT.ring>3
    I=Cshr4>0; Cshr(I)=GT.ring-3; temp=zeros(size(Cmat)); n=find(sum(I)>0); 
    for aa=n, for bb=n, if aa~=bb, temp(aa,bb)=1; end; end; end
    Itemp=temp>0; B(3)=nanmean(Cmat(Itemp));
end

if GT.ring>2
    I=Cshr3>0; Cshr(I)=GT.ring-2; temp=zeros(size(Cmat)); n=find(sum(I)>0); 
    for aa=n, for bb=n, if aa~=bb, temp(aa,bb)=1; end; end; end
    Itemp=temp>0; B(2)=nanmean(Cmat(Itemp)); vB2=Cmat(Itemp);
end

I=Cshr2>0; Cshr(I)=GT.ring-1; temp=zeros(size(Cmat)); n=find(sum(I)>0); 
for aa=n, for bb=n, if aa~=bb, temp(aa,bb)=1; end; end; end
Itemp=temp>0; B(1)=nanmean(Cmat(Itemp)); vB1=Cmat(Itemp);

I=Cshr1>0; Cshr(I)=GT.ring; 

I=Cshr>0; tCshr(I)=Cmat(I);

I=sum(Cshr)>0; 
Cshr(~I,:)=[]; Cshr(:,~I)=[];
tCshr(~I,:)=[]; tCshr(:,~I)=[];

[~,~,~,stats]=ttest2(vB1,vB2);

% Get degree distribution:
A=sum(Cshr>0);

% figure(1);
% subplot(2,2,1); RSA_pdm_GT(tCshr,0,''); title('T-values');
% subplot(2,2,2); RSA_pdm_GT(Cshr,0,'');  title('Ring');
% subplot(2,2,3); 
% bar(A); xlabel('Region'); ylabel('Degree');
% subplot(2,2,4);
% bar(B); xlabel('Degree'); ylabel('Mean Connectivity(t)');
% title(['T = ' num2str(stats.tstat)]);

save([GT.save_to(1:end-4) '_net.mat'],'tCshr','Cshr','I');

% Figure out included nodes
RingVal=max(Cshr);
LabelI=find(strcmp({GT.node.name},'label'));
if ~isempty(LabelI)
   LabelValues=GT.node(LabelI).val;
end
L.labels=LabelValues(I);
L.ring=RingVal;

%=========================================================================%
% Save gelphi network
%=========================================================================%
% Define the ring parameter as well
GT.include=I;

GTn=length(GT.node);
GT.node(GTn+1).name='ring';
GT.node(GTn+1).class='DOUBLE';
GT.node(GTn+1).val=zeros(1,length(Cshr));
GT.node(GTn+1).val(I)=max(Cshr);

write_gelphi_v2(GT,tCshr)





