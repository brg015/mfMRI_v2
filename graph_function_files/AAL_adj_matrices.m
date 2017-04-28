function R=AAL_adj_matrices(flist,fsave)
%=========================================================================%
% Assumptions
%=========================================================================%
% 1) Matrices are sorted 1:90 AAL in alphabetical order and the cerebellum
% is not included
for ii=1:length(flist)
    X=load(flist{ii});
    R(:,:,ii)=X.R;
end
R=mean(R,3);

% Rlabels and Region_code have a 1:1 mapping - there are 45 values here,
% one for each region LR is considered next
Rlabels={'Med. Temp.' 'Subcortical' 'Occipital',...
    'Frontal' 'Temporal' 'Parietal'};
Region_code=[1 6 3 2 4 ...
    6 6 3 4 4 ...
    4 4 4 4 4 ...
    4 4 3 5 1 ...
    5 3 3 3 3 ...
    2 2 1 6 6 ...
    6 6 6 6 2 ...
    4 5 6 6 5 ...
    5 1 1 5 2];
c=1;
for ii=1:2:90
    Rc(ii)=Region_code(c);
    Rc(ii+1)=Region_code(c)+6;
    Rn{ii}=['L. ' Rlabels{Region_code(c)}];
    Rn{ii+1}=['R. ' Rlabels{Region_code(c)}];
    c=c+1;
end
[~,b]=sort(Rc);

% Potential box outline
for ii=1:12
    Box(ii)=sum(Rc==ii);
end
RSA_pdm(R,Box,'',fsave)

% % adj. matrices
% labels=data{4}.col(b);
% R=Adj.N1m(b,b); % Hit adj
% RSA_pdm(R,Box,'')
% 
% R=Adj.N2m(b,b);
% RSA_pdm(R,Box,'')
% 
% 
% 
% h=imagesc(R,[0 1]);
% colormap('jet');
% ax=gca;
% ax.XTickLabel=labels;
% ax.XTickLabelRotation=90;
% ax.XTick=1:length(R);
% ax.YTickLabel=labels;
% ax.YTickLabelRotation=0;
% ax.YTick=1:length(R);
% ax.FontName='Times New Roman';
% ax.FontSize=8;
% ax.GridLineStyle='none';
% ax.TickLength=[0 0];
% colorbar;
% set(gcf,'color','white')





