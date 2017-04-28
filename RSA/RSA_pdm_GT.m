function RSA_pdm_GT(dm,v,name,save_name)

if ~isa(dm,'double'), dm=double(dm); end

dm(end+1,end+1)=0;

Box=ones(length(dm),1);

h=pcolor(dm); set(gcf,'position',[0 0 1280 1024]);
view([0 270]); 
colorbar; title(name,'Interpreter','none');
set(h,'edgealpha',0); box off; hold on;
for jj=1:length(Box)
    x=ones(sum(Box)+1,1)*(sum(Box(1:jj))+1);
    y=1:sum(Box)+1;
    plot(x,y,'k','LineWidth',1); hold on;
    plot(y,x,'k','LineWidth',1); hold on;
end
% plot(1:(sum(Box)+1),1:(sum(Box)+1),'k','LineWidth',3);
set(gcf,'color','w');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
if v~=0
    if min(min(dm))<0
        axis([1 size(dm,2) 1 size(dm,1) 0 0.5 -v v])
    else
        axis([1 size(dm,2) 1 size(dm,1) 0 0.5 0 v])
    end
end
try set(1,'color',[1 1 1]); end % Random fails occasionally - for safety
if exist('save_name','var')
    export_fig(save_name); %'format1','png',h);
    close(1);
end
