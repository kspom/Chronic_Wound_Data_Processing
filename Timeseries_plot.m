% Plotting time-series, non-normalized expression for one cluster
clear
load ChronicWoundData.mat
load ClusterGenes.mat

p=patient;
up=unique(p);
h=healed_label;
w=week_num;
D=log2(1+data);

f18=figure;
% Select cluster to plot (ctp=1, 2, or 3)
ctp=3; 
filename1=append('Clust',num2str(ctp),'.pdf');
clear N
for i=1:20
    N(i)=find(ismember(Gene,Clust(i,ctp)));
end
%plot position in 6x5 plot figure
Np = numel(p);        
r = mod(1:Np, 5);
r(r == 0) = 5;
ih = find(r <= 2); ihc = 1;
in = find(r >= 3); inc = 1;

%remove patient with 1 sample point
tmp=find(ismember(up,'BAART 086'));
tmp1=1:numel(up); tmp1(tmp)=[];

%for i=1:numel(up)
for i=tmp1
    k=find(ismember(p,up(i)));
    %select place to plot - in helaer or non-healer part
    if h(k(1))==1
        ii=ih(ihc);
        ihc=ihc+1;
    else
        ii=in(inc);
        inc=inc+1;
    end
    subplot(9,5,ii)
    clear time X
    time=w(k);
    X=D(N,k);
%    MCV=ClustVal(ctp,k);
    %re-arrange,time must be inceasing
    [tmp,tmp1]=mink(time,numel(time));
    plot(time(tmp1),X(:,tmp1),'.-')
    hold on
%    plot(time(tmp1), MCV(tmp1),'Color','#000000')
    title(up(i),'FontSize',5)
    set(gca,'fontsize', 5) 
    xlim([min(w),max(w)])
end
exportgraphics(f18,filename1,'ContentType','vector');
