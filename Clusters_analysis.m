load ChronicWoundDataBatchCorrected.mat

% Select non-treated samples only, or treated, or all:

%k=find(treatment_label>1);   %   non-treated patients only
%k=find(treatment_label<2);   %   treated patients only
k=find(treatment_label>0);  %   all patients:
p=patient(k);
w=week_num(k);
D=data(:,k);
h=healed_label(k);
clear k


up=unique(p);
for i=1:numel(up)
    k=find(ismember(p,up(i)));
    uh(i,1)=h(k(1));
end
disp([num2str(sum(uh)),' healers vs ',...
    num2str(numel(uh)-sum(uh)),' non-healers'])
clear k i
%% Find gene clusters
KeyGenes={'IL1B','KRT6B','KRT10'};
%KeyGenes={'CXCL8','DSP','FLG2'};
%KeyGenes={'SLC2A3','KRT17','KRT1'};
%KeyGenes={'PLEK','KRT6A','KRT9'};

for i=1:3
    k=find(ismember(Gene,KeyGenes(i)));
    for j=1:numel(Gene)
        tmp=corrcoef(data(k,:),data(j,:));
        c(j)=tmp(1,2);
    end
    [M,N]=maxk(c,20);
    Clust(:,i)=Gene(N);
    ClustExp=D(N,:);
    ClustVal(i,:)=mean(ClustExp);
end
clear M N i j k c
save('ClusterGenes.mat','Clust')
% Scatterplots
ClusterNames={'Inflammatory Cluster','Proliferative Cluster','Proliferative Cluster 2'};
pairs=[2,3;2,1;3,1];
 ps = cellfun(@(s) sscanf(s,'BAART %d'), p);

f1=figure;
for i=1:3
    k1=pairs(i,1);
    k2=pairs(i,2);
    subplot(1,3,i)
    plot(ClustVal(k1,h==0),ClustVal(k2,h==0),'o','MarkerSize',3,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    hold on
    plot(ClustVal(k1,h==1),ClustVal(k2,h==1),'d','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor','#3070B6')
    %plot(ClustVal(k1,h==1),ClustVal(k2,h==1),'d','MarkerSize',3,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    
    % Add Patient number to the plot:
    %text(ClustVal(k1,h==1),ClustVal(k2,h==1),num2str(ps(h==1)))
    %text(ClustVal(k1,h==0),ClustVal(k2,h==0),num2str(ps(h==0)))
    
    xlabel(ClusterNames(k1))
    ylabel(ClusterNames(k2))
    axis square;
end
%legend({'non-healer','healer'})
exportgraphics(f1,'picture1.pdf','ContentType','vector');
%%
bnd1=-0.1; % boundary between proliferative and other regions
bnd2=0.2; % boundary between inflammation and impairment
f7=figure;
i=2;
    k1=pairs(i,1);
    k2=pairs(i,2);
    fill([bnd1 2 2 bnd1],[-4 -4 2 2],[0.9686    0.9961    0.6863],'EdgeColor','none')
    hold on
    fill([-3 bnd1 bnd1 -3],[-4 -4 bnd2 bnd2],[0.9647    0.7765    0.9529],'EdgeColor','none')
    plot([bnd1 bnd1],[-4 bnd2],'k')
    plot([-3 bnd1],[bnd2 bnd2],'--','Color',0.4*[1 1 1])
    plot(ClustVal(k1,h==0),ClustVal(k2,h==0),'o','MarkerSize',6,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    plot(ClustVal(k1,h==1),ClustVal(k2,h==1),'d','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor','#3070B6')    
    
    % plot non-treated samples only
    %plot(ClustVal(k1,h==0&treatment_label>1),ClustVal(k2,h==0&treatment_label>1),'o','MarkerSize',6,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    %plot(ClustVal(k1,h==1&treatment_label>1),ClustVal(k2,h==1&treatment_label>1),'d','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor','#3070B6')    
    
    % plot treated samples only
    %plot(ClustVal(k1,h==0&treatment_label<2),ClustVal(k2,h==0&treatment_label<2),'o','MarkerSize',6,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    %plot(ClustVal(k1,h==1&treatment_label<2),ClustVal(k2,h==1&treatment_label<2),'d','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor','#3070B6')    
    
    xlabel(ClusterNames(k1))
    ylabel(ClusterNames(k2))
    ax = gca; 
    ax.FontSize = 16;
    axis square;

%legend({'non-healer','healer'})
exportgraphics(f7,'Scatterplot_with_regions.pdf','ContentType','vector');
%% Time-series for one cluster

% Select cluster to plot (ctp=1, 2, or 3)
for ctp=1:3
%ctp=1;
f2(ctp)=figure;
filename2=append('norm_cluster',num2str(ctp),'_timeseries.pdf');
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
    MCV=ClustVal(ctp,k);
    %re-arrange,time must be inceasing
    [tmp,tmp1]=mink(time,numel(time));
    plot(time(tmp1),X(:,tmp1),'.-','Color','#AAAAAA')
    hold on
    plot(time(tmp1), MCV(tmp1),'Color','#000000')
    title(up(i),'FontSize',5)
    set(gca,'fontsize', 5) 
    xlim([min(w),max(w)])
end
exportgraphics(f2(ctp),filename2,'ContentType','vector');
end

%% Time-series for 3 clusters
fclust=figure;
for ctp=1:3   
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
    MCV=ClustVal(ctp,k);
    %re-arrange,time must be inceasing
    [tmp,tmp1]=mink(time,numel(time));
    plot(time(tmp1), MCV(tmp1),'.-')
    hold on
    title(up(i),'FontSize',5)
    set(gca,'fontsize', 5) 
    xlim([min(w),max(w)])
end
end
%subplot(9,5,41),legend({'Infl','Prolif','Prolif2'},'Location','eastoutside')
exportgraphics(fclust,'picture2.pdf','ContentType','vector');



%% Scatterplots 3D
ClusterNames={'Inflammatory Cluster','Proliferative Cluster','Proliferative Cluster 2'};
%pairs=[2,3;2,1;3,1];
 ps = cellfun(@(s) sscanf(s,'BAART %d'), p);

f9=figure;
%for i=1:3
    %k1=pairs(i,1);
    %k2=pairs(i,2);
    %subplot(1,3,i)
    plot3(ClustVal(3,h==0),ClustVal(2,h==0),ClustVal(1,h==0),'o','MarkerSize',3,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    hold on
    plot3(ClustVal(3,h==1),ClustVal(2,h==1),ClustVal(1,h==1),'d','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor','#3070B6')
    %plot(ClustVal(k1,h==1),ClustVal(k2,h==1),'d','MarkerSize',3,'MarkerEdgeColor','#000000','LineWidth',0.2,'MarkerFaceColor','#DD3B22')
    
    % Add Patient number to the plot:
    %text(ClustVal(k1,h==1),ClustVal(k2,h==1),num2str(ps(h==1)))
    %text(ClustVal(k1,h==0),ClustVal(k2,h==0),num2str(ps(h==0)))
    
    xlabel(ClusterNames(3))
    ylabel(ClusterNames(2))
    zlabel(ClusterNames(1))
    axis square;
%end
%legend({'non-healer','healer'})
exportgraphics(f9,'picture9.pdf','ContentType','vector');
