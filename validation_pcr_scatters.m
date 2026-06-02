load ChronicWoundDataBatchCorrected.mat
T=readtable('S11_Table_qPCR_results.xlsx');
SampleName=T{:,1};
for i=1:numel(SampleName)
    n(i)=find(ismember(Metadata_filename,SampleName(i)));
end

genes=T.Properties.VariableNames;
genes(1)=[];
genes=upper(genes);
genexp=T{:,2:10};

GeneList={'IL1B', 'KRT6B', 'KRT10'};
ytexts=[0 0 0 -0.2 0 0 .5 0 0 .5;...
        0 0 0 -0.4 0 -0.4 .4 0 0 .5;...
        0.5 -0.5 0 0 0 0 .5 -0.4 0.4 .7];
xtexts=[-1.2 -1.2 -1.2 0.2 -1.5 -1.2 -0.5 -1.5 -1.2 -0.5;...
        -1.2 -1.2 -1.2 0.2 -1.2 -0.5 -0.8 -1.5 -1.2 -0.5;...
        -0.4 -0.4 -1.2 -1.2 -1.5 -1.2 -0.1 0 0  -0.3];

ytxt=[0 0 0 0 0 0.8 0 0 0 -0.9;...
      0 0 0 0 0 0.9 0 0 0 -0.9;...
      0.5 0 0 0 0 0 0 0 .6 0];
xtxt=[0.4 0.4 0.4 0.4 -2.2 -.4 0.4 0.4 0.4 -0.4;...
      0.4 0.4 0.4 0.4 0.4 -0.4 0.4 0.4 0.4 -0.4;...
      -2.2 0.4 0.4 -2.2 -2.2 -2.2 0.4 0.4 0 -2.2];
pxlim=[-7 7;-5.5 5.5; -8 2];
plotorder=[2 1;3 1;2 3];
obj15=figure;
for j=1:3

    k1=plotorder(j,1);
    k2=plotorder(j,2);
    
    n11=find(ismember(genes,GeneList(k1)));
    n12=find(ismember(genes,GeneList(k2)));
    n21=find(ismember(Gene,GeneList(k1)));
    n22=find(ismember(Gene,GeneList(k2)));
    
  
    subplot(2,3,j)
    plot(data(n21,healed_label==0),data(n22,healed_label==0),'o','MarkerSize',3,'MarkerEdgeColor','#000000','LineWidth',0.1,'MarkerFaceColor','#E5A497')
    hold on
    plot(data(n21,healed_label==1),data(n22,healed_label==1),'d','MarkerSize',3,'MarkerEdgeColor','none','MarkerFaceColor','#A0B6D8')    
    plot(data(n21,n),data(n22,n),'ko','LineWidth',1.5, 'MarkerSize',6)
    text(data(n21,n)+xtexts(j,:),data(n22,n)+ytexts(j,:),SampleName,"FontWeight","bold","FontSize",8)
    axis square
    xlabel(GeneList(k1))
    ylabel(GeneList(k2))
    grid on
    
    subplot(2,3,j+3)
    plot(log(genexp(:,n11)),log(genexp(:,n12)),'ko','LineWidth',1.5, 'MarkerSize',6)
    hold on
    text(log(genexp(:,n11))+xtxt(j,:)',log(genexp(:,n12))+ytxt(j,:)',SampleName,"FontWeight","bold","FontSize",8)
    axis square
    xlim(pxlim(k1,:))
    ylim(pxlim(k2,:))
    xlabel(GeneList(k1))
    ylabel(GeneList(k2))
    grid on
end
exportgraphics(obj15,'qpcr_val_scatters.pdf','ContentType','vector')
